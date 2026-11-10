import sys
import gzip
import os

def get_read_id(header):
    # Header format: @READ_ID 1:N:0:... or @READ_ID/1
    header = header.strip()
    if not header.startswith('@'):
        return None
    first_part = header.split()[0]
    if first_part.endswith('/1') or first_part.endswith('/2'):
        return first_part[:-2][1:]
    return first_part[1:]

def open_file(path, mode):
    if path.endswith('.gz'):
        return gzip.open(path, mode + 't')
    return open(path, mode)

def check_and_repair(r1_path, r2_path, out_r1_path, out_r2_path):
    print(f"[FASTQ Repair] Checking paired reads:\n  R1: {r1_path}\n  R2: {r2_path}")
    
    # First pass: check if already in sync
    in_sync = True
    total_pairs = 0
    mismatched_at = -1
    
    with open_file(r1_path, 'r') as f1, open_file(r2_path, 'r') as f2:
        while True:
            h1 = f1.readline()
            seq1 = f1.readline()
            plus1 = f1.readline()
            q1 = f1.readline()
            
            h2 = f2.readline()
            seq2 = f2.readline()
            plus2 = f2.readline()
            q2 = f2.readline()
            
            if not h1 or not h2:
                if h1 or h2:
                    in_sync = False
                    print("[FASTQ Repair] File length mismatch detected (one file ended earlier).")
                break
                
            id1 = get_read_id(h1)
            id2 = get_read_id(h2)
            
            if id1 != id2:
                in_sync = False
                mismatched_at = total_pairs
                print(f"[FASTQ Repair] Header mismatch at pair #{total_pairs+1}:\n  R1 ID: {id1}\n  R2 ID: {id2}")
                break
                
            total_pairs += 1

    if in_sync:
        print(f"[FASTQ Repair] SUCCESS: All {total_pairs:,} read pairs are 100% synchronized!")
        return False # No repair needed

    print(f"[FASTQ Repair] Desynchronization detected at pair #{mismatched_at+1}. Beginning automated pair synchronization...")
    
    # Pass 2: Index R2 headers or build set of shared IDs
    print("[FASTQ Repair] Indexing R2 read identifiers...")
    r2_ids = set()
    with open_file(r2_path, 'r') as f2:
        while True:
            h2 = f2.readline()
            f2.readline()
            f2.readline()
            f2.readline()
            if not h2:
                break
            id2 = get_read_id(h2)
            if id2:
                r2_ids.add(id2)
                
    print(f"[FASTQ Repair] Indexed {len(r2_ids):,} unique IDs from R2. Indexing R1...")
    
    shared_ids = set()
    with open_file(r1_path, 'r') as f1:
        while True:
            h1 = f1.readline()
            f1.readline()
            f1.readline()
            f1.readline()
            if not h1:
                break
            id1 = get_read_id(h1)
            if id1 and id1 in r2_ids:
                shared_ids.add(id1)
                
    print(f"[FASTQ Repair] Found {len(shared_ids):,} matching read pairs shared between R1 and R2.")
    if len(shared_ids) == 0:
        r1_base = os.path.basename(r1_path)
        r2_base = os.path.basename(r2_path)
        for p in (out_r1_path, out_r2_path):
            if p and os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass
        msg = (
            f"[FASTQ Repair] ERROR: R1 ('{r1_base}') and R2 ('{r2_base}') share 0 matching read pair IDs! "
            f"They appear to belong to two different samples. Please upload the matching R1 and R2 FASTQ files for the same sample."
        )
        print(msg, file=sys.stderr)
        print(msg)
        sys.exit(2)
    
    # Filter and write synchronized R1 and R2
    print(f"[FASTQ Repair] Writing synchronized output files...")
    
    r1_records = {}
    with open_file(r1_path, 'r') as f1:
        while True:
            h1 = f1.readline()
            seq1 = f1.readline()
            p1 = f1.readline()
            q1 = f1.readline()
            if not h1:
                break
            id1 = get_read_id(h1)
            if id1 in shared_ids:
                r1_records[id1] = (h1, seq1, p1, q1)
                
    written_pairs = 0
    with open_file(r2_path, 'r') as f2, open_file(out_r1_path, 'w') as out1, open_file(out_r2_path, 'w') as out2:
        while True:
            h2 = f2.readline()
            seq2 = f2.readline()
            p2 = f2.readline()
            q2 = f2.readline()
            if not h2:
                break
            id2 = get_read_id(h2)
            if id2 in shared_ids and id2 in r1_records:
                h1, seq1, p1, q1 = r1_records.pop(id2)
                out1.write(h1 + seq1 + p1 + q1)
                out2.write(h2 + seq2 + p2 + q2)
                written_pairs += 1
                
    print(f"[FASTQ Repair] SUCCESS: Repaired and synchronized {written_pairs:,} paired reads into output files!")
    return True

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print("Usage: python3 repair_fastq_pairs.py <R1.fastq.gz> <R2.fastq.gz> <OUT_R1.fastq.gz> <OUT_R2.fastq.gz>")
        sys.exit(1)
    check_and_repair(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
