#!/usr/bin/env python3
"""
Generate a standalone interactive responsive HTML report from samtools stats, flagstat, and idxstats.
CDFD - MITOCLIN Platform
Zero external CDN dependencies (100% offline-safe standalone SVG + JS).
"""

import sys
import os
import re
import json
import html
from datetime import datetime

def parse_flagstat(flagstat_path):
    if not flagstat_path or not os.path.isfile(flagstat_path):
        return []
    lines = []
    with open(flagstat_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if line:
                lines.append(line)
    return lines

def parse_stats_sn(stats_path):
    if not stats_path or not os.path.isfile(stats_path):
        return {}
    sn = {}
    with open(stats_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("SN\t"):
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    key = parts[1].rstrip(":")
                    val = parts[2]
                    sn[key] = val
    return sn

def parse_idxstats(idxstats_path):
    if not idxstats_path or not os.path.isfile(idxstats_path):
        return []
    rows = []
    with open(idxstats_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                ref, length, mapped, unmapped = parts[0], parts[1], parts[2], parts[3]
                if ref != "*":
                    rows.append({
                        "ref": ref,
                        "length": int(length) if length.isdigit() else length,
                        "mapped": int(mapped) if mapped.isdigit() else mapped,
                        "unmapped": int(unmapped) if unmapped.isdigit() else unmapped
                    })
    return rows

def parse_stats_extended(stats_path):
    data = {
        "is": [],
        "rl": [],
        "gcf": [],
        "gcc": [],
        "ffq": [],
        "lfq": [],
        "id": []
    }
    if not stats_path or not os.path.isfile(stats_path):
        return data

    with open(stats_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("IS\t"):
                parts = line.strip().split("\t")
                if len(parts) >= 3 and parts[1].isdigit():
                    data["is"].append((int(parts[1]), int(parts[2])))
            elif line.startswith("RL\t"):
                parts = line.strip().split("\t")
                if len(parts) >= 3 and parts[1].isdigit():
                    data["rl"].append((int(parts[1]), int(parts[2])))
            elif line.startswith("GCF\t") or line.startswith("GC\t"):
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    try:
                        data["gcf"].append((float(parts[1]), int(parts[2])))
                    except ValueError:
                        pass
            elif line.startswith("GCC\t"):
                parts = line.strip().split("\t")
                if len(parts) >= 6 and parts[1].isdigit():
                    try:
                        data["gcc"].append({
                            "cycle": int(parts[1]),
                            "A": float(parts[2]),
                            "C": float(parts[3]),
                            "G": float(parts[4]),
                            "T": float(parts[5]),
                            "N": float(parts[6]) if len(parts) > 6 else 0.0
                        })
                    except ValueError:
                        pass
            elif line.startswith("FFQ\t"):
                parts = line.strip().split("\t")
                if len(parts) >= 3 and parts[1].isdigit():
                    cycle = int(parts[1])
                    q_counts = [int(x) for x in parts[2:]]
                    total_q = sum(q_counts)
                    if total_q > 0:
                        weighted_q = sum(q * cnt for q, cnt in enumerate(q_counts)) / total_q
                        data["ffq"].append((cycle, round(weighted_q, 2)))
            elif line.startswith("LFQ\t"):
                parts = line.strip().split("\t")
                if len(parts) >= 3 and parts[1].isdigit():
                    cycle = int(parts[1])
                    q_counts = [int(x) for x in parts[2:]]
                    total_q = sum(q_counts)
                    if total_q > 0:
                        weighted_q = sum(q * cnt for q, cnt in enumerate(q_counts)) / total_q
                        data["lfq"].append((cycle, round(weighted_q, 2)))
            elif line.startswith("ID\t"):
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    try:
                        data["id"].append((int(parts[1]), int(parts[2]), int(parts[3])))
                    except ValueError:
                        pass
    return data

# ══════════════════════════════════════════════════════════════
# INTERACTIVE SVG VISUALIZATION GENERATORS
# ══════════════════════════════════════════════════════════════

def make_interactive_donut(segments, chart_id, center_label="TOTAL", total_val=None):
    total = total_val if total_val is not None else sum(s[1] for s in segments)
    if total == 0:
        return '<div style="text-align:center;padding:30px;color:var(--muted)">No data recorded</div>'
    
    r = 56
    cx, cy = 90, 100
    import math
    circ = 2 * math.pi * r
    
    arcs_svg = []
    accum_pct = 0.0
    
    for idx, (label, val, color) in enumerate(segments):
        pct = (val / total) if total > 0 else 0.0
        dash = pct * circ
        offset = -(accum_pct * circ)
        arcs_svg.append(
            f'<circle id="{chart_id}_seg_{idx}" class="donut-slice" '
            f'cx="{cx}" cy="{cy}" r="{r}" fill="none" stroke="{color}" stroke-width="22" '
            f'stroke-dasharray="{dash:.2f} {circ:.2f}" stroke-dashoffset="{offset:.2f}" '
            f'data-chart="{chart_id}" data-idx="{idx}" data-label="{html.escape(label)}" '
            f'data-val="{val}" data-pct="{pct*100.0:.2f}" data-color="{color}" '
            f'style="cursor:pointer;transition:stroke-width 0.2s, stroke-opacity 0.2s"/>'
        )
        accum_pct += pct
        
    legend_html = []
    for idx, (label, val, color) in enumerate(segments):
        pct = (val / total * 100.0) if total > 0 else 0
        legend_html.append(f"""
        <div id="{chart_id}_leg_{idx}" class="donut-legend-item" data-chart="{chart_id}" data-idx="{idx}"
             style="display:flex;align-items:center;justify-content:space-between;margin-bottom:6px;font-size:11.5px;padding:4px 8px;border-radius:6px;cursor:pointer;transition:background 0.15s">
          <div style="display:flex;align-items:center;gap:6px">
            <span style="width:10px;height:10px;border-radius:3px;background:{color};display:inline-block"></span>
            <span style="color:var(--text);font-weight:600">{label}</span>
          </div>
          <div style="font-family:monospace;font-weight:700;color:var(--text)">
            {val:,} <span style="font-weight:400;color:var(--muted);font-size:10px">({pct:.1f}%)</span>
          </div>
        </div>
        """)
        
    legend_str = "".join(legend_html)
    arcs_str = "".join(arcs_svg)
    
    total_disp = f"{total/1e6:.2f}M" if total >= 1e6 else (f"{total/1e3:.1f}k" if total >= 1e4 else f"{total:,}")
    
    return f"""
    <div id="{chart_id}_wrapper" class="interactive-plot-card" data-chart-type="donut" data-chart-id="{chart_id}"
         style="display:flex;align-items:center;justify-content:space-around;flex-wrap:wrap;gap:12px;padding:6px 4px;position:relative">
      <svg id="{chart_id}_svg" width="180" height="200" viewBox="0 0 180 200" style="overflow:visible;flex-shrink:0">
        <g transform="rotate(-90 90 100)">
          <circle cx="{cx}" cy="{cy}" r="{r}" fill="none" stroke="rgba(226,232,240,0.3)" stroke-width="22"/>
          {arcs_str}
        </g>
        <text id="{chart_id}_hub_lbl" x="90" y="94" text-anchor="middle" font-size="9.5" font-weight="700" fill="var(--muted)">{center_label}</text>
        <text id="{chart_id}_hub_val" x="90" y="112" text-anchor="middle" font-size="13" font-weight="800" fill="var(--text)" font-family="monospace">
          {total_disp}
        </text>
        <text id="{chart_id}_hub_pct" x="90" y="126" text-anchor="middle" font-size="9" font-weight="600" fill="var(--teal)" style="display:none">
          100.0%
        </text>
      </svg>
      <div style="flex:1;min-width:180px;padding:4px">
        {legend_str}
      </div>
    </div>
    """

def make_interactive_line_chart(points, chart_id, width=540, height=240, x_label="", y_label="",
                                line_color="#0284c7", v_lines=None, y_max_fixed=None, x_max_fixed=None,
                                unit_x="bp", unit_y="pairs", title="", is_quality=False, is_gc=False):
    pad_left = 60
    pad_right = 20
    pad_top = 25
    pad_bottom = 42
    plot_w = width - pad_left - pad_right
    plot_h = height - pad_top - pad_bottom

    if not points:
        return '<div style="text-align:center;padding:30px;color:var(--muted)">No data points available</div>'

    pts = sorted(points, key=lambda p: p[0])
    x_min = pts[0][0]
    x_max = x_max_fixed if x_max_fixed is not None else pts[-1][0]
    y_min = 0.0
    y_max = y_max_fixed if y_max_fixed is not None else max(p[1] for p in pts)
    if y_max <= 0: y_max = 1.0
    if x_max <= x_min: x_max = x_min + 1.0

    def to_svg(x, y):
        sx = pad_left + ((x - x_min) / (x_max - x_min)) * plot_w
        sy = pad_top + plot_h - ((y - y_min) / (y_max - y_min)) * plot_h
        return sx, sy

    # Quality background zones if requested
    bands_svg = ""
    if is_quality:
        _, q30_y = to_svg(x_min, 30)
        _, q20_y = to_svg(x_min, 20)
        _, q0_y  = to_svg(x_min, 0)
        _, q42_y = to_svg(x_min, 42)
        bands_svg = f"""
        <g class="bg-quality-bands" pointer-events="none">
          <rect x="{pad_left}" y="{q42_y}" width="{plot_w}" height="{max(0, q30_y - q42_y)}" fill="rgba(16,185,129,0.09)"/>
          <rect x="{pad_left}" y="{q30_y}" width="{plot_w}" height="{max(0, q20_y - q30_y)}" fill="rgba(245,158,11,0.08)"/>
          <rect x="{pad_left}" y="{q20_y}" width="{plot_w}" height="{max(0, q0_y - q20_y)}" fill="rgba(239,68,68,0.07)"/>
          <line x1="{pad_left}" y1="{q30_y:.1f}" x2="{width - pad_right}" y2="{q30_y:.1f}" stroke="#10b981" stroke-width="1.2" stroke-dasharray="4,2"/>
          <line x1="{pad_left}" y1="{q20_y:.1f}" x2="{width - pad_right}" y2="{q20_y:.1f}" stroke="#f59e0b" stroke-width="1.2" stroke-dasharray="4,2"/>
          <text x="{width - pad_right - 6}" y="{q30_y - 4:.1f}" text-anchor="end" font-size="8.5" font-weight="700" fill="#10b981">Q30 (≥99.9% Accuracy)</text>
          <text x="{width - pad_right - 6}" y="{q20_y - 4:.1f}" text-anchor="end" font-size="8.5" font-weight="700" fill="#f59e0b">Q20 Threshold</text>
        </g>
        """

    # Horizontal Gridlines
    grid_lines = []
    y_ticks = 4
    for i in range(y_ticks + 1):
        y_val = y_min + i * (y_max - y_min) / y_ticks
        _, gy = to_svg(x_min, y_val)
        y_str = f"{y_val/1e3:.1f}k" if y_val >= 1e4 else (f"{y_val:,.0f}" if y_val >= 10 else f"{y_val:.1f}")
        grid_lines.append(f'<line x1="{pad_left}" y1="{gy:.1f}" x2="{width - pad_right}" y2="{gy:.1f}" stroke="var(--border)" stroke-width="1" stroke-dasharray="3,3"/>')
        grid_lines.append(f'<text x="{pad_left - 8}" y="{gy + 3.5:.1f}" text-anchor="end" font-size="9" fill="var(--muted)" font-family="monospace">{y_str}</text>')

    # Vertical Gridlines / Ticks
    x_ticks = 5
    for i in range(x_ticks + 1):
        x_val = x_min + i * (x_max - x_min) / x_ticks
        gx, _ = to_svg(x_val, y_min)
        x_str = f"{x_val:,.0f}"
        grid_lines.append(f'<text x="{gx:.1f}" y="{height - pad_bottom + 16}" text-anchor="middle" font-size="9" fill="var(--muted)" font-family="monospace">{x_str}</text>')

    path_pts = [to_svg(x, y) for x, y in pts]
    path_d = f"M {path_pts[0][0]:.1f} {path_pts[0][1]:.1f} " + " ".join(f"L {x:.1f} {y:.1f}" for x, y in path_pts[1:])
    fill_d = f"{path_d} L {path_pts[-1][0]:.1f} {pad_top + plot_h:.1f} L {path_pts[0][0]:.1f} {pad_top + plot_h:.1f} Z"

    # Reference Vertical Lines (e.g. Mean = 462.2 bp, or mtDNA GC 44.4%)
    v_svg = []
    if v_lines:
        for vx, vlabel, vcol in v_lines:
            if x_min <= vx <= x_max:
                vsx, _ = to_svg(vx, y_min)
                v_svg.append(f'<line x1="{vsx:.1f}" y1="{pad_top}" x2="{vsx:.1f}" y2="{pad_top + plot_h}" stroke="{vcol}" stroke-width="2" stroke-dasharray="4,3"/>')
                v_svg.append(f'<rect x="{vsx - 38:.1f}" y="{pad_top + 4}" width="76" height="18" rx="4" fill="{vcol}"/>')
                v_svg.append(f'<text x="{vsx:.1f}" y="{pad_top + 16}" text-anchor="middle" font-size="8.5" font-weight="700" fill="#ffffff">{vlabel}</text>')

    grid_str = "".join(grid_lines)
    v_str = "".join(v_svg)

    # Reference marker for mtDNA GC content if applicable
    if is_gc:
        gc_mt_x, _ = to_svg(44.4, y_min)
        v_str += f"""
        <line x1="{gc_mt_x:.1f}" y1="{pad_top}" x2="{gc_mt_x:.1f}" y2="{pad_top + plot_h}" stroke="#0d9488" stroke-width="1.8" stroke-dasharray="4,3"/>
        <rect x="{gc_mt_x - 42:.1f}" y="{pad_top + 4}" width="84" height="18" rx="4" fill="#0d9488"/>
        <text x="{gc_mt_x:.1f}" y="{pad_top + 16}" text-anchor="middle" font-size="8.0" font-weight="700" fill="#ffffff">mtDNA ~44.4%</text>
        """

    # Interactive tracking elements
    crosshair_svg = f"""
    <g id="{chart_id}_crosshair" style="display:none;pointer-events:none">
      <line id="{chart_id}_guide_x" x1="{pad_left}" y1="{pad_top}" x2="{pad_left}" y2="{pad_top + plot_h}" stroke="var(--teal)" stroke-width="1.3" stroke-dasharray="3,3"/>
      <circle id="{chart_id}_guide_dot" cx="{pad_left}" cy="{pad_top + plot_h}" r="5.0" fill="{line_color}" stroke="#ffffff" stroke-width="2"/>
    </g>
    """

    # Embed data points directly into data attribute for fast JS hover search
    data_pts_json = html.escape(json.dumps([[p[0], p[1]] for p in pts]))

    return f"""
    <div id="{chart_id}_wrapper" class="interactive-plot-card" data-chart-type="line" data-chart-id="{chart_id}"
         data-xmin="{x_min}" data-xmax="{x_max}" data-ymin="{y_min}" data-ymax="{y_max}"
         data-pad-left="{pad_left}" data-pad-right="{pad_right}" data-pad-top="{pad_top}" data-pad-bottom="{pad_bottom}"
         data-width="{width}" data-height="{height}" data-unit-x="{unit_x}" data-unit-y="{unit_y}"
         data-label-x="{x_label}" data-label-y="{y_label}" data-pts="{data_pts_json}"
         style="position:relative;width:100%;cursor:crosshair">
      <svg id="{chart_id}_svg" width="100%" height="{height}" viewBox="0 0 {width} {height}" style="overflow:visible;display:block">
        <defs>
          <linearGradient id="{chart_id}_grad" x1="0" y1="0" x2="0" y2="1">
            <stop offset="0%" stop-color="{line_color}" stop-opacity="0.28"/>
            <stop offset="100%" stop-color="{line_color}" stop-opacity="0.01"/>
          </linearGradient>
        </defs>
        {bands_svg}
        {grid_str}
        <path class="plot-fill-path" d="{fill_d}" fill="url(#{chart_id}_grad)"/>
        <path class="plot-line-path" d="{path_d}" fill="none" stroke="{line_color}" stroke-width="2.3" stroke-linejoin="round" stroke-linecap="round"/>
        {v_str}
        {crosshair_svg}
        <text x="{width/2:.1f}" y="{height - 6}" text-anchor="middle" font-size="10.5" font-weight="700" fill="var(--text)">{x_label}</text>
        <text x="14" y="{height/2:.1f}" text-anchor="middle" transform="rotate(-90 14 {height/2:.1f})" font-size="10.5" font-weight="700" fill="var(--text)">{y_label}</text>
      </svg>
    </div>
    """

def make_interactive_acgt_chart(gcc_data, chart_id="plot_acgt", width=540, height=240):
    if not gcc_data:
        return '<div style="text-align:center;padding:30px;color:var(--muted)">No nucleotide cycle data available</div>'

    pad_left = 55
    pad_right = 20
    pad_top = 25
    pad_bottom = 42
    plot_w = width - pad_left - pad_right
    plot_h = height - pad_top - pad_bottom

    cycles = [d["cycle"] for d in gcc_data]
    x_min = min(cycles)
    x_max = max(cycles)
    y_min = 0.0
    y_max = 60.0 # Standard base composition percentage max

    def to_svg(x, y):
        sx = pad_left + ((x - x_min) / (x_max - x_min)) * plot_w
        sy = pad_top + plot_h - ((y - y_min) / (y_max - y_min)) * plot_h
        return sx, sy

    # Gridlines
    grid_lines = []
    for y_val in [0, 15, 30, 45, 60]:
        _, gy = to_svg(x_min, y_val)
        grid_lines.append(f'<line x1="{pad_left}" y1="{gy:.1f}" x2="{width - pad_right}" y2="{gy:.1f}" stroke="var(--border)" stroke-width="1" stroke-dasharray="3,3"/>')
        grid_lines.append(f'<text x="{pad_left - 8}" y="{gy + 3.5:.1f}" text-anchor="end" font-size="9" fill="var(--muted)" font-family="monospace">{y_val}%</text>')

    x_ticks = 5
    for i in range(x_ticks + 1):
        x_val = x_min + i * (x_max - x_min) / x_ticks
        gx, _ = to_svg(x_val, y_min)
        grid_lines.append(f'<text x="{gx:.1f}" y="{height - pad_bottom + 16}" text-anchor="middle" font-size="9" fill="var(--muted)" font-family="monospace">{x_val:,.0f}</text>')

    bases = [
        ("A", "#10b981", "Adenine (A)"),
        ("C", "#0284c7", "Cytosine (C)"),
        ("G", "#f59e0b", "Guanine (G)"),
        ("T", "#ef4444", "Thymine (T)")
    ]

    paths_svg = []
    for base_key, col, b_name in bases:
        pts = [to_svg(d["cycle"], d[base_key]) for d in gcc_data]
        d_str = f"M {pts[0][0]:.1f} {pts[0][1]:.1f} " + " ".join(f"L {x:.1f} {y:.1f}" for x, y in pts[1:])
        paths_svg.append(f'<path id="{chart_id}_path_{base_key}" class="acgt-line acgt-line-{base_key}" d="{d_str}" fill="none" stroke="{col}" stroke-width="2.2" stroke-linecap="round"/>')

    grid_str = "".join(grid_lines)
    paths_str = "".join(paths_svg)

    # Interactive tracking elements (guide line + 4 dots for A, C, G, T)
    crosshair_svg = f"""
    <g id="{chart_id}_crosshair" style="display:none;pointer-events:none">
      <line id="{chart_id}_guide_x" x1="{pad_left}" y1="{pad_top}" x2="{pad_left}" y2="{pad_top + plot_h}" stroke="var(--teal)" stroke-width="1.3" stroke-dasharray="3,3"/>
      <circle id="{chart_id}_dot_A" cx="{pad_left}" cy="{pad_top}" r="4.2" fill="#10b981" stroke="#ffffff" stroke-width="1.5"/>
      <circle id="{chart_id}_dot_C" cx="{pad_left}" cy="{pad_top}" r="4.2" fill="#0284c7" stroke="#ffffff" stroke-width="1.5"/>
      <circle id="{chart_id}_dot_G" cx="{pad_left}" cy="{pad_top}" r="4.2" fill="#f59e0b" stroke="#ffffff" stroke-width="1.5"/>
      <circle id="{chart_id}_dot_T" cx="{pad_left}" cy="{pad_top}" r="4.2" fill="#ef4444" stroke="#ffffff" stroke-width="1.5"/>
    </g>
    """

    data_gcc_json = html.escape(json.dumps(gcc_data))

    return f"""
    <div id="{chart_id}_wrapper" class="interactive-plot-card" data-chart-type="acgt" data-chart-id="{chart_id}"
         data-xmin="{x_min}" data-xmax="{x_max}" data-ymin="{y_min}" data-ymax="{y_max}"
         data-pad-left="{pad_left}" data-pad-right="{pad_right}" data-pad-top="{pad_top}" data-pad-bottom="{pad_bottom}"
         data-width="{width}" data-height="{height}" data-gcc='{data_gcc_json}'
         style="position:relative;width:100%;cursor:crosshair">
      <div style="display:flex;justify-content:center;gap:12px;margin-bottom:8px;font-size:11.5px;font-weight:700">
        <span id="{chart_id}_btn_A" class="acgt-toggle active" data-base="A" onclick="toggleAcgtSeries('{chart_id}', 'A')" style="color:#10b981;cursor:pointer;padding:2px 8px;border-radius:12px;background:rgba(16,185,129,0.12)">
          <span style="display:inline-block;width:8px;height:8px;border-radius:50%;background:#10b981;vertical-align:middle;margin-right:4px"></span>A (%)
        </span>
        <span id="{chart_id}_btn_C" class="acgt-toggle active" data-base="C" onclick="toggleAcgtSeries('{chart_id}', 'C')" style="color:#0284c7;cursor:pointer;padding:2px 8px;border-radius:12px;background:rgba(2,132,199,0.12)">
          <span style="display:inline-block;width:8px;height:8px;border-radius:50%;background:#0284c7;vertical-align:middle;margin-right:4px"></span>C (%)
        </span>
        <span id="{chart_id}_btn_G" class="acgt-toggle active" data-base="G" onclick="toggleAcgtSeries('{chart_id}', 'G')" style="color:#f59e0b;cursor:pointer;padding:2px 8px;border-radius:12px;background:rgba(245,158,11,0.12)">
          <span style="display:inline-block;width:8px;height:8px;border-radius:50%;background:#f59e0b;vertical-align:middle;margin-right:4px"></span>G (%)
        </span>
        <span id="{chart_id}_btn_T" class="acgt-toggle active" data-base="T" onclick="toggleAcgtSeries('{chart_id}', 'T')" style="color:#ef4444;cursor:pointer;padding:2px 8px;border-radius:12px;background:rgba(239,68,68,0.12)">
          <span style="display:inline-block;width:8px;height:8px;border-radius:50%;background:#ef4444;vertical-align:middle;margin-right:4px"></span>T (%)
        </span>
      </div>
      <svg id="{chart_id}_svg" width="100%" height="{height}" viewBox="0 0 {width} {height}" style="overflow:visible;display:block">
        {grid_str}
        {paths_str}
        {crosshair_svg}
        <text x="{width/2:.1f}" y="{height - 6}" text-anchor="middle" font-size="10.5" font-weight="700" fill="var(--text)">Sequencing Cycle Position (bp)</text>
        <text x="14" y="{height/2:.1f}" text-anchor="middle" transform="rotate(-90 14 {height/2:.1f})" font-size="10.5" font-weight="700" fill="var(--text)">Base %</text>
      </svg>
    </div>
    """

def make_interactive_bar_chart(rl_data, chart_id="plot_rl", width=540, height=240,
                               x_label="Read Length (bp)", y_label="Read Count"):
    if not rl_data:
        return '<div style="text-align:center;padding:30px;color:var(--muted)">No read length data available</div>'

    pad_left = 60
    pad_right = 20
    pad_top = 25
    pad_bottom = 42
    plot_w = width - pad_left - pad_right
    plot_h = height - pad_top - pad_bottom

    pts = sorted(rl_data, key=lambda p: p[0])
    x_min = pts[0][0]
    x_max = pts[-1][0]
    y_min = 0
    y_max = max(p[1] for p in pts)
    if y_max <= 0: y_max = 1.0
    if x_max <= x_min: x_max = x_min + 1.0

    def to_svg(x, y):
        sx = pad_left + ((x - x_min) / (x_max - x_min)) * plot_w
        sy = pad_top + plot_h - ((y - y_min) / (y_max - y_min)) * plot_h
        return sx, sy

    # Gridlines
    grid_lines = []
    y_ticks = 4
    for i in range(y_ticks + 1):
        y_val = y_min + i * (y_max - y_min) / y_ticks
        _, gy = to_svg(x_min, y_val)
        y_str = f"{y_val/1e3:.1f}k" if y_val >= 1e4 else f"{y_val:,.0f}"
        grid_lines.append(f'<line x1="{pad_left}" y1="{gy:.1f}" x2="{width - pad_right}" y2="{gy:.1f}" stroke="var(--border)" stroke-width="1" stroke-dasharray="3,3"/>')
        grid_lines.append(f'<text x="{pad_left - 8}" y="{gy + 3.5:.1f}" text-anchor="end" font-size="9" fill="var(--muted)" font-family="monospace">{y_str}</text>')

    x_ticks = 5
    for i in range(x_ticks + 1):
        x_val = x_min + i * (x_max - x_min) / x_ticks
        gx, _ = to_svg(x_val, y_min)
        grid_lines.append(f'<text x="{gx:.1f}" y="{height - pad_bottom + 16}" text-anchor="middle" font-size="9" fill="var(--muted)" font-family="monospace">{x_val:,.0f}</text>')

    bars_svg = []
    bar_w = max(3.5, (plot_w / len(pts)) * 0.75)
    tot_reads = sum(p[1] for p in pts)

    for x, y in pts:
        bx, by = to_svg(x, y)
        bh = max(1.5, (pad_top + plot_h) - by)
        pct = (y / tot_reads * 100.0) if tot_reads > 0 else 0.0
        bars_svg.append(
            f'<rect class="interactive-bar" x="{bx - bar_w/2:.1f}" y="{by:.1f}" width="{bar_w:.1f}" height="{bh:.1f}" '
            f'fill="#0d9488" rx="2" style="cursor:pointer;transition:fill 0.15s, opacity 0.15s" '
            f'data-chart="{chart_id}" data-x="{x}" data-y="{y}" data-pct="{pct:.2f}" data-label="{x_label}"/>'
        )

    grid_str = "".join(grid_lines)
    bars_str = "".join(bars_svg)

    return f"""
    <div id="{chart_id}_wrapper" class="interactive-plot-card" data-chart-type="bar" data-chart-id="{chart_id}"
         data-xmin="{x_min}" data-xmax="{x_max}" data-ymin="{y_min}" data-ymax="{y_max}"
         data-width="{width}" data-height="{height}" style="position:relative;width:100%">
      <svg id="{chart_id}_svg" width="100%" height="{height}" viewBox="0 0 {width} {height}" style="overflow:visible;display:block">
        {grid_str}
        {bars_str}
        <text x="{width/2:.1f}" y="{height - 6}" text-anchor="middle" font-size="10.5" font-weight="700" fill="var(--text)">{x_label}</text>
        <text x="14" y="{height/2:.1f}" text-anchor="middle" transform="rotate(-90 14 {height/2:.1f})" font-size="10.5" font-weight="700" fill="var(--text)">{y_label}</text>
      </svg>
    </div>
    """

def make_interactive_indel_chart(id_data, chart_id="plot_indel", width=540, height=240):
    if not id_data:
        return '<div style="text-align:center;padding:30px;color:var(--muted)">No indel records recorded in alignment</div>'

    pad_left = 60
    pad_right = 20
    pad_top = 25
    pad_bottom = 42
    plot_w = width - pad_left - pad_right
    plot_h = height - pad_top - pad_bottom

    top_indels = sorted(id_data, key=lambda x: x[0])[:12]
    max_count = max(max(x[1], x[2]) for x in top_indels)
    if max_count <= 0: max_count = 1

    grid_lines = []
    y_ticks = 4
    for i in range(y_ticks + 1):
        y_val = i * max_count / y_ticks
        gy = pad_top + plot_h - (i / y_ticks) * plot_h
        y_str = f"{y_val/1e3:.1f}k" if y_val >= 1e4 else f"{y_val:,.0f}"
        grid_lines.append(f'<line x1="{pad_left}" y1="{gy:.1f}" x2="{width - pad_right}" y2="{gy:.1f}" stroke="var(--border)" stroke-width="1" stroke-dasharray="3,3"/>')
        grid_lines.append(f'<text x="{pad_left - 8}" y="{gy + 3.5:.1f}" text-anchor="end" font-size="9" fill="var(--muted)" font-family="monospace">{y_str}</text>')

    bars_svg = []
    slot_w = plot_w / len(top_indels)
    bar_w = slot_w * 0.35

    for idx, (sz, ins, dels) in enumerate(top_indels):
        cx = pad_left + idx * slot_w + slot_w / 2
        
        # Ins bar
        ins_h = max(1.5, (ins / max_count) * plot_h)
        ins_y = pad_top + plot_h - ins_h
        bars_svg.append(
            f'<rect class="interactive-bar indel-bar-ins" x="{cx - bar_w - 1:.1f}" y="{ins_y:.1f}" width="{bar_w:.1f}" height="{ins_h:.1f}" '
            f'fill="#0d9488" rx="2" style="cursor:pointer;transition:fill 0.15s" '
            f'data-chart="{chart_id}" data-type="Insertion" data-sz="{sz}" data-count="{ins}"/>'
        )

        # Del bar
        del_h = max(1.5, (dels / max_count) * plot_h)
        del_y = pad_top + plot_h - del_h
        bars_svg.append(
            f'<rect class="interactive-bar indel-bar-del" x="{cx + 1:.1f}" y="{del_y:.1f}" width="{bar_w:.1f}" height="{del_h:.1f}" '
            f'fill="#ef4444" rx="2" style="cursor:pointer;transition:fill 0.15s" '
            f'data-chart="{chart_id}" data-type="Deletion" data-sz="{sz}" data-count="{dels}"/>'
        )

        # X label
        bars_svg.append(f'<text x="{cx:.1f}" y="{height - pad_bottom + 16}" text-anchor="middle" font-size="9" fill="var(--muted)" font-family="monospace">{sz} bp</text>')

    grid_str = "".join(grid_lines)
    bars_str = "".join(bars_svg)

    return f"""
    <div id="{chart_id}_wrapper" class="interactive-plot-card" data-chart-type="indel" data-chart-id="{chart_id}"
         data-width="{width}" data-height="{height}" style="position:relative;width:100%">
      <div style="display:flex;justify-content:center;gap:16px;margin-bottom:8px;font-size:11.5px;font-weight:700">
        <span class="indel-toggle active" data-type="ins" onclick="toggleIndelSeries('{chart_id}', 'ins')" style="color:#0d9488;cursor:pointer;padding:2px 8px;border-radius:12px;background:rgba(13,148,136,0.12)">
          <span style="display:inline-block;width:10px;height:10px;border-radius:2px;background:#0d9488;vertical-align:middle;margin-right:4px"></span>Insertions
        </span>
        <span class="indel-toggle active" data-type="del" onclick="toggleIndelSeries('{chart_id}', 'del')" style="color:#ef4444;cursor:pointer;padding:2px 8px;border-radius:12px;background:rgba(239,68,68,0.12)">
          <span style="display:inline-block;width:10px;height:10px;border-radius:2px;background:#ef4444;vertical-align:middle;margin-right:4px"></span>Deletions
        </span>
      </div>
      <svg id="{chart_id}_svg" width="100%" height="{height}" viewBox="0 0 {width} {height}" style="overflow:visible;display:block">
        {grid_str}
        {bars_str}
        <text x="{width/2:.1f}" y="{height - 6}" text-anchor="middle" font-size="10.5" font-weight="700" fill="var(--text)">Indel Length (bp)</text>
        <text x="14" y="{height/2:.1f}" text-anchor="middle" transform="rotate(-90 14 {height/2:.1f})" font-size="10.5" font-weight="700" fill="var(--text)">Event Count</text>
      </svg>
    </div>
    """

# ══════════════════════════════════════════════════════════════
def make_interactive_numt_loci_chart(top_loci, total_numt_reads, chart_id="plot_numt_loci", width=540, height=250):
    if not top_loci:
        return '<div style="text-align:center;padding:30px;color:var(--muted)">No nuclear NUMT reads absorbed in this run</div>'

    pad_left = 155
    pad_right = 68
    pad_top = 14
    pad_bottom = 28
    plot_w = width - pad_left - pad_right
    plot_h = height - pad_top - pad_bottom

    loci = top_loci[:10]
    max_reads = max((int(x.get("reads_absorbed", 0)) for x in loci), default=1)
    if max_reads <= 0:
        max_reads = 1

    n_bars = len(loci)
    slot_h = plot_h / max(1, n_bars)
    bar_h = min(16.0, slot_h * 0.68)

    svg_els = []
    # Vertical gridlines
    for i in range(5):
        v = max_reads * i / 4.0
        gx = pad_left + (i / 4.0) * plot_w
        v_str = f"{v/1e3:.1f}k" if v >= 1000 else f"{v:.0f}"
        svg_els.append(f'<line x1="{gx:.1f}" y1="{pad_top}" x2="{gx:.1f}" y2="{pad_top + plot_h}" stroke="var(--border)" stroke-width="1" stroke-dasharray="3,3"/>')
        svg_els.append(f'<text x="{gx:.1f}" y="{height - 8}" text-anchor="middle" font-size="8.5" fill="var(--muted)" font-family="monospace">{v_str}</text>')

    for idx, item in enumerate(loci):
        locus_lbl = str(item.get("locus") or item.get("contig") or f"NUMT_{idx+1}")
        numt_id = str(item.get("numt_id") or "")
        cnt = int(item.get("reads_absorbed", 0))
        pct_numt = (cnt / total_numt_reads * 100.0) if total_numt_reads > 0 else 0.0
        bw = max(2.0, (cnt / max_reads) * plot_w)
        cy = pad_top + idx * slot_h + slot_h / 2.0
        by = cy - bar_h / 2.0
        bar_col = "#dc2626" if idx == 0 else ("#ea580c" if idx < 3 else "#0284c7")

        svg_els.append(
            f'<text x="{pad_left - 6}" y="{cy + 3.2:.1f}" text-anchor="end" font-size="8.8" font-weight="700" fill="var(--text)" font-family="monospace">{html.escape(locus_lbl)}</text>'
        )
        svg_els.append(
            f'<rect class="interactive-bar" x="{pad_left}" y="{by:.1f}" width="{bw:.1f}" height="{bar_h:.1f}" '
            f'fill="{bar_col}" rx="3" style="cursor:pointer;transition:opacity 0.15s" '
            f'data-chart="{chart_id}" data-x="{html.escape(locus_lbl)} ({html.escape(numt_id)})" '
            f'data-y="{cnt}" data-pct="{pct_numt:.2f}" data-label="Nuclear NUMT Locus"/>'
        )
        svg_els.append(
            f'<text x="{pad_left + bw + 5:.1f}" y="{cy + 3.2:.1f}" text-anchor="start" font-size="8.5" font-weight="700" fill="var(--muted)" font-family="monospace">{cnt:,} ({pct_numt:.1f}%)</text>'
        )

    els_str = "".join(svg_els)
    return f"""
    <div id="{chart_id}_wrapper" class="interactive-plot-card" data-chart-type="bar" data-chart-id="{chart_id}"
         data-width="{width}" data-height="{height}" style="position:relative;width:100%">
      <svg id="{chart_id}_svg" width="100%" height="{height}" viewBox="0 0 {width} {height}" style="overflow:visible;display:block">
        {els_str}
      </svg>
    </div>
    """


def load_numt_depletion_json(sample, stats_path, numt_json_path=None):
    candidates = []
    if numt_json_path:
        candidates.append(numt_json_path)
    if stats_path:
        out_root = os.path.dirname(os.path.dirname(os.path.abspath(stats_path)))
        candidates.append(os.path.join(out_root, "03_mapping", f"{sample}_numt_depletion.json"))
    for p in candidates:
        if p and os.path.isfile(p):
            try:
                with open(p, "r", encoding="utf-8") as f:
                    return json.load(f)
            except Exception:
                pass
    return None


# ══════════════════════════════════════════════════════════════
# MAIN REPORT GENERATOR
# ══════════════════════════════════════════════════════════════

def generate_bam_stats_html(sample, stats_path, flagstat_path, idxstats_path, out_html_path, numt_json_path=None):
    sn = parse_stats_sn(stats_path)
    flagstat_lines = parse_flagstat(flagstat_path)
    idxstats_rows = parse_idxstats(idxstats_path)
    ext_data = parse_stats_extended(stats_path)
    numt_data = load_numt_depletion_json(sample, stats_path, numt_json_path)

    total_reads = sn.get("raw total sequences", "0")
    reads_mapped = sn.get("reads mapped", "0")
    reads_unmapped = sn.get("reads unmapped", "0")
    reads_duplicated = sn.get("reads duplicated", "0")
    error_rate = sn.get("error rate", "0")
    avg_length = sn.get("average length", "0")
    avg_quality = sn.get("average quality", "0")
    avg_insert_size = sn.get("insert size average", "0")
    insert_std = sn.get("insert size standard deviation", "0")
    inward_pairs = sn.get("inward oriented pairs", "0")
    outward_pairs = sn.get("outward oriented pairs", "0")
    other_pairs = sn.get("pairs with other orientation", "0")

    # Percentage calculations
    try:
        tot_num = float(total_reads)
        map_num = float(reads_mapped)
        dup_num = float(reads_duplicated)
        unmap_num = float(reads_unmapped)
        pct_mapped = (map_num / tot_num * 100.0) if tot_num > 0 else 0.0
        pct_dup = (dup_num / tot_num * 100.0) if tot_num > 0 else 0.0
    except (ValueError, ZeroDivisionError):
        tot_num = map_num = dup_num = unmap_num = 0
        pct_mapped = pct_dup = 0.0

    is_numt_competitive = bool(numt_data and numt_data.get("alignment_mode") == "numt_aware_competitive")
    numt_tot_sam = int(numt_data.get("total_sam_records", 0)) if is_numt_competitive else int(tot_num)
    numt_tot_mapped = int(numt_data.get("mapped_records", 0)) if is_numt_competitive else int(map_num)
    numt_unmapped = int(numt_data.get("unmapped_records", 0)) if is_numt_competitive else int(unmap_num)
    numt_rcrs_ret = int(numt_data.get("rcrs_retained_reads", 0)) if is_numt_competitive else int(map_num)
    numt_depleted = int(numt_data.get("numt_depleted_reads", 0)) if is_numt_competitive else 0
    numt_chimeric = int(numt_data.get("chimeric_numt_mate_reads", 0)) if is_numt_competitive else 0
    numt_dep_pct = float(numt_data.get("numt_depletion_pct", 0.0)) if is_numt_competitive else 0.0
    numt_loci_hit = int(numt_data.get("unique_numt_loci_hit", 0)) if is_numt_competitive else 0
    numt_threads = int(numt_data.get("threads_used", 4)) if numt_data else 4
    top_numt_loci = (numt_data.get("top_numt_loci") or []) if is_numt_competitive else []

    if is_numt_competitive and numt_tot_sam > 0:
        overall_map_pct = (numt_tot_mapped / numt_tot_sam) * 100.0
    else:
        overall_map_pct = pct_mapped

    try:
        err_f = float(error_rate) * 100.0
        err_str = f"{err_f:.3f}%"
    except ValueError:
        err_str = error_rate

    # Generate Donut 1: Alignment, NUMT Depletion & Duplication
    non_dup_mapped = max(0, int(map_num - dup_num))
    if is_numt_competitive and numt_tot_sam > 0:
        donut1_segments = [
            ("rCRS Unique Mapped", non_dup_mapped, "#10b981"),
            ("rCRS Duplicated", int(dup_num), "#f59e0b"),
            ("Nuclear NUMT Depleted", numt_depleted, "#dc2626"),
            ("Unmapped Reads", numt_unmapped, "#64748b"),
        ]
        svg_donut_alignment = make_interactive_donut(donut1_segments, chart_id="donut_align", center_label="TOTAL SAM", total_val=numt_tot_sam)
    else:
        donut1_segments = [
            ("Mapped Unique", non_dup_mapped, "#10b981"),
            ("Duplicated", int(dup_num), "#f59e0b"),
            ("Unmapped", int(unmap_num), "#ef4444")
        ]
        svg_donut_alignment = make_interactive_donut(donut1_segments, chart_id="donut_align", center_label="TOTAL READS", total_val=int(tot_num))

    # Generate Donut 2: Pair Orientation
    try:
        in_p = int(inward_pairs)
        out_p = int(outward_pairs)
        oth_p = int(other_pairs)
        tot_p = in_p + out_p + oth_p
    except ValueError:
        in_p = out_p = oth_p = tot_p = 0

    donut2_segments = [
        ("Inward Pairs (FR)", in_p, "#0284c7"),
        ("Outward Pairs (RF)", out_p, "#8b5cf6"),
        ("Other / Chimeric", oth_p, "#64748b")
    ]
    svg_donut_orientation = make_interactive_donut(donut2_segments, chart_id="donut_pairs", center_label="PAIRED MATES", total_val=tot_p)

    # Generate Insert Size Plot
    try:
        mean_is = float(avg_insert_size)
    except ValueError:
        mean_is = 0.0
    vlines = [(mean_is, f"Mean: {mean_is:.1f} bp", "#0284c7")] if mean_is > 0 else None
    svg_plot_is = make_interactive_line_chart(
        ext_data["is"][:350] if ext_data["is"] else [],
        chart_id="plot_is",
        width=520, height=230,
        x_label="Insert Size (bp)", y_label="Read Pairs",
        line_color="#0284c7",
        v_lines=vlines,
        unit_x="bp", unit_y="pairs"
    )

    # Generate ACGT Plot
    svg_plot_acgt = make_interactive_acgt_chart(ext_data["gcc"], chart_id="plot_acgt", width=520, height=230)

    # Generate Quality Score Plot
    svg_plot_quality = make_interactive_line_chart(
        ext_data["ffq"],
        chart_id="plot_quality",
        width=520, height=230,
        x_label="Sequencing Cycle (bp)", y_label="Phred Score (Q)",
        line_color="#0284c7",
        y_max_fixed=42.0,
        unit_x="bp", unit_y="Q",
        is_quality=True
    )

    # Generate Read Length Plot
    svg_plot_rl = make_interactive_bar_chart(
        ext_data["rl"],
        chart_id="plot_rl",
        width=520, height=230,
        x_label="Read Length (bp)", y_label="Read Count"
    )

    # Generate GC Content Plot
    svg_plot_gc = make_interactive_line_chart(
        ext_data["gcf"],
        chart_id="plot_gc",
        width=520, height=230,
        x_label="GC Content (%)", y_label="Read Count",
        line_color="#10b981",
        x_max_fixed=100.0,
        unit_x="%", unit_y="reads",
        is_gc=True
    )

    # Generate Indel Plot
    svg_plot_indel = make_interactive_indel_chart(ext_data["id"], chart_id="plot_indel", width=520, height=230)

    # Generate Top 10 Nuclear NUMT Loci Plot (when NUMT competitive alignment is active)
    if is_numt_competitive and top_numt_loci:
        svg_plot_numt_loci = make_interactive_numt_loci_chart(top_numt_loci, numt_depleted, chart_id="plot_numt_loci", width=520, height=250)
        numt_plot_card_html = f"""
      <!-- PLOT 9: TOP 10 NUCLEAR NUMT LOCI READ DEPLETION -->
      <div class="section" style="grid-column: 1 / -1">
        <div class="section-header">
          <span>Competitive NUMT Read Depletion — Top 10 Absorbing Human Nuclear NUMT Loci (UCSC / RHNumtS.2)</span>
          <div class="plot-toolbar">
            <span style="background:rgba(220,38,38,0.14);color:#dc2626;padding:2px 8px;border-radius:12px;font-weight:700;font-size:11px">{numt_depleted:,} NUMT Reads Depleted ({numt_dep_pct:.2f}% &bull; {numt_loci_hit} Loci)</span>
            <button class="plot-btn" onclick="exportPlotPng('plot_numt_loci_svg', '{sample}_numt_loci_depletion.png')" title="Download PNG">PNG</button>
            <button class="plot-btn" onclick="exportPlotSvg('plot_numt_loci_svg', '{sample}_numt_loci_depletion.svg')" title="Download SVG">SVG</button>
          </div>
        </div>
        <div class="section-body">
          <div style="font-size:11.5px;color:var(--muted);margin-bottom:10px;padding:8px 12px;background:var(--card-header-bg);border-left:3.5px solid #0d9488;border-radius:4px">
            <b>Competitive Alignment Summary:</b> Out of <b>{numt_tot_sam:,}</b> total BWA-MEM SAM records (<b>{numt_tot_mapped:,}</b> mapped), <b>{numt_rcrs_ret:,}</b> authentic mitochondrial reads aligned preferentially to <code>NC_012920.1</code> (rCRS) and were retained in the final deduplicated BAM, while <b>{numt_depleted:,} ({numt_dep_pct:.2f}%)</b> nuclear pseudogene reads across <b>{numt_loci_hit}</b> nuclear NUMT contigs (including <b>{numt_chimeric:,}</b> chimeric rCRS&ndash;NUMT mate pairs) were depleted prior to GATK Mutect2 variant calling.
          </div>
          {svg_plot_numt_loci}
        </div>
      </div>"""
    else:
        numt_plot_card_html = ""

    # Safe float / int parsing helper
    def get_int(k, default=0):
        try: return int(sn.get(k, default))
        except: return default
    def get_flt(k, default=0.0):
        try: return float(sn.get(k, default))
        except: return default

    tot_seq = get_int("raw total sequences") or get_int("sequences")
    mapped = get_int("reads mapped")
    unmapped = get_int("reads unmapped")
    dup = get_int("reads duplicated")
    unique_mapped = max(0, mapped - dup)
    supp = get_int("supplementary alignments")
    mq0 = get_int("reads MQ0")

    paired = get_int("reads paired")
    r1 = get_int("1st fragments")
    r2 = get_int("last fragments")
    prop_paired = get_int("reads properly paired")
    inward = get_int("inward oriented pairs")
    outward = get_int("outward oriented pairs")
    other_pairs = get_int("pairs with other orientation")
    
    avg_is = get_flt("insert size average")
    std_is = get_flt("insert size standard deviation")

    tot_len = get_int("total length")
    bases_mapped = get_int("bases mapped")
    bases_dup = get_int("bases duplicated")
    avg_len = get_flt("average length")
    max_len = get_int("maximum length")
    avg_q = get_flt("average quality")
    mismatches = get_int("mismatches")
    err_rate = get_flt("error rate")

    tot_num = float(tot_seq) if tot_seq > 0 else 1.0
    tot_pairs = float(inward + outward + other_pairs) if (inward + outward + other_pairs) > 0 else 1.0

    # Table 1: Primary Alignment, NUMT Depletion & Duplication Summary
    if is_numt_competitive and numt_tot_sam > 0:
        sam_denom = float(numt_tot_sam)
        map_denom = float(numt_tot_mapped) if numt_tot_mapped > 0 else 1.0
        align_tr_html = f"""
    <tr>
      <td style="font-weight:700;color:var(--text)">Total BWA-MEM SAM Records</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{numt_tot_sam:,}</td>
      <td style="font-family:monospace;color:var(--muted)">100.0%</td>
      <td><span class="badge-green">Competitive rCRS + 766 NUMTs</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Total Mapped Records (rCRS + NUMTs)</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{numt_tot_mapped:,}</td>
      <td style="font-family:monospace;color:#10b981;font-weight:700">{(numt_tot_mapped/sam_denom*100.0):.2f}%</td>
      <td><span class="badge-green">Pre-Depletion Mapped</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Authentic rCRS Retained Reads (NC_012920.1)</td>
      <td style="font-family:monospace;font-weight:700;color:#10b981">{tot_seq:,}</td>
      <td style="font-family:monospace;color:#10b981;font-weight:700">{(tot_seq/map_denom*100.0):.2f}% <span style="font-size:10px;color:var(--muted)">(of mapped)</span></td>
      <td><span class="badge-green">100.0% of Final BAM</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Nuclear NUMT Reads Depleted ({numt_loci_hit} Loci)</td>
      <td style="font-family:monospace;font-weight:700;color:#dc2626">{numt_depleted:,} <span style="font-size:10px;color:var(--muted)">({numt_chimeric:,} chimeric)</span></td>
      <td style="font-family:monospace;color:#dc2626;font-weight:700">{numt_dep_pct:.2f}% <span style="font-size:10px;color:var(--muted)">(of mapped)</span></td>
      <td><span style="background:rgba(220,38,38,0.14);color:#dc2626;padding:2px 8px;border-radius:12px;font-weight:700;font-size:11px">Filtered Pre-Mutect2</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Unique Informative rCRS Reads</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{unique_mapped:,}</td>
      <td style="font-family:monospace;color:var(--teal);font-weight:700">{(unique_mapped/tot_num*100.0):.2f}% <span style="font-size:10px;color:var(--muted)">(of rCRS)</span></td>
      <td><span class="badge-blue">Primary Variant Evidence</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Duplicated rCRS Reads (PCR / Optical)</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{dup:,}</td>
      <td style="font-family:monospace;color:#f59e0b;font-weight:700">{(dup/tot_num*100.0):.2f}% <span style="font-size:10px;color:var(--muted)">(of rCRS)</span></td>
      <td><span style="font-size:11px;color:#d97706;font-weight:600">Flagged in GATK Mutect2</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Unmapped Reads (Filtered)</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{numt_unmapped:,}</td>
      <td style="font-family:monospace;color:#64748b;font-weight:700">{(numt_unmapped/sam_denom*100.0):.2f}% <span style="font-size:10px;color:var(--muted)">(of SAM)</span></td>
      <td><span style="font-size:11px;color:var(--muted);font-weight:600">Non-mtDNA / Exogenous</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Supplementary / Low-MQ0 (in rCRS BAM)</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{supp:,} supp / {mq0:,} MQ0</td>
      <td style="font-family:monospace;color:var(--muted)">{(supp/tot_num*100.0):.2f}% / {(mq0/tot_num*100.0):.4f}%</td>
      <td><span class="badge-green">High Mapping Specificity</span></td>
    </tr>
    """
    else:
        align_tr_html = f"""
    <tr>
      <td style="font-weight:700;color:var(--text)">Total Primary Sequences</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{tot_seq:,}</td>
      <td style="font-family:monospace;color:var(--muted)">100.0%</td>
      <td><span class="badge-green">&ge; 1.0M Optimal</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Primary Mapped Reads</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{mapped:,}</td>
      <td style="font-family:monospace;color:#10b981;font-weight:700">{(mapped/tot_num*100.0):.2f}%</td>
      <td><span class="badge-green">&ge; 80.0% High Quality</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Unique Informative Reads</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{unique_mapped:,}</td>
      <td style="font-family:monospace;color:var(--teal);font-weight:700">{(unique_mapped/tot_num*100.0):.2f}%</td>
      <td><span class="badge-blue">Primary Variant Evidence</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Duplicated Reads (PCR / Optical)</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{dup:,}</td>
      <td style="font-family:monospace;color:#f59e0b;font-weight:700">{(dup/tot_num*100.0):.2f}%</td>
      <td><span style="font-size:11px;color:#d97706;font-weight:600">Flagged in GATK Mutect2</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Unmapped Reads</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{unmapped:,}</td>
      <td style="font-family:monospace;color:#ef4444;font-weight:700">{(unmapped/tot_num*100.0):.2f}%</td>
      <td><span style="font-size:11px;color:#ef4444;font-weight:600">&lt; 20.0% Expected</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Supplementary Alignments</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{supp:,}</td>
      <td style="font-family:monospace;color:var(--muted)">{(supp/tot_num*100.0):.2f}%</td>
      <td><span class="badge-blue">&lt; 1.0% Low Chimeric Rate</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Low Mapping Quality (MQ0)</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{mq0:,}</td>
      <td style="font-family:monospace;color:var(--muted)">{(mq0/tot_num*100.0):.4f}%</td>
      <td><span class="badge-green">&lt; 0.1% High Specificity</span></td>
    </tr>
    """

    # Table 2: Paired-End Sequencing & Library Geometry
    pair_tr_html = f"""
    <tr>
      <td style="font-weight:700;color:var(--text)">Total Paired Reads</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{paired:,} <span style="font-size:10.5px;color:var(--muted);font-weight:400">({r1:,} R1 + {r2:,} R2)</span></td>
      <td style="font-family:monospace;color:var(--muted)">100.0%</td>
      <td><span class="badge-green">Paired-End 2&times;150</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Properly Paired Reads</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{prop_paired:,}</td>
      <td style="font-family:monospace;color:#10b981;font-weight:700">{(prop_paired/tot_num*100.0):.2f}%</td>
      <td><span class="badge-green">&ge; 85.0% Optimal</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Inward Oriented (FR) Pairs</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{inward:,}</td>
      <td style="font-family:monospace;color:#0284c7;font-weight:700">{(inward/tot_pairs*100.0):.2f}%</td>
      <td><span class="badge-blue">&gt; 95.0% Standard Topology</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Outward Oriented (RF) Pairs</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{outward:,}</td>
      <td style="font-family:monospace;color:#8b5cf6;font-weight:700">{(outward/tot_pairs*100.0):.2f}%</td>
      <td><span style="font-size:11px;color:var(--muted);font-weight:600">Circular mtDNA Crossings</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Other / Non-Standard Pairs</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{other_pairs:,}</td>
      <td style="font-family:monospace;color:var(--muted)">{(other_pairs/tot_pairs*100.0):.3f}%</td>
      <td><span class="badge-green">&lt; 0.1% Rare Artifacts</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Mean Insert Size</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{avg_is:.1f} &plusmn; {std_is:.1f} bp</td>
      <td style="font-family:monospace;color:var(--muted)">&mdash;</td>
      <td><span class="badge-green">300&ndash;450 bp Optimal</span></td>
    </tr>
    """

    # Table 3: Base Yield & Sequencing Accuracy
    tot_gb = tot_len / 1e9
    map_gb = bases_mapped / 1e9
    yield_tr_html = f"""
    <tr>
      <td style="font-weight:700;color:var(--text)">Total Sequenced Base Yield</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{tot_gb:.2f} Gb <span style="font-size:10.5px;color:var(--muted);font-weight:400">({tot_len:,} bp)</span></td>
      <td style="font-family:monospace;color:var(--muted)">100.0%</td>
      <td><span class="badge-green">High Throughput</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Aligned Base Yield</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{map_gb:.2f} Gb <span style="font-size:10.5px;color:var(--muted);font-weight:400">({bases_mapped:,} bp)</span></td>
      <td style="font-family:monospace;color:#10b981;font-weight:700">{(bases_mapped/tot_len*100.0 if tot_len>0 else 0):.2f}%</td>
      <td><span class="badge-green">&ge; 80.0% High Yield</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Average Read Length</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{avg_len:.1f} bp <span style="font-size:10.5px;color:var(--muted);font-weight:400">(Max: {max_len} bp)</span></td>
      <td style="font-family:monospace;color:var(--muted)">&mdash;</td>
      <td><span class="badge-blue">147&ndash;151 bp (2&times;150 PE)</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Mean Base Phred Quality</td>
      <td style="font-family:monospace;font-weight:700;color:#10b981">Q{avg_q:.1f}</td>
      <td style="font-family:monospace;color:#10b981;font-weight:700">&gt; 99.95% Accuracy</td>
      <td><span class="badge-green">&ge; Q30 High Quality</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Total Alignment Mismatches</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{mismatches:,} bp</td>
      <td style="font-family:monospace;color:var(--muted)">{(mismatches/bases_mapped*100.0 if bases_mapped>0 else 0):.3f}%</td>
      <td><span class="badge-green">Expected Polymorphisms</span></td>
    </tr>
    <tr>
      <td style="font-weight:700;color:var(--text)">Sequencing Error Rate</td>
      <td style="font-family:monospace;font-weight:700;color:var(--text)">{(err_rate*100.0):.3f}% <span style="font-size:10.5px;color:var(--muted);font-weight:400">({err_rate:.3e})</span></td>
      <td style="font-family:monospace;color:var(--muted)">&mdash;</td>
      <td><span class="badge-green">&lt; 1.0% High Fidelity</span></td>
    </tr>
    """

    # Table 4: Idxstats + Competitive NUMT Contig Distribution rows
    idx_tr_html = ""
    for r in idxstats_rows:
        idx_tr_html += (
            f"<tr>"
            f"<td style='font-weight:700;color:var(--teal)'>{html.escape(str(r['ref']))} <span style='font-size:10.5px;color:var(--muted);font-weight:600'>(Mitochondrial rCRS)</span></td>"
            f"<td style='font-family:monospace'>{r['length']:,} bp</td>"
            f"<td><span class='badge-green'>{r['mapped']:,} reads</span></td>"
            f"<td><span class='badge-blue'>Retained in BAM</span></td>"
            f"</tr>\n"
        )
    if is_numt_competitive:
        idx_tr_html += (
            f"<tr style='background:rgba(220,38,38,0.05)'>"
            f"<td style='font-weight:700;color:#dc2626'>All {numt_loci_hit} Active Nuclear NUMT Contigs <span style='font-size:10.5px;color:var(--muted);font-weight:600'>(766 UCSC/RHNumtS.2 Loci)</span></td>"
            f"<td style='font-family:monospace'>629,517 bp</td>"
            f"<td><span style='background:rgba(220,38,38,0.14);color:#dc2626;padding:2px 8px;border-radius:12px;font-weight:700;font-size:11px'>{numt_depleted:,} reads ({numt_dep_pct:.2f}%)</span></td>"
            f"<td><span style='font-size:11px;font-weight:700;color:#dc2626'>Depleted Pre-Calling</span></td>"
            f"</tr>\n"
        )
        for rank_i, nl in enumerate(top_numt_loci, start=1):
            loc_s = str(nl.get("locus") or nl.get("contig") or "")
            nid_s = str(nl.get("numt_id") or "")
            r_abs = int(nl.get("reads_absorbed", 0))
            m_len = re.search(r':(\d+)-(\d+)', loc_s)
            span_bp = f"{(int(m_len.group(2)) - int(m_len.group(1))):,} bp" if m_len else "—"
            pct_of_numt = (r_abs / numt_depleted * 100.0) if numt_depleted > 0 else 0.0
            idx_tr_html += (
                f"<tr>"
                f"<td style='font-family:monospace;font-weight:700;color:var(--text)'>#{rank_i} {html.escape(loc_s)} <span style='font-size:10px;color:var(--muted);font-weight:400'>[{html.escape(nid_s)}]</span></td>"
                f"<td style='font-family:monospace'>{span_bp}</td>"
                f"<td style='font-family:monospace;font-weight:700;color:#ea580c'>{r_abs:,} reads <span style='font-size:10px;color:var(--muted);font-weight:400'>({pct_of_numt:.1f}% of NUMT)</span></td>"
                f"<td><span style='font-size:10.5px;color:#d97706;font-weight:600'>NUMT Contig Depleted</span></td>"
                f"</tr>\n"
            )

    date_str = datetime.now().strftime("%d %b %Y, %H:%M")

    sample_esc = html.escape(sample)
    total_seq_disp = f"{int(float(total_reads)):,}"
    avg_len_disp = f"{float(avg_length):.1f} bp"
    avg_is_disp = f"{float(avg_insert_size):.1f} bp"
    is_std_disp = f"{float(insert_std):.1f}"
    align_mode_disp = (
        f"Competitive rCRS + 766 Nuclear NUMTs ({numt_threads} Threads/Job)"
        if is_numt_competitive
        else f"Canonical rCRS (NC_012920.1) ({numt_threads} Threads/Job)"
    )
    numt_kpi_card_html = (
        f"""
      <div class="stat-card" style="border-top-color:#dc2626">
        <div class="stat-val" style="color:#dc2626">{numt_depleted:,} <span style="font-size:13px;font-weight:700">({numt_dep_pct:.2f}%)</span></div>
        <div class="stat-lbl">NUMT Reads Depleted ({numt_loci_hit} Loci)</div>
      </div>"""
        if is_numt_competitive
        else ""
    )

    html_template = """<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>MITOCLIN BAM QC Statistics &amp; Interactive Plots — __SAMPLE__</title>
<script>
  (function() {
    const saved = localStorage.getItem('mitoclin_theme');
    if (saved === 'dark' || (!saved && window.matchMedia('(prefers-color-scheme: dark)').matches)) {
      document.documentElement.setAttribute('data-theme', 'dark');
    }
  })();
</script>
<style>
  :root, [data-theme="light"] {
    --navy: #0b2135;
    --navy2: #162f48;
    --teal: #0d9488;
    --cyan: #0284c7;
    --blue2: #0284c7;
    --border: #e2e8f0;
    --bg: #f8fafc;
    --card-bg: #ffffff;
    --card-header-bg: #f8fafc;
    --text: #0f172a;
    --text-heading: #0f172a;
    --muted: #64748b;
    --shadow: 0 4px 20px rgba(0,0,0,0.06);
    --tooltip-bg: #0f172a;
    --tooltip-text: #ffffff;
    --tooltip-border: #334155;
  }

  [data-theme="dark"] {
    --navy: #090d16;
    --navy2: #0f172a;
    --teal: #2dd4bf;
    --cyan: #38bdf8;
    --blue2: #38bdf8;
    --border: #1e293b;
    --bg: #090d16;
    --card-bg: #131d31;
    --card-header-bg: #0f172a;
    --text: #f1f5f9;
    --text-heading: #ffffff;
    --muted: #94a3b8;
    --shadow: 0 4px 20px rgba(0,0,0,0.4);
    --tooltip-bg: #0b1329;
    --tooltip-text: #f8fafc;
    --tooltip-border: #1e293b;
  }

  * { box-sizing: border-box; margin: 0; padding: 0; }
  body {
    font-family: 'Segoe UI', -apple-system, BlinkMacSystemFont, Roboto, sans-serif;
    background: var(--bg);
    color: var(--text);
    padding: 24px;
    line-height: 1.5;
    transition: background 0.2s, color 0.2s;
  }
  .container {
    max-width: 1200px;
    margin: 0 auto;
    background: var(--card-bg);
    border-radius: 12px;
    border: 1px solid var(--border);
    box-shadow: var(--shadow);
    overflow: hidden;
  }
  .header {
    background: linear-gradient(135deg, #0b2135 0%, #162f48 100%);
    padding: 24px 28px;
    color: #ffffff;
    display: flex;
    justify-content: space-between;
    align-items: center;
    border-bottom: 3px solid var(--teal);
  }
  .brand-title { font-size: 20px; font-weight: 800; letter-spacing: 0.5px; }
  .brand-title span { color: #10b981; }
  .header-sub { font-size: 11px; color: #94a3b8; margin-top: 3px; }
  
  .meta-bar {
    background: var(--card-header-bg);
    padding: 10px 28px;
    font-size: 12px;
    color: var(--muted);
    border-bottom: 1px solid var(--border);
    display: flex;
    justify-content: space-between;
    align-items: center;
    flex-wrap: wrap;
    gap: 12px;
  }
  .meta-bar b { color: var(--text-heading); }

  .theme-toggle-btn {
    background: rgba(13,148,136,0.12);
    border: 1px solid var(--border);
    color: var(--text);
    padding: 5px 12px;
    border-radius: 16px;
    font-size: 11.5px;
    font-weight: 700;
    cursor: pointer;
    display: inline-flex;
    align-items: center;
    gap: 6px;
  }

  .content { padding: 24px 28px; }
  
  /* STAT CARDS */
  .stat-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
    gap: 12px;
    margin-bottom: 24px;
  }
  .stat-card {
    background: var(--card-bg);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 14px 16px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.04);
    position: relative;
    border-top: 3px solid var(--teal);
  }
  .stat-val {
    font-size: 22px;
    font-weight: 800;
    color: var(--text-heading);
    font-family: 'Consolas', monospace;
    line-height: 1.1;
  }
  .stat-lbl {
    font-size: 11px;
    color: var(--muted);
    margin-top: 5px;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.5px;
  }

  /* PLOTS GRID */
  .plots-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(500px, 1fr));
    gap: 18px;
    margin-bottom: 24px;
  }
  @media(max-width: 600px) {
    .plots-grid { grid-template-columns: 1fr; }
  }
  
  /* SECTIONS */
  .section {
    margin-bottom: 24px;
    background: var(--card-bg);
    border: 1px solid var(--border);
    border-radius: 8px;
    overflow: hidden;
  }
  .section-header {
    background: var(--card-header-bg);
    padding: 12px 18px;
    font-size: 13px;
    font-weight: 700;
    color: var(--text-heading);
    border-bottom: 1px solid var(--border);
    display: flex;
    align-items: center;
    justify-content: space-between;
  }
  .section-body { padding: 16px 18px; }

  /* PLOT TOOLBAR BUTTONS */
  .plot-toolbar {
    display: inline-flex;
    align-items: center;
    gap: 6px;
  }
  .plot-btn {
    background: var(--card-bg);
    border: 1px solid var(--border);
    color: var(--muted);
    padding: 2px 7px;
    border-radius: 4px;
    font-size: 10.5px;
    font-weight: 700;
    cursor: pointer;
    display: inline-flex;
    align-items: center;
    gap: 4px;
    transition: all 0.15s;
  }
  .plot-btn:hover {
    color: var(--teal);
    border-color: var(--teal);
    background: rgba(13,148,136,0.06);
  }
  
  /* TABLES */
  table {
    width: 100%;
    border-collapse: collapse;
    font-size: 12px;
  }
  th {
    background: var(--card-header-bg);
    text-align: left;
    padding: 8px 12px;
    font-weight: 700;
    color: var(--muted);
    border-bottom: 2px solid var(--border);
  }
  td {
    padding: 8px 12px;
    border-bottom: 1px solid var(--border);
    vertical-align: middle;
    color: var(--text);
  }
  tr:hover td { background: rgba(13,148,136,0.05); }
  .badge-green {
    background: rgba(16,185,129,0.15);
    color: #10b981;
    padding: 2px 8px;
    border-radius: 12px;
    font-weight: 700;
    font-size: 11px;
  }
  .badge-blue {
    background: rgba(2,132,199,0.15);
    color: #0284c7;
    padding: 2px 8px;
    border-radius: 12px;
    font-weight: 700;
    font-size: 11px;
  }
  .footer {
    text-align: center;
    padding: 16px;
    font-size: 11px;
    color: var(--muted);
    border-top: 1px solid var(--border);
    background: var(--card-header-bg);
  }

  /* INTERACTIVE TOOLTIP */
  #bamPlotTooltip {
    position: absolute;
    z-index: 10000;
    display: none;
    background: var(--tooltip-bg);
    color: var(--tooltip-text);
    border: 1px solid var(--tooltip-border);
    padding: 8px 12px;
    border-radius: 6px;
    font-size: 11px;
    line-height: 1.45;
    box-shadow: 0 8px 24px rgba(0,0,0,0.35);
    pointer-events: none;
    max-width: 280px;
    backdrop-filter: blur(4px);
    transition: opacity 0.08s ease-out;
  }

  .interactive-bar:hover {
    opacity: 0.82;
    filter: brightness(1.15);
  }
  .donut-slice:hover {
    stroke-width: 26px !important;
    filter: drop-shadow(0 2px 6px rgba(0,0,0,0.25));
  }
  .acgt-toggle.inactive, .indel-toggle.inactive {
    opacity: 0.45;
    text-decoration: line-through;
    filter: grayscale(0.8);
  }
</style>
</head>
<body>

<!-- GLOBAL FLOATING TOOLTIP -->
<div id="bamPlotTooltip"></div>

<div class="container">
  <div class="header">
    <div>
      <div class="brand-title">MITO<span>CLIN</span> &mdash; BAM QC &amp; Interactive Alignment Plots</div>
      <div class="header-sub">Centre for DNA Fingerprinting and Diagnostics (CDFD) &bull; Comprehensive samtools Metrics</div>
    </div>
    <div style="text-align:right">
      <div style="font-size:15px;font-weight:800;letter-spacing:0.3px">__SAMPLE__</div>
      <div style="font-size:11px;color:#94a3b8">BAM Quality Control Report</div>
    </div>
  </div>

  <div class="meta-bar">
    <div style="display:flex;gap:18px;flex-wrap:wrap">
      <span>Sample: <b>__SAMPLE__</b></span>
      <span>Generated: <b>__DATE_STR__</b></span>
      <span>Pipeline: <b>MITOCLIN v1.0</b></span>
      <span>Alignment Mode: <b style="color:var(--teal)">__ALIGN_MODE_DISP__</b></span>
      <span>Canonical Target: <b style="font-family:monospace">rCRS (NC_012920.1)</b></span>
    </div>
    <div>
      <button class="theme-toggle-btn" onclick="toggleTheme()">
        <span id="theme-btn-icon">🌓</span> <span id="theme-btn-label">Theme</span>
      </button>
    </div>
  </div>

  <div class="content">
    <!-- STATS OVERVIEW -->
    <div class="stat-grid">
      <div class="stat-card">
        <div class="stat-val">__TOTAL_SEQ__</div>
        <div class="stat-lbl">rCRS Retained Reads</div>
      </div>
      __NUMT_KPI_CARD__
      <div class="stat-card" style="border-top-color:#10b981">
        <div class="stat-val">__PCT_MAPPED__%</div>
        <div class="stat-lbl">Overall Mapping Rate</div>
      </div>
      <div class="stat-card" style="border-top-color:#f59e0b">
        <div class="stat-val">__PCT_DUP__%</div>
        <div class="stat-lbl">Duplication Rate</div>
      </div>
      <div class="stat-card" style="border-top-color:#0284c7">
        <div class="stat-val">__AVG_LEN__</div>
        <div class="stat-lbl">Average Read Length</div>
      </div>
      <div class="stat-card" style="border-top-color:#8b5cf6">
        <div class="stat-val">__AVG_IS__</div>
        <div class="stat-lbl">Mean Insert Size</div>
      </div>
      <div class="stat-card" style="border-top-color:#ef4444">
        <div class="stat-val">__ERR_STR__</div>
        <div class="stat-lbl">Error Rate</div>
      </div>
    </div>

    <!-- ═══ VISUAL PLOTS SECTION ═══ -->
    <div style="margin-bottom:16px;display:flex;align-items:center;justify-content:space-between">
      <div style="display:flex;align-items:center;gap:8px">
        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="#0d9488" stroke-width="2.5"><polyline points="22 12 18 12 15 21 9 3 6 12 2 12"/></svg>
        <h2 style="font-size:16px;font-weight:800;color:var(--text-heading)">Interactive Alignment &amp; Quality Control Plots</h2>
      </div>
      <div style="font-size:11px;color:var(--muted)">
        <i>Hover over plots to inspect exact values &bull; Click legends to filter series</i>
      </div>
    </div>

    <div class="plots-grid">
      <!-- PLOT 1: ALIGNMENT BREAKDOWN -->
      <div class="section">
        <div class="section-header">
          <span>Alignment, NUMT Depletion &amp; Duplication Breakdown</span>
          <div class="plot-toolbar">
            <span class="badge-green">__PCT_MAPPED_1DEC__% Mapped</span>
            <button class="plot-btn" onclick="exportPlotPng('donut_align_svg', '__SAMPLE___alignment_breakdown.png')" title="Download PNG">PNG</button>
            <button class="plot-btn" onclick="exportPlotSvg('donut_align_svg', '__SAMPLE___alignment_breakdown.svg')" title="Download SVG">SVG</button>
          </div>
        </div>
        <div class="section-body">
          __SVG_DONUT_ALIGN__
        </div>
      </div>

      <!-- PLOT 2: PAIR ORIENTATION BREAKDOWN -->
      <div class="section">
        <div class="section-header">
          <span>Paired-End Mate Orientation Breakdown</span>
          <div class="plot-toolbar">
            <span class="badge-blue">Inward (FR) Pairs</span>
            <button class="plot-btn" onclick="exportPlotPng('donut_pairs_svg', '__SAMPLE___pair_orientation.png')" title="Download PNG">PNG</button>
            <button class="plot-btn" onclick="exportPlotSvg('donut_pairs_svg', '__SAMPLE___pair_orientation.svg')" title="Download SVG">SVG</button>
          </div>
        </div>
        <div class="section-body">
          __SVG_DONUT_PAIRS__
        </div>
      </div>

      __NUMT_PLOT_CARD__

      <!-- PLOT 3: INSERT SIZE DISTRIBUTION -->
      <div class="section">
        <div class="section-header">
          <span>Insert Size Distribution (IS)</span>
          <div class="plot-toolbar">
            <span style="font-size:11px;color:var(--muted);font-family:monospace">Mean: __AVG_IS__ &plusmn; __IS_STD__ bp</span>
            <button class="plot-btn" onclick="exportPlotPng('plot_is_svg', '__SAMPLE___insert_size.png')" title="Download PNG">PNG</button>
            <button class="plot-btn" onclick="exportPlotSvg('plot_is_svg', '__SAMPLE___insert_size.svg')" title="Download SVG">SVG</button>
          </div>
        </div>
        <div class="section-body">
          __SVG_PLOT_IS__
        </div>
      </div>

      <!-- PLOT 4: PER-CYCLE PHRED QUALITY SCORES -->
      <div class="section">
        <div class="section-header">
          <span>Per-Cycle Mean Base Quality Score (Phred Q-Score)</span>
          <div class="plot-toolbar">
            <span class="badge-green">&ge; Q30 (High Quality)</span>
            <button class="plot-btn" onclick="exportPlotPng('plot_quality_svg', '__SAMPLE___quality_score.png')" title="Download PNG">PNG</button>
            <button class="plot-btn" onclick="exportPlotSvg('plot_quality_svg', '__SAMPLE___quality_score.svg')" title="Download SVG">SVG</button>
          </div>
        </div>
        <div class="section-body">
          __SVG_PLOT_QUALITY__
        </div>
      </div>

      <!-- PLOT 5: PER-CYCLE NUCLEOTIDE COMPOSITION -->
      <div class="section">
        <div class="section-header">
          <span>Per-Cycle Base Composition (% ACGTN)</span>
          <div class="plot-toolbar">
            <span style="font-size:11px;color:var(--muted)">Cycle Equilibrium</span>
            <button class="plot-btn" onclick="exportPlotPng('plot_acgt_svg', '__SAMPLE___base_composition.png')" title="Download PNG">PNG</button>
            <button class="plot-btn" onclick="exportPlotSvg('plot_acgt_svg', '__SAMPLE___base_composition.svg')" title="Download SVG">SVG</button>
          </div>
        </div>
        <div class="section-body">
          __SVG_PLOT_ACGT__
        </div>
      </div>

      <!-- PLOT 6: READ LENGTH DISTRIBUTION -->
      <div class="section">
        <div class="section-header">
          <span>Read Length Distribution (RL)</span>
          <div class="plot-toolbar">
            <span style="font-size:11px;color:var(--muted);font-family:monospace">Mean: __AVG_LEN__</span>
            <button class="plot-btn" onclick="exportPlotPng('plot_rl_svg', '__SAMPLE___read_length.png')" title="Download PNG">PNG</button>
            <button class="plot-btn" onclick="exportPlotSvg('plot_rl_svg', '__SAMPLE___read_length.svg')" title="Download SVG">SVG</button>
          </div>
        </div>
        <div class="section-body">
          __SVG_PLOT_RL__
        </div>
      </div>

      <!-- PLOT 7: GC CONTENT DISTRIBUTION -->
      <div class="section">
        <div class="section-header">
          <span>GC Content Distribution (GCF)</span>
          <div class="plot-toolbar">
            <span class="badge-green">mtDNA Profile</span>
            <button class="plot-btn" onclick="exportPlotPng('plot_gc_svg', '__SAMPLE___gc_content.png')" title="Download PNG">PNG</button>
            <button class="plot-btn" onclick="exportPlotSvg('plot_gc_svg', '__SAMPLE___gc_content.svg')" title="Download SVG">SVG</button>
          </div>
        </div>
        <div class="section-body">
          __SVG_PLOT_GC__
        </div>
      </div>

      <!-- PLOT 8: INDEL SIZE DISTRIBUTION -->
      <div class="section">
        <div class="section-header">
          <span>Indel Size Distribution (Insertions vs Deletions)</span>
          <div class="plot-toolbar">
            <span style="font-size:11px;color:var(--muted)">Alignment Indels</span>
            <button class="plot-btn" onclick="exportPlotPng('plot_indel_svg', '__SAMPLE___indel_distribution.png')" title="Download PNG">PNG</button>
            <button class="plot-btn" onclick="exportPlotSvg('plot_indel_svg', '__SAMPLE___indel_distribution.svg')" title="Download SVG">SVG</button>
          </div>
        </div>
        <div class="section-body">
          __SVG_PLOT_INDEL__
        </div>
      </div>
    </div>

    <!-- ═══ STRUCTURED METRIC TABLES SECTION ═══ -->
    <div style="margin:28px 0 16px 0;display:flex;align-items:center;gap:8px">
      <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="#0d9488" stroke-width="2.5"><line x1="8" y1="6" x2="21" y2="6"/><line x1="8" y1="12" x2="21" y2="12"/><line x1="8" y1="18" x2="21" y2="18"/><line x1="3" y1="6" x2="3.01" y2="6"/><line x1="3" y1="12" x2="3.01" y2="12"/><line x1="3" y1="18" x2="3.01" y2="18"/></svg>
      <h2 style="font-size:16px;font-weight:800;color:var(--text-heading)">Alignment &amp; Sequencing Quality Metrics Summary</h2>
    </div>

    <div style="display:grid;grid-template-columns:repeat(auto-fit, minmax(540px, 1fr));gap:20px;margin-bottom:24px">
      <!-- TABLE 1: PRIMARY ALIGNMENT & DUPLICATION -->
      <div class="section" style="margin-bottom:0">
        <div class="section-header">
          <div style="display:flex;align-items:center;gap:6px">
            <svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><polyline points="22 12 18 12 15 21 9 3 6 12 2 12"/></svg>
            <span>Primary Alignment, NUMT Depletion &amp; Duplication Summary</span>
          </div>
        </div>
        <div class="section-body" style="padding:0">
          <table>
            <thead>
              <tr><th>Metric</th><th>Count</th><th>Rate (%)</th><th>Benchmark / Status</th></tr>
            </thead>
            <tbody>
              __ALIGN_ROWS__
            </tbody>
          </table>
        </div>
      </div>

      <!-- TABLE 2: PAIRED-END SEQUENCING & GEOMETRY -->
      <div class="section" style="margin-bottom:0">
        <div class="section-header">
          <div style="display:flex;align-items:center;gap:6px">
            <svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="12" cy="12" r="10"/><path d="M12 8v8M8 12h8"/></svg>
            <span>Paired-End Sequencing &amp; Library Topology</span>
          </div>
        </div>
        <div class="section-body" style="padding:0">
          <table>
            <thead>
              <tr><th>Metric</th><th>Value</th><th>Rate (%)</th><th>Benchmark</th></tr>
            </thead>
            <tbody>
              __PAIR_ROWS__
            </tbody>
          </table>
        </div>
      </div>

      <!-- TABLE 3: BASE YIELD & ACCURACY -->
      <div class="section" style="margin-bottom:0">
        <div class="section-header">
          <div style="display:flex;align-items:center;gap:6px">
            <svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M14 2H6a2 2 0 00-2 2v16a2 2 0 002 2h12a2 2 0 002-2V8z"/><polyline points="14 2 14 8 20 8"/><line x1="16" y1="13" x2="8" y2="13"/><line x1="16" y1="17" x2="8" y2="17"/></svg>
            <span>Base Yield, Length &amp; Error Rate</span>
          </div>
        </div>
        <div class="section-body" style="padding:0">
          <table>
            <thead>
              <tr><th>Metric</th><th>Value</th><th>Proportion</th><th>Benchmark</th></tr>
            </thead>
            <tbody>
              __YIELD_ROWS__
            </tbody>
          </table>
        </div>
      </div>

      <!-- TABLE 4: CHROMOSOME & NUCLEAR NUMT CONTIG DISTRIBUTION -->
      <div class="section" style="margin-bottom:0">
        <div class="section-header" style="flex-wrap:wrap;gap:8px">
          <div style="display:flex;align-items:center;gap:6px">
            <svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 16V8a2 2 0 0 0-1-1.73l-7-4a2 2 0 0 0-2 0l-7 4A2 2 0 0 0 3 8v8a2 2 0 0 0 1 1.73l7 4a2 2 0 0 0 2 0l7-4A2 2 0 0 0 21 16z"/></svg>
            <span>Chromosome &amp; Nuclear NUMT Contig Distribution</span>
          </div>
          <input type="text" id="contigSearchInput" oninput="filterContigRows()" placeholder="Filter contig / chr locus..." style="padding:3px 8px;font-size:11px;border:1px solid var(--border);border-radius:4px;background:var(--card-bg);color:var(--text);width:170px">
        </div>
        <div class="section-body" style="padding:0">
          <table id="contigDistTable">
            <thead>
              <tr><th>Reference / Nuclear NUMT Contig</th><th>Contig Span</th><th>Mapped / Absorbed Reads</th><th>Pipeline Action</th></tr>
            </thead>
            <tbody>
              __IDXSTATS_ROWS__
            </tbody>
          </table>
        </div>
      </div>
    </div>

    <!-- ══════════════════════════════════════════════════════════════
         QUALITY METRICS GLOSSARY & REFERENCE STANDARDS
         ══════════════════════════════════════════════════════════════ -->
    <div style="margin:28px 0 16px 0;display:flex;align-items:center;gap:8px">
      <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="#0d9488" stroke-width="2.5"><path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"/><path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"/></svg>
      <h2 style="font-size:16px;font-weight:800;color:var(--text-heading)">Quality Metrics Glossary &amp; Clinical Sequencing Definitions</h2>
    </div>

    <div style="display:grid;grid-template-columns:repeat(auto-fit, minmax(340px, 1fr));gap:16px;margin-bottom:28px">
      <!-- 1. TOTAL SEQUENCES -->
      <div style="background:var(--card-bg);border:1px solid var(--border);border-left:4px solid #0d9488;border-radius:6px;padding:14px">
        <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:6px">
          <span style="font-weight:800;font-size:13px;color:var(--text-heading)">rCRS Retained Reads</span>
          <span class="badge-green">__TOTAL_SEQ__ Reads</span>
        </div>
        <div style="font-size:11.5px;color:var(--muted);line-height:1.5">
          <b>Definition:</b> The total number of authentic mitochondrial DNA reads retained on <code>NC_012920.1</code> (rCRS) after competitive nuclear NUMT depletion.
          <br/><b>Clinical Utility:</b> High read counts (&gt; 1M reads) ensure ultra-deep depth across the 16,569 bp mitochondrial genome to detect low-level heteroplasmic variants down to &le; 5% allele frequency.
        </div>
      </div>

      <!-- 2. COMPETITIVE NUMT DEPLETION -->
      <div style="background:var(--card-bg);border:1px solid var(--border);border-left:4px solid #dc2626;border-radius:6px;padding:14px">
        <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:6px">
          <span style="font-weight:800;font-size:13px;color:var(--text-heading)">Competitive NUMT Depletion</span>
          <span style="background:rgba(220,38,38,0.14);color:#dc2626;padding:2px 8px;border-radius:12px;font-weight:700;font-size:11px">766 Nuclear NUMT Loci</span>
        </div>
        <div style="font-size:11.5px;color:var(--muted);line-height:1.5">
          <b>Definition:</b> Simultaneous BWA-MEM alignment against <code>NC_012920.1</code> (rCRS) and 543 merged human nuclear NUMT chromosomal contigs (766 UCSC/RHNumtS.2 loci across <code>chr1..chrY</code>).
          <br/><b>Clinical Utility:</b> Depletes nuclear pseudogene carryover and chimeric rCRS&ndash;NUMT mate pairs prior to GATK Mutect2 calling, preventing false-positive low-VAF heteroplasmy artifacts.
        </div>
      </div>

      <!-- 3. ALIGNMENT RATE -->
      <div style="background:var(--card-bg);border:1px solid var(--border);border-left:4px solid #10b981;border-radius:6px;padding:14px">
        <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:6px">
          <span style="font-weight:800;font-size:13px;color:var(--text-heading)">Overall Mapping Rate</span>
          <span class="badge-green">__PCT_MAPPED__% Mapped</span>
        </div>
        <div style="font-size:11.5px;color:var(--muted);line-height:1.5">
          <b>Definition:</b> The percentage of sequenced reads that successfully align to the mitochondrial reference (and homologous nuclear NUMT contigs in competitive mode).
          <br/><b>Clinical Utility:</b> High mapping rates confirm strong mtDNA enrichment efficiency and minimal exogenous bacterial or adapter-dimer contamination.
        </div>
      </div>

      <!-- 4. DUPLICATION RATE -->
      <div style="background:var(--card-bg);border:1px solid var(--border);border-left:4px solid #f59e0b;border-radius:6px;padding:14px">
        <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:6px">
          <span style="font-weight:800;font-size:13px;color:var(--text-heading)">Duplication Rate</span>
          <span style="background:rgba(245,158,11,0.15);color:#d97706;padding:2px 8px;border-radius:12px;font-weight:700;font-size:11px">__PCT_DUP__% Duplicate</span>
        </div>
        <div style="font-size:11.5px;color:var(--muted);line-height:1.5">
          <b>Definition:</b> The fraction of mapped reads sharing identical 5&prime; and 3&prime; genomic start/stop coordinates due to PCR amplification during library preparation.
          <br/><b>Clinical Utility:</b> PCR duplicates are flagged by GATK MarkDuplicates so variant allele fraction (VAF) calculations and heteroplasmy calls are not skewed by amplification bias.
        </div>
      </div>

      <!-- 5. AVERAGE READ LENGTH -->
      <div style="background:var(--card-bg);border:1px solid var(--border);border-left:4px solid #0284c7;border-radius:6px;padding:14px">
        <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:6px">
          <span style="font-weight:800;font-size:13px;color:var(--text-heading)">Average Read Length</span>
          <span class="badge-blue">__AVG_LEN__</span>
        </div>
        <div style="font-size:11.5px;color:var(--muted);line-height:1.5">
          <b>Definition:</b> The mean physical length of sequenced fragments across all sequencing cycles.
          <br/><b>Clinical Utility:</b> Standard Illumina 2&times;150 bp paired-end chemistry yields reads averaging ~135&ndash;151 bp after adapter trimming, providing optimal span across repetitive mtDNA control regions.
        </div>
      </div>

      <!-- 6. MEAN INSERT SIZE -->
      <div style="background:var(--card-bg);border:1px solid var(--border);border-left:4px solid #8b5cf6;border-radius:6px;padding:14px">
        <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:6px">
          <span style="font-weight:800;font-size:13px;color:var(--text-heading)">Mean Insert Size</span>
          <span style="background:rgba(139,92,246,0.15);color:#7c3aed;padding:2px 8px;border-radius:12px;font-weight:700;font-size:11px">__AVG_IS__</span>
        </div>
        <div style="font-size:11.5px;color:var(--muted);line-height:1.5">
          <b>Definition:</b> The physical genomic distance between the outer adapters of paired forward and reverse reads.
          <br/><b>Clinical Utility:</b> Consistent insert sizes ensure stable bridge amplification on the flowcell, uninhibited cluster generation, and clean mate-pair resolution.
        </div>
      </div>

      <!-- 7. ERROR RATE -->
      <div style="background:var(--card-bg);border:1px solid var(--border);border-left:4px solid #ef4444;border-radius:6px;padding:14px">
        <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:6px">
          <span style="font-weight:800;font-size:13px;color:var(--text-heading)">Sequencing Error Rate</span>
          <span style="background:rgba(239,68,68,0.15);color:#b91c1c;padding:2px 8px;border-radius:12px;font-weight:700;font-size:11px">__ERR_STR__</span>
        </div>
        <div style="font-size:11.5px;color:var(--muted);line-height:1.5">
          <b>Definition:</b> The proportion of mismatch bases relative to the aligned reference genome across all aligned CIGAR strings.
          <br/><b>Clinical Utility:</b> An error rate &lt; 1.0% verifies high optical fidelity and ensures true somatic/germline heteroplasmies are cleanly separated from sequencing noise.
        </div>
      </div>
    </div>

    <!-- ══════════════════════════════════════════════════════════════
         VISUAL PLOT INTERPRETATIONS & DIAGNOSTIC GUIDANCE
         ══════════════════════════════════════════════════════════════ -->
    <div style="margin:28px 0 16px 0;display:flex;align-items:center;gap:8px">
      <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="#0d9488" stroke-width="2.5"><circle cx="12" cy="12" r="10"/><path d="M12 16v-4"/><path d="M12 8h.01"/></svg>
      <h2 style="font-size:16px;font-weight:800;color:var(--text-heading)">Plot Interpretations &amp; Diagnostic Quality Guidance</h2>
    </div>

    <div style="display:grid;grid-template-columns:repeat(auto-fit, minmax(360px, 1fr));gap:16px;margin-bottom:32px">
      <!-- PLOT 1 -->
      <div style="background:var(--card-bg);border:1px solid var(--border);border-radius:8px;padding:16px">
        <div style="display:flex;align-items:center;gap:8px;margin-bottom:8px">
          <span style="font-weight:800;font-size:12.5px;color:var(--text-heading)">1. Alignment, NUMT Depletion &amp; Duplication Breakdown</span>
        </div>
        <div style="font-size:11.5px;color:var(--muted);line-height:1.55">
          <b>Purpose:</b> Quantifies the library composition into unique informative rCRS reads, marked PCR duplicates, depleted nuclear NUMT reads, and unmapped reads.
          <br/><b>Interpretation:</b> Unique rCRS mapped reads (green) provide independent mitochondrial evidence for variant calling, while nuclear NUMT reads (red) and duplicates (amber) are prevented from biasing VAF estimates.
        </div>
      </div>

      <!-- PLOT 2 -->
      <div style="background:var(--card-bg);border:1px solid var(--border);border-radius:8px;padding:16px">
        <div style="display:flex;align-items:center;gap:8px;margin-bottom:8px">
          <span style="font-weight:800;font-size:12.5px;color:var(--text-heading)">2. Paired-End Mate Orientation Breakdown</span>
        </div>
        <div style="font-size:11.5px;color:var(--muted);line-height:1.55">
          <b>Purpose:</b> Verifies the structural topology and pairing geometry of forward (R1) and reverse (R2) reads.
          <br/><b>Interpretation:</b> Standard Illumina paired-end libraries show inward-facing (FR) pairs (&gt; 95%). Outward (RF) pairs (&lt; 3%) represent reads spanning the circular mtDNA origin (m.16569/m.1 boundary).
        </div>
      </div>

      <!-- PLOT 3 -->
      <div style="background:var(--card-bg);border:1px solid var(--border);border-radius:8px;padding:16px">
        <div style="display:flex;align-items:center;gap:8px;margin-bottom:8px">
          <span style="font-weight:800;font-size:12.5px;color:var(--text-heading)">3. Top 10 Absorbing Nuclear NUMT Loci</span>
        </div>
        <div style="font-size:11.5px;color:var(--muted);line-height:1.55">
          <b>Purpose:</b> Identifies which human nuclear chromosomal NUMT segments (such as <code>chr1:564314-570454</code> or <code>chr17:19501722-19509396</code>) absorbed the most non-mitochondrial reads during competitive alignment.
          <br/><b>Interpretation:</b> Removing these nuclear pseudogene reads prior to GATK Mutect2 calling eliminates false-positive heteroplasmy calls in high-homology regions (e.g., <i>MT-CO1</i>, <i>MT-ATP6</i>, <i>MT-ND2</i>).
        </div>
      </div>

      <!-- PLOT 4 -->
      <div style="background:var(--card-bg);border:1px solid var(--border);border-radius:8px;padding:16px">
        <div style="display:flex;align-items:center;gap:8px;margin-bottom:8px">
          <span style="font-weight:800;font-size:12.5px;color:var(--text-heading)">4. Insert Size Distribution (IS)</span>
        </div>
        <div style="font-size:11.5px;color:var(--muted);line-height:1.55">
          <b>Purpose:</b> Assesses the consistency and uniformity of DNA fragmentation during library construction.
          <br/><b>Interpretation:</b> A smooth bell-shaped curve demonstrates uniform fragmentation and size selection without adapter-dimer contamination (&lt; 50 bp).
        </div>
      </div>

      <!-- PLOT 5 -->
      <div style="background:var(--card-bg);border:1px solid var(--border);border-radius:8px;padding:16px">
        <div style="display:flex;align-items:center;gap:8px;margin-bottom:8px">
          <span style="font-weight:800;font-size:12.5px;color:var(--text-heading)">5. Per-Cycle Mean Base Quality Scores (Phred Q)</span>
        </div>
        <div style="font-size:11.5px;color:var(--muted);line-height:1.55">
          <b>Purpose:</b> Tracks base-calling accuracy across all sequencing cycles from cycle 1 to 151.
          <br/><b>Interpretation:</b> Average Phred scores remaining well above Q30 (green zone, &ge; 99.9% accuracy) across the entire run confirm optical and fluorophore stability.
        </div>
      </div>

      <!-- PLOT 6 -->
      <div style="background:var(--card-bg);border:1px solid var(--border);border-radius:8px;padding:16px">
        <div style="display:flex;align-items:center;gap:8px;margin-bottom:8px">
          <span style="font-weight:800;font-size:12.5px;color:var(--text-heading)">6. Per-Cycle Base Composition (% ACGTN)</span>
        </div>
        <div style="font-size:11.5px;color:var(--muted);line-height:1.55">
          <b>Purpose:</b> Evaluates nucleotide diversity and detects sequence bias along read positions.
          <br/><b>Interpretation:</b> Constant, parallel lines matching mtDNA nucleotide frequencies (A: ~31%, C: ~31%, T: ~25%, G: ~13%) indicate an unbiased sequencing library.
        </div>
      </div>

      <!-- PLOT 7 -->
      <div style="background:var(--card-bg);border:1px solid var(--border);border-radius:8px;padding:16px">
        <div style="display:flex;align-items:center;gap:8px;margin-bottom:8px">
          <span style="font-weight:800;font-size:12.5px;color:var(--text-heading)">7. Read Length &amp; GC Content Distributions</span>
        </div>
        <div style="font-size:11.5px;color:var(--muted);line-height:1.55">
          <b>Purpose:</b> Monitors read length consistency following adapter trimming and evaluates amplification uniformity against expected human mitochondrial GC content (~44.4%).
          <br/><b>Interpretation:</b> A dominant peak near full read length and a Gaussian GC curve centered near ~44.4% confirm uniform mtDNA capture.
        </div>
      </div>

      <!-- PLOT 8 -->
      <div style="background:var(--card-bg);border:1px solid var(--border);border-radius:8px;padding:16px">
        <div style="display:flex;align-items:center;gap:8px;margin-bottom:8px">
          <span style="font-weight:800;font-size:12.5px;color:var(--text-heading)">8. Indel Size Distribution (Insertions vs Deletions)</span>
        </div>
        <div style="font-size:11.5px;color:var(--muted);line-height:1.55">
          <b>Purpose:</b> Evaluates the frequency and length distribution of alignment insertion and deletion events.
          <br/><b>Interpretation:</b> An exponential decay from 1 bp to longer indels reflects natural biological variations and polymerase slippage across homopolymer tracts (e.g., poly-C tracts in D-loop), verifying robust alignment gap handling.
        </div>
      </div>
    </div>
  </div>

  <div class="footer">
    Generated automatically by MITOCLIN Pipeline &bull; Centre for DNA Fingerprinting and Diagnostics &bull; For Research &amp; Clinical Analysis
  </div>
</div>

<!-- ══════════════════════════════════════════════════════════════
     CLIENT-SIDE INTERACTIVE PLOT ENGINE
     ══════════════════════════════════════════════════════════════ -->
<script>
  function toggleTheme() {
    const current = document.documentElement.getAttribute('data-theme') || 'light';
    const next = current === 'dark' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', next);
    try { localStorage.setItem('mitoclin_theme', next); } catch(e) {}
  }

  const tooltip = document.getElementById('bamPlotTooltip');

  function showTooltip(htmlContent, clientX, clientY) {
    if (!tooltip) return;
    tooltip.innerHTML = htmlContent;
    tooltip.style.display = 'block';

    const ttWidth = tooltip.offsetWidth || 220;
    const ttHeight = tooltip.offsetHeight || 80;
    const pad = 15;

    let x = clientX + pad;
    let y = clientY - ttHeight - 10;

    if (x + ttWidth > window.innerWidth - 10) {
      x = clientX - ttWidth - pad;
    }
    if (y < 10) {
      y = clientY + pad;
    }

    tooltip.style.left = Math.max(10, x + window.scrollX) + 'px';
    tooltip.style.top = Math.max(10, y + window.scrollY) + 'px';
  }

  function hideTooltip() {
    if (tooltip) tooltip.style.display = 'none';
  }

  // 1. DONUT INTERACTION
  document.querySelectorAll('.donut-slice, .donut-legend-item').forEach(el => {
    el.addEventListener('mouseenter', function(e) {
      const chartId = this.getAttribute('data-chart');
      const idx = this.getAttribute('data-idx');
      const slice = document.getElementById(`${chartId}_seg_${idx}`);
      if (!slice) return;

      const label = slice.getAttribute('data-label');
      const val = parseInt(slice.getAttribute('data-val') || 0);
      const pct = slice.getAttribute('data-pct');
      const color = slice.getAttribute('data-color');

      // Update center hub
      const hubLbl = document.getElementById(`${chartId}_hub_lbl`);
      const hubVal = document.getElementById(`${chartId}_hub_val`);
      const hubPct = document.getElementById(`${chartId}_hub_pct`);
      if (hubLbl) hubLbl.textContent = label.toUpperCase();
      if (hubVal) hubVal.textContent = val >= 1e6 ? (val/1e6).toFixed(2)+'M' : (val >= 1e4 ? (val/1e3).toFixed(1)+'k' : val.toLocaleString());
      if (hubPct) {
        hubPct.textContent = `${pct}%`;
        hubPct.style.display = 'block';
      }

      // Dim other slices
      document.querySelectorAll(`[id^="${chartId}_seg_"]`).forEach(s => {
        if (s.getAttribute('data-idx') === idx) {
          s.style.strokeWidth = '26px';
          s.style.strokeOpacity = '1';
        } else {
          s.style.strokeWidth = '22px';
          s.style.strokeOpacity = '0.35';
        }
      });

      // Show floating tooltip
      const tipHtml = `
        <div style="font-weight:700;color:${color};font-size:12px;margin-bottom:3px">${label}</div>
        <div><b>Count:</b> <span style="font-family:monospace;font-weight:700">${val.toLocaleString()} reads</span></div>
        <div><b>Share:</b> <span style="color:#10b981;font-weight:700">${pct}%</span> of total</div>
      `;
      showTooltip(tipHtml, e.clientX, e.clientY);
    });

    el.addEventListener('mouseleave', function() {
      const chartId = this.getAttribute('data-chart');
      document.querySelectorAll(`[id^="${chartId}_seg_"]`).forEach(s => {
        s.style.strokeWidth = '22px';
        s.style.strokeOpacity = '1';
      });
      const hubPct = document.getElementById(`${chartId}_hub_pct`);
      if (hubPct) hubPct.style.display = 'none';
      hideTooltip();
    });
  });

  // 2. CONTINUOUS LINE CHARTS (INSERT SIZE, GC, QUALITY) INTERACTION
  document.querySelectorAll('.interactive-plot-card[data-chart-type="line"]').forEach(card => {
    const chartId = card.getAttribute('data-chart-id');
    const svg = card.querySelector('svg');
    const xMin = parseFloat(card.getAttribute('data-xmin'));
    const xMax = parseFloat(card.getAttribute('data-xmax'));
    const yMin = parseFloat(card.getAttribute('data-ymin'));
    const yMax = parseFloat(card.getAttribute('data-ymax'));
    const padL = parseFloat(card.getAttribute('data-pad-left'));
    const padR = parseFloat(card.getAttribute('data-pad-right'));
    const padT = parseFloat(card.getAttribute('data-pad-top'));
    const padB = parseFloat(card.getAttribute('data-pad-bottom'));
    const width = parseFloat(card.getAttribute('data-width'));
    const height = parseFloat(card.getAttribute('data-height'));
    const unitX = card.getAttribute('data-unit-x') || '';
    const unitY = card.getAttribute('data-unit-y') || '';
    const lblX = card.getAttribute('data-label-x') || 'X';
    const lblY = card.getAttribute('data-label-y') || 'Y';
    const pts = JSON.parse(card.getAttribute('data-pts') || '[]');

    const crosshair = document.getElementById(`${chartId}_crosshair`);
    const guideX = document.getElementById(`${chartId}_guide_x`);
    const guideDot = document.getElementById(`${chartId}_guide_dot`);

    const plotW = width - padL - padR;
    const plotH = height - padT - padB;

    card.addEventListener('mousemove', function(e) {
      if (!pts || !pts.length) return;
      const rect = svg.getBoundingClientRect();
      const mouseX = ((e.clientX - rect.left) / rect.width) * width;
      const mouseY = ((e.clientY - rect.top) / rect.height) * height;

      if (mouseX < padL || mouseX > width - padR || mouseY < padT || mouseY > height - padB) {
        if (crosshair) crosshair.style.display = 'none';
        hideTooltip();
        return;
      }

      // Calculate corresponding X in data space
      const xVal = xMin + ((mouseX - padL) / plotW) * (xMax - xMin);

      // Find closest data point
      let closest = pts[0];
      let minDist = Math.abs(pts[0][0] - xVal);
      for (let i = 1; i < pts.length; i++) {
        const d = Math.abs(pts[i][0] - xVal);
        if (d < minDist) {
          minDist = d;
          closest = pts[i];
        }
      }

      // Convert closest point back to SVG coordinates
      const cx = padL + ((closest[0] - xMin) / (xMax - xMin)) * plotW;
      const cy = padT + plotH - ((closest[1] - yMin) / (yMax - yMin)) * plotH;

      if (crosshair) crosshair.style.display = 'block';
      if (guideX) {
        guideX.setAttribute('x1', cx);
        guideX.setAttribute('x2', cx);
      }
      if (guideDot) {
        guideDot.setAttribute('cx', cx);
        guideDot.setAttribute('cy', cy);
      }

      // Build Tooltip HTML
      let qualityInfo = '';
      if (chartId === 'plot_quality') {
        const qVal = closest[1];
        const acc = (100 * (1 - Math.pow(10, -qVal / 10))).toFixed(3);
        const qColor = qVal >= 30 ? '#10b981' : (qVal >= 20 ? '#f59e0b' : '#ef4444');
        const qStatus = qVal >= 30 ? 'High Quality (≥99.9% base accuracy)' : (qVal >= 20 ? 'Acceptable (99% - 99.9%)' : 'Sub-optimal (<99%)');
        qualityInfo = `
          <div style="margin-top:4px;padding-top:4px;border-top:1px dashed #334155;color:${qColor}">
            <b>Accuracy:</b> ${acc}%<br/>
            <span style="font-size:10px">${qStatus}</span>
          </div>
        `;
      }

      const tipHtml = `
        <div style="font-weight:700;color:var(--teal);margin-bottom:3px">${lblX}: ${closest[0].toLocaleString()} ${unitX}</div>
        <div><b>${lblY}:</b> <span style="font-family:monospace;font-weight:700;color:#ffffff">${closest[1].toLocaleString()} ${unitY}</span></div>
        ${qualityInfo}
      `;
      showTooltip(tipHtml, e.clientX, e.clientY);
    });

    card.addEventListener('mouseleave', function() {
      if (crosshair) crosshair.style.display = 'none';
      hideTooltip();
    });
  });

  // 3. ACGT MULTI-LINE COMPOSITION INTERACTION
  document.querySelectorAll('.interactive-plot-card[data-chart-type="acgt"]').forEach(card => {
    const chartId = card.getAttribute('data-chart-id');
    const svg = card.querySelector('svg');
    const xMin = parseFloat(card.getAttribute('data-xmin'));
    const xMax = parseFloat(card.getAttribute('data-xmax'));
    const yMin = parseFloat(card.getAttribute('data-ymin'));
    const yMax = parseFloat(card.getAttribute('data-ymax'));
    const padL = parseFloat(card.getAttribute('data-pad-left'));
    const padR = parseFloat(card.getAttribute('data-pad-right'));
    const padT = parseFloat(card.getAttribute('data-pad-top'));
    const padB = parseFloat(card.getAttribute('data-pad-bottom'));
    const width = parseFloat(card.getAttribute('data-width'));
    const height = parseFloat(card.getAttribute('data-height'));
    const gccData = JSON.parse(card.getAttribute('data-gcc') || '[]');

    const crosshair = document.getElementById(`${chartId}_crosshair`);
    const guideX = document.getElementById(`${chartId}_guide_x`);
    const dotA = document.getElementById(`${chartId}_dot_A`);
    const dotC = document.getElementById(`${chartId}_dot_C`);
    const dotG = document.getElementById(`${chartId}_dot_G`);
    const dotT = document.getElementById(`${chartId}_dot_T`);

    const plotW = width - padL - padR;
    const plotH = height - padT - padB;

    function toSvg(x, y) {
      return {
        x: padL + ((x - xMin) / (xMax - xMin)) * plotW,
        y: padT + plotH - ((y - yMin) / (yMax - yMin)) * plotH
      };
    }

    card.addEventListener('mousemove', function(e) {
      if (!gccData || !gccData.length) return;
      const rect = svg.getBoundingClientRect();
      const mouseX = ((e.clientX - rect.left) / rect.width) * width;
      const mouseY = ((e.clientY - rect.top) / rect.height) * height;

      if (mouseX < padL || mouseX > width - padR || mouseY < padT || mouseY > height - padB) {
        if (crosshair) crosshair.style.display = 'none';
        hideTooltip();
        return;
      }

      const cycleVal = xMin + ((mouseX - padL) / plotW) * (xMax - xMin);
      let closest = gccData[0];
      let minDist = Math.abs(gccData[0].cycle - cycleVal);
      for (let i = 1; i < gccData.length; i++) {
        const d = Math.abs(gccData[i].cycle - cycleVal);
        if (d < minDist) {
          minDist = d;
          closest = gccData[i];
        }
      }

      const posA = toSvg(closest.cycle, closest.A);
      const posC = toSvg(closest.cycle, closest.C);
      const posG = toSvg(closest.cycle, closest.G);
      const posT = toSvg(closest.cycle, closest.T);

      if (crosshair) crosshair.style.display = 'block';
      if (guideX) { guideX.setAttribute('x1', posA.x); guideX.setAttribute('x2', posA.x); }
      if (dotA) { dotA.setAttribute('cx', posA.x); dotA.setAttribute('cy', posA.y); }
      if (dotC) { dotC.setAttribute('cx', posC.x); dotC.setAttribute('cy', posC.y); }
      if (dotG) { dotG.setAttribute('cx', posG.x); dotG.setAttribute('cy', posG.y); }
      if (dotT) { dotT.setAttribute('cx', posT.x); dotT.setAttribute('cy', posT.y); }

      const gcContent = (closest.G + closest.C).toFixed(1);
      const tipHtml = `
        <div style="font-weight:700;color:var(--teal);margin-bottom:4px">Cycle: ${closest.cycle} bp</div>
        <div style="display:grid;grid-template-columns:1fr 1fr;gap:4px;font-family:monospace;font-size:11px">
          <span style="color:#10b981"><b>A:</b> ${closest.A.toFixed(1)}%</span>
          <span style="color:#0284c7"><b>C:</b> ${closest.C.toFixed(1)}%</span>
          <span style="color:#f59e0b"><b>G:</b> ${closest.G.toFixed(1)}%</span>
          <span style="color:#ef4444"><b>T:</b> ${closest.T.toFixed(1)}%</span>
        </div>
        <div style="margin-top:4px;padding-top:4px;border-top:1px dashed #334155;color:#38bdf8">
          <b>GC Content:</b> ${gcContent}%
        </div>
      `;
      showTooltip(tipHtml, e.clientX, e.clientY);
    });

    card.addEventListener('mouseleave', function() {
      if (crosshair) crosshair.style.display = 'none';
      hideTooltip();
    });
  });

  // 4. BAR CHARTS (READ LENGTH, INDELS) INTERACTION
  document.querySelectorAll('.interactive-bar').forEach(bar => {
    bar.addEventListener('mouseenter', function(e) {
      const x = this.getAttribute('data-x');
      const y = parseInt(this.getAttribute('data-y') || 0);
      const pct = this.getAttribute('data-pct');
      const type = this.getAttribute('data-type');
      const sz = this.getAttribute('data-sz');
      const count = parseInt(this.getAttribute('data-count') || 0);
      const lbl = this.getAttribute('data-label') || 'Length';

      let tipHtml = '';
      if (type) {
        const color = type === 'Insertion' ? '#0d9488' : '#ef4444';
        tipHtml = `
          <div style="font-weight:700;color:${color};font-size:12px;margin-bottom:3px">${type} (${sz} bp)</div>
          <div><b>Event Count:</b> <span style="font-family:monospace;font-weight:700">${count.toLocaleString()} events</span></div>
        `;
      } else {
        tipHtml = `
          <div style="font-weight:700;color:var(--teal);margin-bottom:3px">${lbl}: ${x} bp</div>
          <div><b>Read Count:</b> <span style="font-family:monospace;font-weight:700">${y.toLocaleString()} reads</span></div>
          ${pct ? `<div><b>Share:</b> <span style="color:#10b981;font-weight:700">${pct}%</span> of total</div>` : ''}
        `;
      }
      showTooltip(tipHtml, e.clientX, e.clientY);
    });

    bar.addEventListener('mouseleave', hideTooltip);
  });

  // 5. SERIES TOGGLING (ACGT & INDELS)
  function toggleAcgtSeries(chartId, base) {
    const line = document.getElementById(`${chartId}_path_${base}`);
    const btn = document.getElementById(`${chartId}_btn_${base}`);
    const dot = document.getElementById(`${chartId}_dot_${base}`);
    if (!line || !btn) return;

    const isActive = btn.classList.contains('active');
    if (isActive) {
      btn.classList.remove('active');
      btn.classList.add('inactive');
      line.style.display = 'none';
      if (dot) dot.style.display = 'none';
    } else {
      btn.classList.remove('inactive');
      btn.classList.add('active');
      line.style.display = 'block';
      if (dot) dot.style.display = 'block';
    }
  }

  function toggleIndelSeries(chartId, type) {
    const bars = document.querySelectorAll(`#${chartId}_wrapper .indel-bar-${type}`);
    const btn = document.querySelector(`#${chartId}_wrapper .indel-toggle[data-type="${type}"]`);
    if (!bars.length || !btn) return;

    const isActive = btn.classList.contains('active');
    if (isActive) {
      btn.classList.remove('active');
      btn.classList.add('inactive');
      bars.forEach(b => b.style.display = 'none');
    } else {
      btn.classList.remove('inactive');
      btn.classList.add('active');
      bars.forEach(b => b.style.display = 'block');
    }
  }

  // 6. SVG / PNG EXPORT ENGINE
  function exportPlotSvg(svgId, filename) {
    const svg = document.getElementById(svgId);
    if (!svg) return;
    const serializer = new XMLSerializer();
    const source = serializer.serializeToString(svg);
    const blob = new Blob([source], { type: 'image/svg+xml;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = filename || 'plot.svg';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  }

  function exportPlotPng(svgId, filename) {
    const svg = document.getElementById(svgId);
    if (!svg) return;
    const serializer = new XMLSerializer();
    const source = serializer.serializeToString(svg);
    const svgBlob = new Blob([source], { type: 'image/svg+xml;charset=utf-8' });
    const URL_obj = window.URL || window.webkitURL || window;
    const blobURL = URL_obj.createObjectURL(svgBlob);

    const image = new Image();
    image.onload = function() {
      const canvas = document.createElement('canvas');
      const scale = 2.0; // High-res retina scale
      canvas.width = (svg.viewBox.baseVal.width || svg.clientWidth || 540) * scale;
      canvas.height = (svg.viewBox.baseVal.height || svg.clientHeight || 240) * scale;
      const ctx = canvas.getContext('2d');
      ctx.scale(scale, scale);
      
      // Draw background matching theme
      const isDark = document.documentElement.getAttribute('data-theme') === 'dark';
      ctx.fillStyle = isDark ? '#131d31' : '#ffffff';
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      
      ctx.drawImage(image, 0, 0);
      URL_obj.revokeObjectURL(blobURL);

      const pngUrl = canvas.toDataURL('image/png');
      const link = document.createElement('a');
      link.href = pngUrl;
      link.download = filename || 'plot.png';
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    };
    image.src = blobURL;
  }

  function filterContigRows() {
    const q = (document.getElementById('contigSearchInput')?.value || '').trim().toLowerCase();
    const rows = document.querySelectorAll('#contigDistTable tbody tr');
    rows.forEach(r => {
      const txt = r.textContent.toLowerCase();
      r.style.display = (!q || txt.includes(q)) ? '' : 'none';
    });
  }
</script>

</body>
</html>
"""
    # Replace placeholders
    html_content = html_template.replace("__SAMPLE__", sample_esc) \
                                .replace("__DATE_STR__", date_str) \
                                .replace("__ALIGN_MODE_DISP__", html.escape(align_mode_disp)) \
                                .replace("__NUMT_KPI_CARD__", numt_kpi_card_html) \
                                .replace("__NUMT_PLOT_CARD__", numt_plot_card_html) \
                                .replace("__TOTAL_SEQ__", total_seq_disp) \
                                .replace("__PCT_MAPPED__", f"{overall_map_pct:.2f}") \
                                .replace("__PCT_MAPPED_1DEC__", f"{overall_map_pct:.1f}") \
                                .replace("__PCT_DUP__", f"{pct_dup:.2f}") \
                                .replace("__AVG_LEN__", avg_len_disp) \
                                .replace("__AVG_IS__", avg_is_disp) \
                                .replace("__IS_STD__", is_std_disp) \
                                .replace("__ERR_STR__", err_str) \
                                .replace("__SVG_DONUT_ALIGN__", svg_donut_alignment) \
                                .replace("__SVG_DONUT_PAIRS__", svg_donut_orientation) \
                                .replace("__SVG_PLOT_IS__", svg_plot_is) \
                                .replace("__SVG_PLOT_QUALITY__", svg_plot_quality) \
                                .replace("__SVG_PLOT_ACGT__", svg_plot_acgt) \
                                .replace("__SVG_PLOT_RL__", svg_plot_rl) \
                                .replace("__SVG_PLOT_GC__", svg_plot_gc) \
                                .replace("__SVG_PLOT_INDEL__", svg_plot_indel) \
                                .replace("__ALIGN_ROWS__", align_tr_html if align_tr_html else "<tr><td colspan='4' style='color:#888'>No alignment data available</td></tr>") \
                                .replace("__PAIR_ROWS__", pair_tr_html if pair_tr_html else "<tr><td colspan='4' style='color:#888'>No pairing data available</td></tr>") \
                                .replace("__YIELD_ROWS__", yield_tr_html if yield_tr_html else "<tr><td colspan='4' style='color:#888'>No yield data available</td></tr>") \
                                .replace("__IDXSTATS_ROWS__", idx_tr_html if idx_tr_html else "<tr><td colspan='4' style='color:#888'>No contig mapping records</td></tr>")

    os.makedirs(os.path.dirname(os.path.abspath(out_html_path)), exist_ok=True)
    with open(out_html_path, "w", encoding="utf-8") as f:
        f.write(html_content)
    print(f"Generated BAM stats HTML with interactive visual plots: {out_html_path}")
    return out_html_path

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python3 generate_bam_stats_html.py <SAMPLE> <stats.txt> <flagstat.txt> <idxstats.txt> <out.html> [numt_depletion.json]")
        sys.exit(1)
    numt_arg = sys.argv[6] if len(sys.argv) > 6 else None
    generate_bam_stats_html(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], numt_arg)
