#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ssa_plot.py — SSA 轨迹 CSV 的通用可视化（Origin 风格）
用法示例：
  python ssa_plot.py trajectory.csv
  python ssa_plot.py trajectory.csv --x step --ylog --smooth 50 --hist
  python ssa_plot.py trajectory.csv --species NO,O2 --out fig.tiff --dpi 300
"""
import argparse
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

META_COLS = {'step', 'time', 'reaction', 'reaction_name', 'dt',
             'N2O2_eq', 'shot', 'slow', 'r_theory', 'r_actual', 'eq_dev'}

# Origin 默认分组配色
ORIGIN_COLORS = ['#000000', '#E8000B', '#023EFF', '#00A08B',
                 '#FF7F0E', '#8C0A90', '#8B5E00', '#0072B2',
                 '#D55E00', '#009E73']
MARKERS = ['s', 'o', '^', 'D', 'v', 'p', '<', 'h', '>', 'x']


def apply_origin_style():
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 10.5,
        'mathtext.fontset': 'stix',
        'axes.linewidth': 1.4,
        'axes.labelsize': 12,
        'axes.formatter.use_mathtext': True,
        'xtick.direction': 'in', 'ytick.direction': 'in',
        'xtick.top': True, 'ytick.right': True,
        'xtick.major.size': 5.0, 'ytick.major.size': 5.0,
        'xtick.minor.size': 2.8, 'ytick.minor.size': 2.8,
        'xtick.major.width': 1.4, 'ytick.major.width': 1.4,
        'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0,
        'xtick.labelsize': 10.5, 'ytick.labelsize': 10.5,
        'legend.frameon': True, 'legend.fancybox': False,
        'legend.edgecolor': 'black', 'legend.framealpha': 1.0,
        'legend.fontsize': 10,
        'figure.dpi': 110, 'savefig.dpi': 300,
        'savefig.bbox': 'tight', 'savefig.pad_inches': 0.05,
    })


def smooth(y, n):
    if n <= 1:
        return y
    k = np.ones(n) / n
    ypad = np.pad(y, (n // 2, n - 1 - n // 2), mode='edge')
    return np.convolve(ypad, k, mode='valid')


def main():
    p = argparse.ArgumentParser(description='SSA 轨迹可视化（Origin 风格）')
    p.add_argument('csv', help='轨迹 CSV（ssa_simulator.py 输出）')
    p.add_argument('--x', choices=['time', 'step'], default='time',
                   help='横轴：模拟时间或事件步（默认 time）')
    p.add_argument('--species', default=None,
                   help='只画这些物种（逗号分隔，默认全部）')
    p.add_argument('--ylog', action='store_true', help='纵轴对数')
    p.add_argument('--smooth', type=int, default=1,
                   help='滑动平均窗口（默认 1 不平滑）')
    p.add_argument('--style', choices=['line', 'both'], default='line',
                   help='line=仅线；both=线+符号（配 --marker-every）')
    p.add_argument('--marker-every', type=int, default=1,
                   help='both 模式下每隔 N 个点画一个符号')
    p.add_argument('--hist', action='store_true',
                   help='附加子图：各反应触发次数柱状图（横向条形）')
    p.add_argument('--hist-top', type=int, default=12,
                   help='柱状图只显示触发次数前 N 的反应，其余合并为 others')
    p.add_argument('--figsize', default='120,80', help='图尺寸 宽,高（毫米）')
    p.add_argument('--dpi', type=int, default=300)
    p.add_argument('--xlabel', default=None)
    p.add_argument('--ylabel', default='Molecules')
    p.add_argument('--title', default=None)
    p.add_argument('--out', default=None,
                   help='输出文件，支持 png/pdf/tiff/svg（默认 轨迹名_plot.png）')
    a = p.parse_args()

    df = pd.read_csv(a.csv)
    xcol = 'time' if a.x == 'time' else 'step'
    species = [c for c in df.columns if c not in META_COLS]
    if a.species:
        species = [s for s in a.species.split(',') if s in df.columns]
    if not species:
        raise SystemExit('未找到可绘制的物种列')

    w, h = (float(v) for v in a.figsize.split(','))
    n_panel = 2 if a.hist else 1
    fig, axes = plt.subplots(
        1, n_panel, figsize=(w / 25.4, h / 25.4),
        gridspec_kw={'width_ratios': [3, 1]} if a.hist else None)
    if n_panel == 1:
        axes = [axes]
    ax = axes[0]

    x = df[xcol].to_numpy()
    for i, sp in enumerate(species):
        y = smooth(df[sp].to_numpy(dtype=float), a.smooth)
        c = ORIGIN_COLORS[i % len(ORIGIN_COLORS)]
        if a.style == 'both':
            ax.plot(x, y, '-', color=c, lw=1.4, markevery=a.marker_every,
                    marker=MARKERS[i % len(MARKERS)], ms=4.5,
                    mfc='white' if i % 2 else c, mec=c, mew=1.0, label=sp)
        else:
            ax.plot(x, y, '-', color=c, lw=1.6, label=sp)

    ax.set_xlabel(a.xlabel if a.xlabel else
                  ('Time' if a.x == 'time' else 'Event step'))
    ax.set_ylabel(a.ylabel)
    if a.ylog:
        ax.set_yscale('log')
    if a.title:
        ax.set_title(a.title, fontsize=12)
    ax.legend(loc='best', handlelength=2.2, borderpad=0.6, labelspacing=0.4)
    ax.minorticks_on()

    if a.hist:
        src_col = df['reaction_name'] if 'reaction_name' in df.columns else df['reaction']
        cnt = src_col.value_counts()
        if len(cnt) > a.hist_top:
            n_rest = len(cnt) - a.hist_top
            rest = cnt.iloc[a.hist_top:].sum()
            cnt = cnt.iloc[:a.hist_top].copy()
            cnt[f'others ({n_rest} rxns)'] = rest
        cnt = cnt.iloc[::-1]                      # 最大者在顶部
        ypos = np.arange(len(cnt))
        axh = axes[1]
        axh.barh(ypos, cnt.values, height=0.65,
                 color=[ORIGIN_COLORS[i % len(ORIGIN_COLORS)]
                        for i in range(len(cnt))],
                 edgecolor='black', linewidth=0.8)
        axh.set_yticks(ypos)
        axh.set_yticklabels(cnt.index, fontsize=7.5)
        axh.set_xlabel('Events')
        axh.minorticks_on()
        axh.tick_params(axis='y', length=0)

    fig.tight_layout()
    out = a.out or a.csv.rsplit('.', 1)[0] + '_plot.png'
    fig.savefig(out, dpi=a.dpi)
    print(f'已保存 {out}（{w:.0f}x{h:.0f} mm, {a.dpi} dpi）')


if __name__ == '__main__':
    apply_origin_style()
    main()
