#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gillespie SSA 反应模拟器（原型）
输入：反应表 CSV + 命令行初始分子数
核心：倾向 = k × Π C(n_i, |ν_i|)（组合数），加权随机选反应，
      物理时间 Δt = -ln(u)/a0，滑动窗口均值判断稳态终止。
"""
import argparse, csv, math, random, sys, time
from math import comb

def parse_stoich(s):
    return int(s)

def load_reactions(path):
    """CSV 表头：反应ID,反应名称,反应物A,反应物B,反应物C,反应物D,反应速率
       （反应物X 列实际存化学计量系数：反应物为负、产物为正）"""
    reactions, species = [], []
    with open(path, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            r = {'id': (row.get('反应ID') or row.get('id') or row.get('ID') or '').strip(),
                 'name': (row.get('反应名称') or row.get('name') or '').strip(),
                 'k': float(row.get('反应速率') or row.get('rate') or row.get('k'))}
            # 找到所有化学计量列（反应物A/B/C/D 或通用 stoich_*）
            st = {}
            reserved = {'反应ID', 'id', 'ID', '反应名称', 'name', '反应速率', 'rate', 'k'}
            for k, v in row.items():
                if v.strip() == '' or k in reserved:
                    continue
                for pre in ('反应物', 'stoich_'):
                    if k.startswith(pre):
                        k = k[len(pre):]
                        break
                st[k] = parse_stoich(v)
            r['stoich'] = st
            reactions.append(r)
            for sp in st:
                if sp not in species:
                    species.append(sp)
    for sp in species:  # 每个反应补齐所有物种的系数
        for r in reactions:
            r['stoich'].setdefault(sp, 0)
    return reactions, species

def propensity(r, counts):
    """a = k × Π C(n_i, |ν_i|) 只对 ν<0（反应物）；不足时组合数为0 → 倾向0 → 不会被选中"""
    a = r['k']
    for sp, nu in r['stoich'].items():
        if nu < 0:
            n = counts.get(sp, 0)
            need = -nu
            if n < need:
                return 0.0
            a *= comb(n, need)
    return a

def simulate(reactions, species, init, max_steps, max_time,
             window, rel_eps, seed):
    rng = random.Random(seed)
    counts = {sp: init.get(sp, 0) for sp in species}
    t, step = 0.0, 0
    traj = []          # 详细轨迹记录
    prop_cache = [propensity(r, counts) for r in reactions]
    window_sums = {sp: 0.0 for sp in species}
    prev_means, terminate = None, None

    while True:
        if step >= max_steps:
            terminate = f'达到最大步数 {max_steps}'
            break
        a0 = sum(prop_cache)
        if a0 <= 0:
            terminate = '总倾向为 0（反应物耗尽，无反应可发生）'
            break
        if t >= max_time:
            terminate = f'达到最大时间 {max_time}'
            break

        # 1) 按倾向加权选一个反应
        u, acc = rng.random() * a0, 0.0
        for j, a in enumerate(prop_cache):
            acc += a
            if u <= acc:
                chosen = j
                break
        r = reactions[chosen]

        # 2) 推进物理时间（事件间隔 ~ Exp(a0)）
        t += -math.log(max(rng.random(), 1e-300)) / a0

        # 3) 执行反应：更新分子数，只重算受影响反应的倾向
        affected = {sp for sp, nu in r['stoich'].items() if nu != 0}
        old = dict(counts)
        for sp, nu in r['stoich'].items():
            if nu:
                counts[sp] = counts.get(sp, 0) + nu
        for j, rj in enumerate(reactions):
            if any(rj['stoich'].get(sp, 0) != 0 for sp in affected):
                prop_cache[j] = propensity(rj, counts)

        step += 1
        traj.append({'step': step, 'time': t, 'reaction': r['id'],
                     'reaction_name': r['name'], **counts})

        # 4) 稳态检测：滑动窗口均值变化
        if rel_eps is not None and step >= window:
            for sp in species:
                window_sums[sp] += counts[sp]
            if step % window == 0:
                means = {sp: window_sums[sp] / window for sp in species}
                window_sums = {sp: 0.0 for sp in species}
                if prev_means is not None:
                    stable = all(abs(means[sp] - prev_means[sp]) <=
                                 rel_eps * max(1.0, abs(prev_means[sp]))
                                 for sp in species)
                    if stable:
                        terminate = (f'稳态：连续窗口均值相对变化 ≤ {rel_eps} '
                                     f'（窗口 {window}）')
                        break
                prev_means = means

    return counts, t, step, traj, terminate

def main():
    p = argparse.ArgumentParser(description='Gillespie SSA 反应模拟器')
    p.add_argument('--csv', required=True, help='反应表 CSV 路径')
    p.add_argument('--init', required=True,
                   help='初始分子数，如 "A=100,B=200,C=0,D=0"')
    p.add_argument('--max-steps', type=int, default=100000, help='最大事件步数')
    p.add_argument('--max-time', type=float, default=1e9, help='最大模拟时间')
    p.add_argument('--window', type=int, default=200, help='稳态检测滑动窗口大小')
    p.add_argument('--rel-eps', type=float, default=0.01,
                   help='窗口均值相对变化阈值（设 0 关闭稳态检测）')
    p.add_argument('--seed', type=int, default=None)
    p.add_argument('--out', default='trajectory.csv', help='轨迹输出 CSV')
    a = p.parse_args()

    init = {}
    for item in a.init.split(','):
        sp, n = item.split('=')
        init[sp.strip()] = int(n)

    reactions, species = load_reactions(a.csv)
    print(f'载入 {len(reactions)} 个反应，物种：{species}')
    print(f'初始：{init}')

    t0 = time.time()
    final, t_end, steps, traj, reason = simulate(
        reactions, species, init, a.max_steps, a.max_time,
        a.window, a.rel_eps or None, a.seed)
    wall = time.time() - t0

    with open(a.out, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['step', 'time', 'reaction',
                                          'reaction_name'] + species)
        w.writeheader()
        w.writerows(traj)

    print(f'\n=== 模拟结束（用时 {wall:.2f}s 真实时间）===')
    print(f'终止原因：{reason}')
    print(f'事件步数：{steps}，模拟时间：{t_end:.4g}')
    print('分子数变化：')
    for sp in species:
        print(f'  {sp}: {init.get(sp,0)} -> {final.get(sp,0)}')
    print(f'详细轨迹已写入 {a.out}（{len(traj)} 条记录）')

if __name__ == '__main__':
    main()
