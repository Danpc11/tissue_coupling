#!/usr/bin/env python3
"""Figura del contraste clamp frente a salina a las mismas horas (GSE40589).

Lee las tablas que escribe python/analyses/28_clamp_vs_saline.py y dibuja cuatro paneles:
  A  cambio del brazo de salida del reloj en cada par, por brazo
  B  cambio del conjunto de respuesta a insulina (control positivo del ensayo)
  C  cambio de las dianas glucocorticoides
  D  gen a gen: cambio bajo salina frente a cambio bajo clamp euglucemico

Cada punto de A-C es una persona (un par de biopsias, -40 y 240 min). Los pares que no pasan la
regla de identidad del paso 28 se dibujan huecos y no entran en ninguna estadistica. Los valores
de P se escriben; no se usan asteriscos.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from paper_style import use_style, figure, panel_label, OI, WIDTH_DOUBLE
OUT = os.environ.get("T2D_OUT", "results")

ARMS = [("NACL_INF", "Saline", OI["grey"]),
        ("NACL_EU", "Euglycaemic\nclamp", OI["orange"]),
        ("NACL_HC", "Hyperglycaemic\nclamp", "#A8660A")]
PANELS = [("A", "Clock output", "Change in clock-output score\n(240 min − baseline, z)"),
          ("B", "Insulin core", "Change in insulin-response score\n(240 min − baseline, z)"),
          ("C", "Glucocorticoid", "Change in glucocorticoid-target score\n(240 min − baseline, z)")]
SHOW = {"C1orf51": "CIART", "ARNTL": "BMAL1"}          # simbolos vigentes de HGNC


def fmt_p(p):
    return f"P = {p:.3f}" if p >= 0.001 else f"P = {p:.1e}"


def strip_panel(ax, P, letter, gset, ylabel):
    rng = np.random.default_rng(0)
    keep = ~P.excluded.values
    vals = [P.loc[(P.arm == a).values & keep, f"d_{gset}"].values for a, _, _ in ARMS]
    out = [P.loc[(P.arm == a).values & ~keep, f"d_{gset}"].values for a, _, _ in ARMS]
    allv = np.concatenate(vals + out)
    lo, hi = allv.min(), allv.max(); span = hi - lo
    for k, ((a, lab, col), v, o) in enumerate(zip(ARMS, vals, out)):
        ax.scatter(k + rng.uniform(-0.13, 0.13, len(v)), v, s=22, color=col, edgecolor="white", linewidth=0.5, zorder=3)
        if len(o):                                         # excluidos por la regla de identidad
            ax.scatter(k + rng.uniform(-0.13, 0.13, len(o)), o, s=22, facecolor="none", edgecolor=col, linewidth=0.8, zorder=3)
        ax.hlines(np.median(v), k - 0.26, k + 0.26, color="black", lw=1.3, zorder=4)
        ax.text(k, hi + 0.07 * span, f"{int((v < 0).sum())}/{len(v)} down\n{fmt_p(stats.wilcoxon(v).pvalue)}",
                ha="center", va="bottom", fontsize=6.5)
    # clamp frente a salina: una llave por comparacion, apiladas por encima de las anotaciones
    for j, k in enumerate((1, 2)):
        y = hi + (0.36 + 0.17 * j) * span
        p = stats.mannwhitneyu(vals[k], vals[0]).pvalue
        ax.plot([0, 0, k, k], [y - 0.02 * span, y, y, y - 0.02 * span], color="black", lw=0.7)
        ax.text(k / 2, y + 0.01 * span, f"vs saline, {fmt_p(p)}", ha="center", va="bottom", fontsize=6.5)
    ax.axhline(0, color="black", lw=0.6, ls=":")
    ax.set_xticks(range(len(ARMS)))
    ax.set_xticklabels([f"{lab}\n(n = {len(v)})" for (_, lab, _), v in zip(ARMS, vals)])
    ax.set_xlim(-0.6, len(ARMS) - 0.4); ax.set_ylim(lo - 0.08 * span, hi + 0.72 * span)
    ax.set_ylabel(ylabel)
    panel_label(ax, letter, dx=-0.16)


def gene_panel(ax, G, letter):
    cols = {"Clock output": OI["blue"], "Core oscillator": OI["skyblue"]}
    for s, c in cols.items():
        q = G[G.set == s]
        ax.scatter(q.saline_median, q.clamp_median, s=26, color=c, edgecolor="white", linewidth=0.5,
                   zorder=3, label=s.replace("Clock output", "Output arm"))
    x, y = G.saline_median.values, G.clamp_median.values
    lim = [min(x.min(), y.min()) - 0.25, max(x.max(), y.max()) + 0.25]
    ax.plot(lim, lim, color=OI["grey"], lw=0.8, ls=":", zorder=1)
    b = np.polyfit(x, y, 1); xx = np.linspace(x.min(), x.max(), 50)
    ax.plot(xx, np.polyval(b, xx), color="black", lw=1.0, ls="--", zorder=2)
    r = stats.spearmanr(x, y)
    ax.text(0.97, 0.04, f"slope = {b[0]:.2f}\n$\\rho$ = {r.statistic:.2f}, {fmt_p(r.pvalue)}\ndotted: identity",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5)
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
    ax.set_xlabel("Change under saline (log$_2$)"); ax.set_ylabel("Change under euglycaemic clamp (log$_2$)")
    ax.legend(frameon=False, loc="upper left", fontsize=6.5, handletextpad=0.2, borderaxespad=0.2)
    panel_label(ax, letter, dx=-0.16)
    place_labels(ax, x, y, [SHOW.get(g, g) for g in G.gene])


def place_labels(ax, x, y, names):
    """Coloca cada etiqueta en la primera posicion libre junto a su punto. Prueba primero a los
    lados, luego arriba y abajo y luego en diagonal, a distancias crecientes, y descarta las
    posiciones que pisan otra etiqueta u otro punto. Si tiene que alejarse traza una linea fina."""
    fig = ax.figure; fig.canvas.draw(); rend = fig.canvas.get_renderer()
    box = ax.get_window_extent(rend)
    pts = ax.transData.transform(np.column_stack([x, y]))
    half = 3.2 * fig.dpi / 72                                           # radio del marcador, en pixeles
    taken = [matplotlib.transforms.Bbox.from_extents(px - half, py - half, px + half, py + half) for px, py in pts]
    dirs = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, 1), (1, -1), (-1, -1)]
    for i in np.argsort(-np.hypot(x - x.mean(), y - y.mean())):        # primero los aislados
        best = None
        for rad in (4.5, 9, 14, 20):
            for dx, dy in dirs:
                k = rad * (0.75 if dx and dy else 1)
                t = ax.annotate(names[i], (x[i], y[i]), xytext=(dx * k, dy * k), textcoords="offset points",
                                ha="left" if dx > 0 else "right" if dx < 0 else "center",
                                va="bottom" if dy > 0 else "top" if dy < 0 else "center",
                                fontsize=6, style="italic", zorder=5,
                                bbox=dict(fc="white", ec="none", pad=0.3, alpha=0.85))
                bb = t.get_window_extent(rend).expanded(1.04, 1.10)
                ok = box.contains(bb.x0, bb.y0) and box.contains(bb.x1, bb.y1) \
                    and not any(bb.overlaps(o) for j, o in enumerate(taken) if j != i)
                if ok:
                    best = (t, bb, rad); break
                t.remove()
            if best: break
        if best is None:                                                # sin hueco: posicion por defecto
            t = ax.annotate(names[i], (x[i], y[i]), xytext=(4, 4), textcoords="offset points",
                            fontsize=6, style="italic", zorder=5)
            best = (t, t.get_window_extent(rend), 4.5)
        t, bb, rad = best; taken.append(bb)
        if rad > 4.5:
            ax.annotate("", (x[i], y[i]), xytext=t.xyann, textcoords="offset points", zorder=2,
                        arrowprops=dict(arrowstyle="-", lw=0.4, color=OI["grey"], shrinkA=0, shrinkB=2))


if __name__ == "__main__":
    use_style()
    P = pd.read_csv(f"{OUT}/new/clamp_vs_saline_pairs.tsv", sep="\t")
    if "excluded" not in P.columns: P["excluded"] = False
    G = pd.read_csv(f"{OUT}/new/clamp_vs_saline_genes.tsv", sep="\t")
    fig = figure(width=WIDTH_DOUBLE, height=WIDTH_DOUBLE * 0.86)
    ax = fig.subplots(2, 2)
    for (letter, gset, ylabel), a in zip(PANELS, (ax[0, 0], ax[0, 1], ax[1, 0])):
        strip_panel(a, P, letter, gset, ylabel)
    gene_panel(ax[1, 1], G, "D")
    os.makedirs(f"{OUT}/figs", exist_ok=True)
    fig.savefig(f"{OUT}/figs/clamp_vs_saline.pdf"); fig.savefig(f"{OUT}/figs/clamp_vs_saline.png", dpi=400)
    print(f"ok -> {OUT}/figs/clamp_vs_saline.pdf")
