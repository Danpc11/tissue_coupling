#!/usr/bin/env python3
"""Figura diagnostica del paso 05: distribuciones por persona, P explicitas."""
import os, sys, importlib.util
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
os.environ.setdefault("T2D_OUT", "results")
HERE = os.path.dirname(os.path.abspath(__file__))
ld = lambda n, f: (lambda sp: (sp.loader.exec_module(m := importlib.util.module_from_spec(sp)) or m))(importlib.util.spec_from_file_location(n, f))
c5 = ld("c5", os.path.join(HERE, "../analyses/05_clamp_within_person.py"))
p4 = ld("p4", os.path.join(HERE, "../analyses/04_perturbations.py"))
E, d = c5.load()
t = pd.read_csv("results/template_2harm.tsv", sep="\t").set_index("gene")
allscore = set(sum(c5.SCORE_SETS.values(), []))
ranked = [x for x in t.sort_values("r2", ascending=False).index if x in E.index]
g800 = ranked[:800]; gclean = [x for x in ranked if x not in allscore][:800]
rng = np.random.default_rng(0)
R = lambda h: float(abs(np.mean(np.exp(1j * np.asarray(h) * 2 * np.pi / 24))))
BLUE, ORANGE, GREY = "#1f77b4", "#d95f02", "#777777"

mus = (d.tissue == "muscle") & (d.tzd == "pre")
pre, post = c5.pairs(d, mus & (d.clamp == "pre"), mus & (d.clamp == "post"))
a, b = list(pre.gsm), list(post.gsm); s = pre.sensitivity.values
# basal: centrado dentro de las basales
base = c5.phase_amplitude_signed(E[a], t, g800, p4).set_index("sample")
# pre+post: centrado conjunto
both = c5.phase_amplitude_signed(E[a + b], t, g800, p4).set_index("sample")
clean = c5.phase_amplitude_signed(E[a + b], t, gclean, p4).set_index("sample")

fig, ax = plt.subplots(2, 3, figsize=(14, 8.5))
plt.rcParams["axes.titlesize"] = 10; ax = ax.ravel()
# a) fase basal por grupo
ms = s == "sensitive"
Rs, Rr = R(base.phase[ms]), R(base.phase[~ms])
null = []
for _ in range(5000):
    pm = rng.permutation(ms); null.append(R(base.phase[pm]) - R(base.phase[~pm]))
pR = (np.sum(np.abs(null) >= abs(Rs - Rr)) + 1) / 5001
keep = (base.fit >= np.percentile(base.fit, 5)) & (base.fit <= np.percentile(base.fit, 95))
null2 = []
for _ in range(5000):
    pm = rng.permutation(ms[keep.values]); null2.append(R(base.phase[keep][pm]) - R(base.phase[keep][~pm]))
obs2 = R(base.phase[keep][ms[keep.values]]) - R(base.phase[keep][~ms[keep.values]])
pR2 = (np.sum(np.abs(null2) >= abs(obs2)) + 1) / 5001
print("filtro 5-95:", R(base.phase[keep][ms[keep.values]]), R(base.phase[keep][~ms[keep.values]]), pR2)
bins = np.arange(0, 25, 2)
ax[0].hist(base.phase[ms], bins=bins, color=BLUE, alpha=.7, label=f"sensibles n={ms.sum()}, R={Rs:.2f}")
ax[0].hist(base.phase[~ms], bins=bins, color=ORANGE, alpha=.5, label=f"resistentes n={(~ms).sum()}, R={Rr:.2f}")
ax[0].set(xlabel="fase inferida (h), musculo basal", ylabel="personas", title=f"a  Fase basal por grupo\nP(perm, ΔR) = {pR:.3f}")
ax[0].legend(fontsize=8, frameon=False)
# b) amplitud basal por persona
for i, (m, col, lab) in enumerate(((ms, BLUE, "sensibles"), (~ms, ORANGE, "resistentes"))):
    y = base.amplitude[m].values
    ax[1].scatter(i + rng.uniform(-.18, .18, len(y)), y, color=col, s=22, alpha=.8)
    ax[1].hlines(np.median(y), i - .3, i + .3, color="k")
pA = stats.mannwhitneyu(base.amplitude[ms], base.amplitude[~ms]).pvalue
ax[1].set(xticks=[0, 1], xticklabels=["sensibles", "resistentes"], ylabel="amplitud (proyeccion)",
          title=f"b  Amplitud basal por persona\nP(Mann-Whitney) = {pA:.3f}")
# c) desplazamiento de fase por persona tras el clamp
dphi = ((both.phase[b].values - both.phase[a].values + 12) % 24) - 12
ax[2].hist(dphi[ms], bins=np.arange(-12, 13, 2), color=BLUE, alpha=.7, label="sensibles")
ax[2].hist(dphi[~ms], bins=np.arange(-12, 13, 2), color=ORANGE, alpha=.5, label="resistentes")
ax[2].axvline(0, color="k", lw=.8)
ax[2].set(xlabel="fase post − pre (h)", ylabel="personas",
          title=f"c  Desplazamiento de fase tras el clamp\nR de los desplazamientos = {R(dphi):.2f}")
ax[2].legend(fontsize=8, frameon=False)
# d) cambio de scores por persona
names = ["Fasting/FOXO3", "Glucocorticoid", "Clock output", "Core oscillator", "Glucose uptake"]
Et = E[a + b]; txt = []
for i, name in enumerate(names):
    sc, _ = c5.set_scores(Et, c5.SCORE_SETS[name]); dS = sc[b].values - sc[a].values
    for m, col, off in ((ms, BLUE, -.17), (~ms, ORANGE, .17)):
        ax[3].scatter(i + off + rng.uniform(-.1, .1, m.sum()), dS[m], color=col, s=14, alpha=.7)
        ax[3].hlines(np.median(dS[m]), i + off - .14, i + off + .14, color="k", lw=1.5)
    txt.append(f"{name.split('/')[0]}: P0={stats.wilcoxon(dS).pvalue:.1e}, Pgrupo={stats.mannwhitneyu(dS[ms], dS[~ms]).pvalue:.2f}")
ax[3].axhline(0, color="k", lw=.8)
ax[3].set(xticks=range(len(names)), xticklabels=["Ayuno/\nFOXO3", "Gluco-\ncorticoide", "Salida\ndel reloj", "Oscilador\ncentral", "Captacion\nglucosa"], ylabel="Δ score (z), post − pre",
          title="d  Respuesta transcripcional al clamp\npor persona")
ax[3].tick_params(axis="x", labelsize=7)
ax[3].text(.01, .02, "\n".join(txt), transform=ax[3].transAxes, fontsize=6.5, va="bottom")
# e) dGC vs |dphi| (plantilla sin genes de score)
sc, _ = c5.set_scores(Et, c5.SCORE_SETS["Glucocorticoid"]); dG = sc[b].values - sc[a].values
adphi = np.abs(((clean.phase[b].values - clean.phase[a].values + 12) % 24) - 12)
r, p = c5.perm_spearman(dG, adphi, rng, 5000)
ax[4].scatter(dG[ms], adphi[ms], color=BLUE, s=22, alpha=.8, label="sensibles")
ax[4].scatter(dG[~ms], adphi[~ms], color=ORANGE, s=22, alpha=.8, label="resistentes")
ax[4].set(xlabel="Δ score glucocorticoide (post − pre)", ylabel="|Δ fase| (h), plantilla sin genes de score",
          title=f"e  Δ eje glucocorticoide frente a |Δ fase|\nρ = {r:.2f}, P(perm) = {p:.3f}")
ax[4].legend(fontsize=8, frameon=False)
# f) sensibilidad a la regla de signo y al centrado
S = p4.waveform(t, g800); GRID = p4.GRID
def est(cols, cc, signed):
    X = E.loc[g800, cols].sub(E.loc[g800, cc].mean(1), axis=0).values
    proj = S @ X; idx = np.argmax(proj if signed else np.abs(proj), axis=0); return GRID[idx] * 24 / (2 * np.pi)
rows = []
for signed, sl in ((False, "|proj| (04)"), (True, "con signo")):
    for cc, cl in ((a, "centrado en basales"), (a + b, "centrado pre+post")):
        ph = est(a + b, cc, signed); P, Q = ph[:len(a)], ph[len(a):]
        rows.append((f"{sl}\n{cl}", R(P[ms]), R(P[~ms]), R(Q)))
x = np.arange(len(rows)); w = .25
ax[5].bar(x - w, [r_[1] for r_ in rows], w, color=BLUE, label="R basal sensibles")
ax[5].bar(x, [r_[2] for r_ in rows], w, color=ORANGE, label="R basal resistentes")
ax[5].bar(x + w, [r_[3] for r_ in rows], w, color=GREY, label="R post-clamp (todos)")
ax[5].set(xticks=x, xticklabels=[r_[0] for r_ in rows], ylabel="concentracion circular R",
          title="f  R segun regla de signo y centrado")
ax[5].tick_params(axis="x", labelsize=7); ax[5].legend(fontsize=7, frameon=False)
for k in ax:
    k.spines[["top", "right"]].set_visible(False); k.title.set_fontsize(10)
fig.suptitle("GSE13070, musculo preTZD: 68 personas con biopsia antes y despues del clamp (plantilla de dos armonicos, 800 genes)", fontsize=10)
fig.tight_layout()
os.makedirs("results/figs", exist_ok=True)
fig.savefig("results/figs/F_clamp_within_person.png", dpi=170); fig.savefig("results/figs/F_clamp_within_person.pdf")
print("ok", pR, pA, R(dphi), r, p)
