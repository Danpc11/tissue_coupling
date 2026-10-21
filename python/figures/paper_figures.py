#!/usr/bin/env python3
"""The five manuscript figures (Nature style, English labels, uppercase panel letters).

F1 The problem and the tool   F2 The axis   F3 Cortisol yes, insulin no
F4 Regulation and specificity F5 The model

Schematic panels are left as labelled empty frames: they are drawn by hand, because
AI-generated images are not accepted by the publishers.
"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../lib")); sys.path.insert(0, os.path.join(HERE, "../model"))
from gene_sets import SETS, resolve
OUT = os.environ.get("T2D_OUT", "results"); RAW = os.environ.get("T2D_RAW", "data/raw")
AD, MU, IS, GC, GR = "#2a9d8f", "#d95f02", "#6a3d9a", "#1f77b4", "#9a9a9a"
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7,
                     "ytick.labelsize": 7, "legend.fontsize": 6.5, "axes.titlesize": 8.5})

def panel(ax, letter, _title=None):
    """Solo la letra del panel. La descripcion va en el pie de figura (FIGURE_CAPTIONS.md),
    que es la convencion de Nature: los paneles no llevan titulo."""
    ax.text(-0.17, 1.13, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")

def frame(ax, letter, txt):
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.03, 0.04), 0.94, 0.86, fill=False, ls="--", lw=.9,
                               ec="#bbbbbb", transform=ax.transAxes))
    ax.text(.5, .5, txt, ha="center", va="center", fontsize=7.5, color="#888888",
            transform=ax.transAxes, style="italic")
    ax.text(-0.05, 1.06, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")

def src(ax, txt, y=1.03, x=0.0):
    """Procedencia del panel: cohorte, tejido, n y metrica, dentro del area del panel."""
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=6.3, color="#666666",
            style="italic", va="bottom", ha="left")

def trend(ax, x, y, color="k"):
    """Línea de tendencia por mínimos cuadrados, con la rho de Spearman ya calculada aparte."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3: return
    b = np.polyfit(x[ok], y[ok], 1); xs = np.linspace(x[ok].min(), x[ok].max(), 50)
    ax.plot(xs, np.polyval(b, xs), color=color, lw=1.1, ls="--", zorder=1)

def sc(X, g):
    gg = [r for r in (resolve(x, X.index) for x in g) if r]
    Z = X.loc[gg]; return ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T.mean(0)

def center(X): return X.sub(X.mean(0), axis=1)

def tidy(fig):
    for a in fig.axes: a.spines[["top", "right"]].set_visible(False)

# ------------------------------------------------------------------ F1
def fig1():
    fig = plt.figure(figsize=(17.5, 14)); gs = fig.add_gridspec(3, 4, hspace=.78, wspace=.52)
    D = pd.read_csv(f"{OUT}/fig1_de_per_cohort.tsv", sep="\t")
    S18 = pd.read_csv(f"{OUT}/fig1_signature_genes.tsv", sep="\t")
    frame(fig.add_subplot(gs[0, 0]), "A",
          "Study design\n\nfive cohorts, three tissues,\nwhat each design can answer\n\n(hand-drawn)")
    # B  classic differential expression, by tissue and cohort
    ax = fig.add_subplot(gs[0, 1])
    d = D.copy(); d["tis"] = d.tejido.map({"músculo": "muscle", "adiposo": "adipose", "islotes": "islets"})
    isl = f"{OUT}/new/GSE159984_islets_T2D_vs_ND.tsv"
    for ext in ("", ".gz"):
        if os.path.exists(isl + ext):
            r = pd.read_csv(isl + ext, sep="\t", index_col=0)
            d = pd.concat([d, pd.DataFrame([dict(cohorte="GSE159984", tejido="islotes", tis="islets",
                 casos=27, controles=58, genes_FDR01=int((r.padj < 0.1).sum()))])], ignore_index=True); break
    d = d.sort_values(["tis", "genes_FDR01"])
    col = {"muscle": MU, "adipose": AD, "islets": IS}
    ax.barh(range(len(d)), d.genes_FDR01.clip(lower=.25), color=[col[t] for t in d.tis])
    for i, (n, c, co) in enumerate(zip(d.genes_FDR01, d.casos, d.controles)):
        ax.text(max(n, .25) * 1.25, i, f"{int(n)}   ({int(c)} vs {int(co)})", va="center", fontsize=6.5)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([c.replace(" · músculo", "").replace(" · adiposo", "").replace(" · islotes", "")
                        for c in d.cohorte], fontsize=7)
    ax.set_xscale("log"); ax.set_xlim(.2, 4000); ax.set_xlabel("genes at FDR < 0.1")
    h = [plt.Line2D([], [], color=col[t], lw=5) for t in ("muscle", "adipose", "islets")]
    ax.legend(h, ["muscle", "adipose", "islets"], frameon=False, loc="lower right")
    src(ax, "case vs control within each cohort; moderated t, Benjamini-Hochberg")
    panel(ax, "B", "Classic differential expression\nby tissue and cohort")
    # C  dot plot of the 18 signature genes
    ax = fig.add_subplot(gs[0, 2:4])
    cohs = list(d.cohorte); genes = [g for g in SETS["Clock output"]["genes"] + SETS["Glucocorticoid"]["genes"]
                                     if g in set(S18.gen)]
    for i, c in enumerate(cohs):
        q = S18[S18.cohorte == c].set_index("gen")
        for j, g in enumerate(genes):
            if g in q.index:
                ax.scatter(j, i, s=10 + 80 * min(-np.log10(max(q.loc[g, "P"], 1e-6)) / 3, 1),
                           c=[q.loc[g, "t"]], cmap="RdBu_r", vmin=-4, vmax=4, edgecolors="k", linewidths=.3)
    ax.axvline(len(SETS["Clock output"]["genes"]) - .5, color="k", lw=.8, ls="--")
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=90, fontsize=6.5)
    ax.set_yticks(range(len(cohs)))
    ax.set_yticklabels([c.replace(" · ", "\n").replace("músculo", "muscle").replace("adiposo", "adipose")
                        for c in cohs], fontsize=6.5)
    sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-4, 4))
    cb = plt.colorbar(sm, ax=ax, shrink=.8, pad=.02); cb.set_label("moderated t (case − control)", fontsize=7)
    ax.set_ylim(-1.5, len(cohs) - .4)
    ax.text(1.5, -1.0, "clock output", fontsize=7, ha="center", style="italic")
    ax.text(len(genes) - 4.5, -1.0, "glucocorticoid", fontsize=7, ha="center", style="italic")
    src(ax, "five cohorts \u00b7 muscle, adipose and islets")
    panel(ax, "C", "No signature gene reaches FDR < 0.1 in any cohort   (dot size: −log$_{10}$ P, uncorrected)")
    # D  set-level enrichment
    Z = pd.read_csv(f"{OUT}/fig1_set_enrichment.tsv", sep="\t")
    p = Z.pivot_table(index="cohorte", columns="conjunto", values="z").reindex(
        columns=["Clock output", "Glucocorticoid", "Insulin core", "Oxidative"]).astype(float)
    ax = fig.add_subplot(gs[1, 0]); im = ax.imshow(p.values, cmap="RdBu_r", vmin=-7, vmax=7, aspect="auto")
    ax.set_xticks(range(4)); ax.set_xticklabels(["clock\noutput", "gluco-\ncorticoid", "insulin", "oxidative"], fontsize=6.5)
    ax.set_yticks(range(len(p)))
    ax.set_yticklabels([c.replace(" · ", "\n").replace("músculo", "muscle").replace("adiposo", "adipose")
                        for c in p.index], fontsize=6.5)
    for i in range(p.shape[0]):
        for j in range(p.shape[1]):
            v = p.values[i, j]
            ax.text(j, i, "—" if np.isnan(v) else f"{v:.1f}", ha="center", va="center", fontsize=6.5,
                    color="w" if (not np.isnan(v) and abs(v) > 4) else "k")
    cb = plt.colorbar(im, ax=ax, shrink=.75, pad=.03); cb.set_label("z", fontsize=7)
    src(ax, "set-level z against 2,000 size-matched random sets")
    panel(ax, "D", "The clock set is already altered;\nthe glucocorticoid set never is")
    frame(fig.add_subplot(gs[1, 1]), "E",
          "Single-sample estimator\n\none biopsy → phase and amplitude,\nvalidated against the time series\nof the same person\n\n(hand-drawn)")
    # F  myotubes
    ax = fig.add_subplot(gs[1, 2])
    M = pd.read_csv(f"{OUT}/new/GSE182117_inferred_free_vs_measured.tsv", sep="\t")
    for trt, c in (("CTL", GC), ("HGI", MU)):
        q = M[M.trt == trt]; r = stats.spearmanr(q.amp_free, q.amp_mean_log2)
        ax.scatter(q.amp_mean_log2, q.amp_free, s=22, color=c, label=f"{trt}   ρ = {r.statistic:.2f}")
        trend(ax, q.amp_mean_log2, q.amp_free, c)
    ax.set(xlabel="measured amplitude, cosinor (log$_2$)", ylabel="inferred amplitude")
    ax.legend(frameon=False, loc="upper left")
    src(ax, "GSE182117 \u00b7 human myotubes \u00b7 11 donors \u00d7 2 conditions \u00b7 8 time points")
    panel(ax, "F", "Synchronised myotubes:\ninferred amplitude is the measured one")
    # G  in vivo
    ax = fig.add_subplot(gs[1, 3])
    G = pd.read_csv(f"{OUT}/new/GSE129843_amplitude_per_person.tsv", sep="\t")
    r = stats.spearmanr(G.amp_measured, G.amp_free)
    for cond, c in (("EXF", GC), ("TRF", AD)):
        q = G[G.cond == cond]; ax.scatter(q.amp_measured, q.amp_free, s=26, color=c, label=cond)
    trend(ax, G.amp_measured, G.amp_free)
    ax.set(xlabel="measured amplitude, cosinor (6 points)", ylabel="inferred amplitude")
    ax.legend(frameon=False, loc="upper left", title=f"ρ = {r.statistic:.2f}, P = {r.pvalue:.0e}",
              title_fontsize=7)
    src(ax, "GSE129843 \u00b7 skeletal muscle in vivo \u00b7 11 men \u00b7 6 biopsies per person and condition")
    panel(ax, "G", f"Skeletal muscle in vivo:\nρ = {r.statistic:.2f}, P = {r.pvalue:.1e}, 22 series")
    # H  per-person distributions
    ax = fig.add_subplot(gs[2, 0]); rng = np.random.default_rng(0)
    COH = []
    X = pd.read_pickle(f"{OUT}/new/GSE25462_entrez.pkl"); md = pd.read_pickle(f"{OUT}/new/GSE25462_meta.pkl").set_index("gsm").reindex(X.columns)
    k = md.index[md.fh.isin(["DM", "Family history negative"])]; COH.append(("GSE25462", X[list(k)], (md.fh.reindex(k) == "DM").values, MU))
    X = pd.read_pickle(f"{OUT}/new/GSE26637_entrez.pkl"); md = pd.read_pickle(f"{OUT}/new/GSE26637_meta.pkl").set_index("gsm").reindex(X.columns)
    k = md.index[md.stim == "fasting"]; COH.append(("GSE26637", X[list(k)], (md.status.reindex(k) == "resistant").values, AD))
    X = pd.read_pickle(f"{OUT}/new/GSE78721_log2.pkl"); md = pd.read_pickle(f"{OUT}/new/GSE78721_meta.pkl").set_index("gsm").reindex(X.columns)
    COH.append(("GSE78721", X, (md.dx != "Normal").values, AD))
    for i, (n, X, grp, c) in enumerate(COH):
        s = np.asarray(sc(center(X), SETS["Clock output"]["genes"]), float)
        for k2, m in enumerate((~grp, grp)):
            v = s[m]; xx = i * 2.3 + k2 * .8
            ax.scatter(np.full(v.shape, xx) + rng.uniform(-.16, .16, len(v)), v, s=8, alpha=.6,
                       color=c, marker="o" if k2 else "s")
            ax.hlines(np.median(v), xx - .26, xx + .26, color="k", lw=1.5)
        ax.text(i * 2.3 + .4, 2.45, f"P = {stats.mannwhitneyu(s[grp], s[~grp]).pvalue:.3f}", fontsize=6.5, ha="center")
    ax.set_xticks([i * 2.3 + .4 for i in range(3)])
    ax.set_xticklabels([f"{c[0]}\n{'muscle' if c[3]==MU else 'adipose'}" for c in COH], fontsize=6.5)
    ax.set(ylabel="clock-output score (z)", ylim=(-2.7, 2.9))
    ax.legend([plt.Line2D([], [], marker="s", ls="", color="k"), plt.Line2D([], [], marker="o", ls="", color="k")],
              ["control", "case"], frameon=False, loc="lower right")
    src(ax, "clock-output score per person, sample-centred expression")
    panel(ax, "H", "Per person: the distributions\noverlap almost completely")
    # I  P-value distribution
    ax = fig.add_subplot(gs[2, 1])
    for c, col_, lab in (("GSE25462_músculo", MU, "GSE25462 muscle"), ("GSE78721_adiposo", AD, "GSE78721 adipose"),
                         ("GSE26637_adiposo", IS, "GSE26637 adipose")):
        for ext in (".tsv", ".tsv.gz"):
            f = f"{OUT}/fig1_de_{c}{ext}"
            if os.path.exists(f):
                r = pd.read_csv(f, sep="\t", index_col=0)
                ax.hist(r.P.dropna(), bins=40, histtype="step", lw=1.3, color=col_, density=True, label=lab); break
    ax.axhline(1, color="k", lw=.8, ls=":")
    ax.set(xlabel="uncorrected P, all genes", ylabel="density")
    ax.legend(frameon=False, loc="upper right")
    src(ax, "all genes, uncorrected P from the moderated t")
    panel(ax, "I", "Only the largest cohort departs from a\nuniform P distribution, and still yields 3 genes")
    # J  rho between individuals
    ax = fig.add_subplot(gs[2, 2])
    R = pd.read_csv(f"{OUT}/new/meta_per_cohort.tsv", sep="\t")
    s2 = R[(R.contraste == "NR3C1~reloj") & (R.cohorte.isin(["GSE25462", "GSE182120", "GSE13070"]))].sort_values("rho")
    lo = np.tanh(s2.z - 1.96 * s2.se); hi = np.tanh(s2.z + 1.96 * s2.se)
    for i, (_, r) in enumerate(s2.iterrows()):
        c = AD if "adiposo" in r.tejido else MU
        ax.plot([lo.iloc[i], hi.iloc[i]], [i, i], color=c, lw=1.4)
        ax.scatter(r.rho, i, s=18 + r.n / 4, color=c, zorder=3)
    ax.axvline(0, color="k", lw=.8); ax.set_yticks(range(len(s2)))
    ax.set_yticklabels([f"{r.cohorte} · {r.tejido.replace('m\u00fasculo','muscle').replace('adiposo','adipose')}  (n = {int(r.n)})"
                        for _, r in s2.iterrows()], fontsize=6.5)
    ax.set_xlabel("\u03c1  receptor ~ clock output"); ax.set_xlim(-.05, .8); ax.set_ylim(-.6, len(s2) - .4)
    src(ax, "NR3C1 vs clock-output score, across individuals; 95% CI")
    panel(ax, "J", "The missing signal lies between\nindividuals, not between groups")
    frame(fig.add_subplot(gs[2, 3]), "K",
          "From comparing groups\nto correlating within individuals\n\n(hand-drawn)")
    tidy(fig)
    fig.savefig(f"{OUT}/figs/F1.png", dpi=400, bbox_inches="tight")
    fig.savefig(f"{OUT}/figs/F1.pdf", bbox_inches="tight"); plt.close(fig)

# ------------------------------------------------------------------ F3
def fig3():
    fig = plt.figure(figsize=(14.5, 14)); gs = fig.add_gridspec(3, 3, hspace=.72, wspace=.46)
    frame(fig.add_subplot(gs[0, 0]), "A",
          "The two perturbations\n\nclamp: 0 h and 3 h, same morning\nhydrocortisone: 07:00 vs 07:00,\ncrossed with placebo\n\n(hand-drawn)")
    C = pd.read_csv(f"{OUT}/new/time_calibration_vs_clamp.tsv", sep="\t", index_col=0)
    gl = ["NR1D1", "DBP", "CIART", "PER1", "PER3", "TEF"]; x = np.arange(len(gl)); w = .2
    ax = fig.add_subplot(gs[0, 1])
    labs = ["GSE129843 EXF", "GSE129843 TRF", "GSE108539", "GSE13070 clamp"]
    for k, (c, col_) in enumerate(zip(C.columns, [GR, AD, "#b5a642", MU])):
        ax.bar(x + (k - 1.5) * w, C.loc[gl, c], w, color=col_, label=labs[k])
    ax.axhline(0, color="k", lw=.8); ax.set_xticks(x); ax.set_xticklabels(gl, rotation=45, fontsize=7)
    ax.set_ylabel("change (log$_2$)"); ax.legend(frameon=False, loc="lower right")
    src(ax, "clock genes \u00b7 clamp vs cohorts with recorded times", y=1.06)
    panel(ax, "B", "Three hours of clamp equal\nthree hours of morning")
    D = pd.read_csv(f"{OUT}/new/GSE26637_clamp_adipose_per_person.tsv", sep="\t", index_col=0) \
        if os.path.exists(f"{OUT}/new/GSE26637_clamp_adipose_per_person.tsv") \
        else pd.read_csv(f"{OUT}/new/GSE26637_clamp_adipose_per_person.tsv.gz", sep="\t", index_col=0)
    ax = fig.add_subplot(gs[0, 2])
    gl2 = ["PDK4", "TXNIP", "DDIT4", "PER1", "NR1D2", "CIART", "PER3", "HK2"]
    for i, g in enumerate(gl2):
        if g in D.index:
            v = D.loc[g].values
            ax.scatter(np.full(v.shape, i) + np.random.default_rng(i).uniform(-.15, .15, len(v)), v, s=11, color=AD, alpha=.75)
            ax.hlines(np.median(v), i - .28, i + .28, color="k", lw=1.5)
    ax.axhline(0, color="k", lw=.8); ax.set_xticks(range(len(gl2))); ax.set_xticklabels(gl2, rotation=45, fontsize=7)
    ax.set_ylabel("hyperinsulinaemia − fasting (log$_2$)")
    src(ax, "GSE26637 \u00b7 adipose \u00b7 10 people, paired", y=1.06)
    panel(ax, "C", "Adipose, paired clamp (n = 10):\nthe output arm falls in 10/10")
    ax = fig.add_subplot(gs[1, 0])
    H = pd.read_csv(f"{OUT}/new/GSE148640_HC_response_stats.tsv", sep="\t", index_col=0)
    HC = pd.read_pickle(f"{OUT}/GSE148640_diff_genes.pkl"); rng = np.random.default_rng(3)
    for i, (v, c) in enumerate(((HC.loc["PER1"].values, IS), (D.loc["PER1"].values, AD))):
        ax.scatter(np.full(v.shape, i) + rng.uniform(-.13, .13, len(v)), v, s=22, color=c, alpha=.85)
        ax.hlines(np.median(v), i - .26, i + .26, color="k", lw=1.8)
        ax.text(i, max(v) + .22, f"{int((v > 0).sum())}/{len(v)} up", ha="center", fontsize=6.5)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["hydrocortisone\nvs placebo", "clamp\nvs fasting"], fontsize=7)
    ax.axhline(0, color="k", lw=.8); ax.set_ylabel("PER1, per-person change (log$_2$)")
    src(ax, "PER1 \u00b7 adipose \u00b7 GSE148640 and GSE26637, n = 10 each", y=1.06)
    panel(ax, "D", "Two opposite perturbations\nof the same gene")
    ax = fig.add_subplot(gs[1, 1])
    gl3 = [g for g in ["PER1", "PER3", "BMAL1", "CIART", "NR1D1", "DBP", "TEF", "ZBTB16", "FKBP5", "KLF15"] if g in HC.index]
    rng2 = np.random.default_rng(7)
    for i, g in enumerate(gl3):
        v = HC.loc[g].values
        ax.scatter(np.full(v.shape, i) + rng2.uniform(-.15, .15, len(v)), v, s=11,
                   color=AD, alpha=.75)   # un solo color: el signo lo da la posicion respecto al cero
        ax.hlines(np.median(v), i - .28, i + .28, color="k", lw=1.5)
        ax.text(i, 3.45, f"{int((v > 0).sum())}/10", ha="center", fontsize=6)
    ax.axhline(0, color="k", lw=.8); ax.set_xticks(range(len(gl3)))
    ax.set_xticklabels(gl3, rotation=45, fontsize=7)
    ax.set_ylabel("hydrocortisone \u2212 placebo (log$_2$)"); ax.set_ylim(-1.8, 3.9)
    src(ax, "GSE148640 \u00b7 adipose \u00b7 both biopsies at 07:00", y=1.06)
    panel(ax, "E", "Same time of day: the action is\ngene-specific, not a clock rotation")
    ax = fig.add_subplot(gs[1, 2])
    L = pd.read_csv(f"{OUT}/new/GSE17992_dex_log2ratio_muscle.tsv", sep="\t", index_col=0) \
        if os.path.exists(f"{OUT}/new/GSE17992_dex_log2ratio_muscle.tsv") \
        else pd.read_csv(f"{OUT}/new/GSE17992_dex_log2ratio_muscle.tsv.gz", sep="\t", index_col=0)
    gl4 = [g for g in ["KLF15", "DBP", "CRY2", "CLOCK", "TEF", "HLF", "PER3", "BMAL1", "ZBTB16"] if g in L.index]
    for i, g in enumerate(gl4):
        v = L.loc[g].dropna().values
        ax.scatter(np.full(v.shape, i) + np.random.default_rng(i).uniform(-.15, .15, len(v)), v, s=11, color=MU, alpha=.75)
        ax.hlines(np.median(v), i - .28, i + .28, color="k", lw=1.5)
    ax.axhline(0, color="k", lw=.8); ax.set_xticks(range(len(gl4))); ax.set_xticklabels(gl4, rotation=45, fontsize=7)
    ax.set_ylabel("dexamethasone, 4 days (log$_2$)")
    src(ax, "GSE17992 \u00b7 muscle \u00b7 10 men \u00b7 4 d dexamethasone", y=1.06)
    panel(ax, "F", "Dexamethasone in muscle (n = 10)")
    ax = fig.add_subplot(gs[2, 0])
    Xc_ = pd.read_pickle(f"{OUT}/new/GSE66446_log2cpm.pkl"); Xc_.index = [{"ARNTL": "BMAL1"}.get(i, i) for i in Xc_.index]
    mdc = pd.read_pickle(f"{OUT}/new/GSE66446_meta.pkl").set_index("gsm").reindex(Xc_.columns)
    cu = mdc.c0.str.contains("Cushing").values; Xc2 = center(Xc_); rng = np.random.default_rng(4)
    gl5 = [g for g in ["PER3", "DBP", "PER2", "TEF", "HLF", "NR1D1", "CIART", "BMAL1"] if g in Xc2.index]
    for i, g in enumerate(gl5):
        for k, (m, mk, c) in enumerate(((~cu, "s", GR), (cu, "o", IS))):
            v = Xc2.loc[g, m].values - Xc2.loc[g, ~cu].median(); xx = i + k * .3 - .15
            ax.scatter(np.full(v.shape, xx) + rng.uniform(-.06, .06, len(v)), v, s=11, color=c, marker=mk, alpha=.85)
            ax.hlines(np.median(v), xx - .12, xx + .12, color="k", lw=1.2)
    ax.axhline(0, color="k", lw=.8); ax.set_xticks(range(len(gl5))); ax.set_xticklabels(gl5, rotation=45, fontsize=7)
    ax.set_ylabel("log$_2$ (centred on control)")
    ax.legend([plt.Line2D([], [], marker="s", ls="", color=GR), plt.Line2D([], [], marker="o", ls="", color=IS)],
              ["control (n = 11)", "Cushing (n = 5)"], frameon=False, loc="upper left")
    src(ax, "GSE66446 \u00b7 adipose \u00b7 5 Cushing vs 11 control", y=1.06)
    panel(ax, "G", "Cushing syndrome, adipose")
    ax = fig.add_subplot(gs[2, 1])
    S = pd.read_csv(f"{OUT}/new/signatures_log2_common_genes.tsv", sep="\t", index_col=0) \
        if os.path.exists(f"{OUT}/new/signatures_log2_common_genes.tsv") \
        else pd.read_csv(f"{OUT}/new/signatures_log2_common_genes.tsv.gz", sep="\t", index_col=0)
    cols = {"insulin-resistant adipose": "IR adipose (GSE13070, 35 vs 6)", "Cushing": "Cushing adipose (GSE66446)",
            "hydrocortisone": "hydrocortisone adipose (GSE148640)", "TZD, 3 months": "TZD 3mo adipose (GSE13070)"}
    gl6 = [g for g in SETS["Clock output"]["genes"] if g in S.index]   # registro, no redefinir
    T = S.loc[gl6, list(cols.values())]; T.columns = list(cols.keys()); x = np.arange(len(gl6)); w = .2
    for k, (c, col_) in enumerate(zip(T.columns, [MU, IS, GC, AD])):
        ax.bar(x + (k - 1.5) * w, T[c], w, color=col_, label=c)
    ax.axhline(0, color="k", lw=.8); ax.set_xticks(x); ax.set_xticklabels(gl6, rotation=45, fontsize=7)
    ax.set_ylabel("log$_2$"); ax.legend(frameon=False, loc="upper left", ncol=1)
    r1 = stats.spearmanr(T["insulin-resistant adipose"], T["Cushing"]).statistic
    src(ax, "adipose \u00b7 four contrasts, same nine clock genes", y=1.06)
    panel(ax, "H")
    # I  la TZD frente a la firma, panel movido desde la Fig. 4
    ax = fig.add_subplot(gs[2, 2])
    gl7 = [g for g in SETS["Clock output"]["genes"] if g in S.index]   # el mismo conjunto que el panel H
    ir = S.loc[gl7, "IR adipose (GSE13070, 35 vs 6)"]; tz = S.loc[gl7, "TZD 3mo adipose (GSE13070)"]
    ax.scatter(ir, tz, s=34, color=AD); trend(ax, ir, tz)
    for g in gl7: ax.annotate(g, (ir[g], tz[g]), fontsize=6, xytext=(3, 3), textcoords="offset points")
    ax.axhline(0, color="k", lw=.8); ax.axvline(0, color="k", lw=.8)
    rz = stats.spearmanr(ir, tz); sl = np.polyfit(ir, tz, 1)[0]
    ax.set(xlabel="insulin-resistant − sensitive (log$_2$)", ylabel="TZD, 3 months (log$_2$)")
    src(ax, f"GSE13070 \u00b7 adipose \u00b7 {len(gl7)} clock-output genes \u00b7 slope = {sl:+.3f}", y=1.06)
    panel(ax, "I")
    tidy(fig); fig.savefig(f"{OUT}/figs/F3.png", dpi=400, bbox_inches="tight")
    fig.savefig(f"{OUT}/figs/F3.pdf", bbox_inches="tight"); plt.close(fig)

# ------------------------------------------------------------------ F4
def fig4():
    fig = plt.figure(figsize=(17, 12.5)); gs = fig.add_gridspec(3, 3, hspace=.62, wspace=.72)
    T = pd.read_csv(f"{OUT}/new/AAGMEx_TF_activity_vs_clock.tsv", sep="\t")
    ax = fig.add_subplot(gs[0, 0])
    P = pd.concat([T.nsmallest(8, "rho"), T.nlargest(8, "rho")])
    # lollipop: el segmento ancla cada valor al cero, que es lo que una nube de puntos pierde
    sz = 18 + 110 * (P.n_dianas - P.n_dianas.min()) / (P.n_dianas.max() - P.n_dianas.min() + 1e-9)
    for i, (rho_, n_) in enumerate(zip(P.rho, sz)):
        ax.plot([0, rho_], [i, i], color=GC if rho_ > 0 else MU, lw=1.1, zorder=1)
    ax.scatter(P.rho, range(len(P)), s=sz, c=[GC if v > 0 else MU for v in P.rho],
               edgecolors="k", linewidths=.4, zorder=2)
    ax.set_yticks(range(len(P))); ax.set_yticklabels(P.TF, fontsize=7); ax.axvline(0, color="k", lw=.8)
    ax.set_xlabel("ρ  TF activity ~ clock output"); ax.set_xlim(-.55, .62)
    hs = [plt.Line2D([], [], marker="o", ls="", color="#bbbbbb", markeredgecolor="k",
                     markersize=np.sqrt(18 + 110 * f) / 1.6) for f in (0, 1)]
    ax.legend(hs, [f"{int(P.n_dianas.min())} targets", f"{int(P.n_dianas.max())} targets"],
              frameon=False, loc="upper left", handletextpad=1.0, labelspacing=1.3,
              borderpad=.2, bbox_to_anchor=(.55, .42))
    src(ax, "AAGMEx \u00b7 adipose \u00b7 n = 232 \u00b7 134 TFs with a DoRothEA A\u2013C regulon")
    panel(ax, "A", "TF activity (DoRothEA A–C regulons):\nFOXO3 leads 134 factors, NF-κB trails")
    H = pd.read_csv(f"{OUT}/new/AAGMEx_NR3C1_enrichment_extra.tsv", sep="\t"); Hh = H[H.lib == "h_all"].copy()
    Hh["lab"] = Hh.set.str.replace("HALLMARK_", "").str.replace("_", " ").str.capitalize()
    Hh = pd.concat([Hh.nlargest(3, "z"), Hh.nsmallest(11, "z")]).sort_values("z")
    ax = fig.add_subplot(gs[0, 1])
    sz = 20 + 130 * (Hh.n - Hh.n.min()) / (Hh.n.max() - Hh.n.min() + 1e-9)
    s_ = ax.scatter(Hh.z, range(len(Hh)), s=sz, c=Hh.z, cmap="RdBu_r", vmin=-9, vmax=9, edgecolors="k", linewidths=.4)
    ax.set_yticks(range(len(Hh))); ax.set_yticklabels(Hh.lab, fontsize=6.5); ax.axvline(0, color="k", lw=.8)
    cb = plt.colorbar(s_, ax=ax, shrink=.75, pad=.03); cb.set_label("z", fontsize=7)
    ax.set_xlabel("z vs size-matched random sets")
    src(ax, "AAGMEx \u00b7 adipose \u00b7 n = 232 \u00b7 genome-wide ranking by partial \u03c1 with NR3C1")
    panel(ax, "B", "Hallmark: the receptor axis is anti-inflammatory\nand runs against OXPHOS and adipogenesis")
    # C  HPO x tissue: aqui si hay matriz, de modo que el mapa de calor es la representacion
    # correcta; con un solo ranking seria una columna unica y se usaria un lollipop.
    HP = pd.read_csv(f"{OUT}/new/HPO_by_tissue.tsv", sep="\t")
    P2 = HP.pivot_table(index="term", columns="tissue", values="z").reindex(columns=["adipose", "muscle", "islets"])
    P2 = P2.loc[P2.adipose.sort_values().index]
    P2.index = [t.replace("HP_", "").replace("_", " ").lower().capitalize() for t in P2.index]
    ax = fig.add_subplot(gs[0, 2])
    im = ax.imshow(P2.values, cmap="RdBu_r", vmin=-7, vmax=7, aspect="auto")
    ax.set_xticks(range(3)); ax.set_xticklabels(P2.columns, fontsize=7)
    ax.set_yticks(range(len(P2))); ax.set_yticklabels(P2.index, fontsize=6.5)
    for i in range(P2.shape[0]):
        for j in range(P2.shape[1]):
            v = P2.values[i, j]
            ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=6.5,
                    color="w" if abs(v) > 4 else "k")
    cb = plt.colorbar(im, ax=ax, shrink=.8, pad=.03); cb.set_label("z", fontsize=7)
    src(ax, "receptor axis of each tissue \u00b7 AAGMEx adipose and muscle, GSE159984 islets")
    panel(ax, "C")
    ax = fig.add_subplot(gs[1, 0])
    G = pd.read_csv(f"{OUT}/new/AAGMEx_GCreceptor_vs_adiposity.tsv", sep="\t")
    x = np.arange(len(G)); w = .38
    ax.bar(x - w / 2, G.rho_vs_LEP, w, color=GR, label="ρ with LEP (adiposity)")
    ax.bar(x + w / 2, G.rho_vs_clock, w, color=AD, label="ρ with clock output")
    ax.axhline(0, color="k", lw=.8); ax.set_xticks(x); ax.set_xticklabels(G.gene, rotation=45, fontsize=7)
    ax.set_ylabel("ρ"); ax.legend(frameon=False, loc="lower left")
    src(ax, "AAGMEx \u00b7 adipose \u00b7 n = 232 \u00b7 adjusted for sex and age")
    panel(ax, "D")
    ax = fig.add_subplot(gs[1, 1])
    R = pd.read_csv(f"{OUT}/new/AAGMEx_rank_vs_NR3C1.tsv", sep="\t", index_col=0)["rho"] \
        if os.path.exists(f"{OUT}/new/AAGMEx_rank_vs_NR3C1.tsv") \
        else pd.read_csv(f"{OUT}/new/AAGMEx_rank_vs_NR3C1.tsv.gz", sep="\t", index_col=0)["rho"]
    ax.hist(R.dropna(), bins=60, color="#dddddd", edgecolor="none")
    for g, c in (("FOXO3", AD), ("NFKB1", MU), ("HSD11B1", IS)):
        if g in R.index:
            ax.axvline(R[g], color=c, lw=1.5); ax.text(R[g], ax.get_ylim()[1] * .9, g, rotation=90, fontsize=6.5, color=c, ha="right", va="top")
    ax.set_xlim(-.6, .6); ax.set(xlabel="ρ with NR3C1 (adjusted)", ylabel="genes")
    src(ax, "AAGMEx \u00b7 adipose \u00b7 34,020 genes")
    panel(ax, "E")
    E2 = pd.read_csv(f"{OUT}/new/AAGMEx_NR3C1_enrichment.tsv", sep="\t")
    for k, (lib, let, ttl) in enumerate(((("c2_cp_reactome"), "F", "Reactome"), (("c5_go_bp"), "G", "GO biological process"))):
        q = E2[E2.lib == lib].copy(); q = pd.concat([q.nlargest(4, "z"), q.nsmallest(10, "z")]).sort_values("z")
        q["lab"] = q.set.str.split("_").str[1:].str.join(" ").str.capitalize().str.slice(0, 38)
        a = fig.add_subplot(gs[2, k])
        sz = 20 + 130 * (q.n - q.n.min()) / (q.n.max() - q.n.min() + 1e-9)
        s3 = a.scatter(q.z, range(len(q)), s=sz, c=q.z, cmap="RdBu_r", vmin=-9, vmax=9, edgecolors="k", linewidths=.4)
        a.set_yticks(range(len(q))); a.set_yticklabels(q.lab, fontsize=6); a.axvline(0, color="k", lw=.8)
        cb = plt.colorbar(s3, ax=a, shrink=.75, pad=.03); cb.set_label("z", fontsize=7)
        a.set_xlabel("z vs size-matched random sets")
        src(a, "AAGMEx \u00b7 adipose \u00b7 n = 232 \u00b7 " + ttl + " \u00b7 dot size: genes in the set")
        panel(a, let)
    tidy(fig); fig.savefig(f"{OUT}/figs/F4.png", dpi=400, bbox_inches="tight")
    fig.savefig(f"{OUT}/figs/F4.pdf", bbox_inches="tight"); plt.close(fig)

def fig2():
    import subprocess; subprocess.run([sys.executable, os.path.join(HERE, "fig3_eje.py")], check=True)
    for ext in ("png", "pdf"): os.replace(f"{OUT}/figs/F3_eje.{ext}", f"{OUT}/figs/F2.{ext}")

def fig5():
    import subprocess; subprocess.run([sys.executable, os.path.join(HERE, "fig7_modelo.py")], check=True)
    for ext in ("png", "pdf"): os.replace(f"{OUT}/figs/F7_modelo.{ext}", f"{OUT}/figs/F5.{ext}")

if __name__ == "__main__":
    fig1(); print("F1 ok"); fig2(); print("F2 ok"); fig3(); print("F3 ok")
    fig4(); print("F4 ok"); fig5(); print("F5 ok")
