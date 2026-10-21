#!/usr/bin/env python3
"""Extended Data figures ED1-ED7 (English, Nature style, uppercase panel letters).

Descriptions live in FIGURE_CAPTIONS.md; panels carry only their letter and a provenance line.
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
                     "ytick.labelsize": 7, "legend.fontsize": 6.5})

def panel(ax, L): ax.text(-0.18, 1.13, L, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")
def src(ax, t): ax.text(0, 1.03, t, transform=ax.transAxes, fontsize=6.3, color="#666666", style="italic", va="bottom")
def tidy(f):
    for a in f.axes: a.spines[["top", "right"]].set_visible(False)
def save(fig, name):
    fig.savefig(f"{OUT}/figs/{name}.pdf", bbox_inches="tight")
    fig.savefig(f"{OUT}/figs/{name}.png", dpi=400, bbox_inches="tight"); plt.close(fig)

def sc(X, g):
    gg = [r for r in (resolve(x, X.index) for x in g) if r]
    Z = X.loc[gg]; return ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T.mean(0)

# ---------------------------------------------------------------- ED2
def ed2():
    """Muscle amplitude does not replicate."""
    D = pd.read_csv(f"{OUT}/new/amplitude_phase_distributions_two_cohorts.tsv", sep="\t", index_col=0)
    D["coh"] = D.cohort.str.replace(" músculo basal", "", regex=False).str.replace(" músculo", "", regex=False).str.replace(" M", "", regex=False).str.strip()
    fig = plt.figure(figsize=(14, 8.4)); gs = fig.add_gridspec(2, 3, hspace=.6, wspace=.62)
    ests = [("amp_800", "800-gene projection"), ("amp_clock", "clock genes, anchored"), ("amp_free", "clock genes, free ordering")]
    rng = np.random.default_rng(0)
    for k, (col, lab) in enumerate(ests):
        ax = fig.add_subplot(gs[0, k]); x0 = 0
        for coh, q in D.groupby("coh"):
            for gi, (grp, qq) in enumerate(q.groupby("group")):
                v = qq[col].dropna().values
                ax.scatter(np.full(v.shape, x0 + gi * .8) + rng.uniform(-.15, .15, len(v)), v,
                           s=9, color=MU if gi else GR, alpha=.65)
                ax.hlines(np.median(v), x0 + gi * .8 - .26, x0 + gi * .8 + .26, color="k", lw=1.5)
            gs_ = [qq[col].dropna().values for _, qq in q.groupby("group")]
            if len(gs_) == 2:
                ax.text(x0 + .4, D[col].max() * 1.02, f"P = {stats.mannwhitneyu(*gs_).pvalue:.3f}", ha="center", fontsize=6.5)
            x0 += 2.4
        ax.set_xticks([.4, 2.8]); ax.set_xticklabels(sorted(D.coh.unique()), fontsize=6.5)
        ax.set_ylabel("inferred amplitude"); src(ax, lab); panel(ax, "ABC"[k])
    ax = fig.add_subplot(gs[1, 0])
    for k, (col, lab) in enumerate(ests):
        rows = []
        for coh, q in D.groupby("coh"):
            gs_ = [qq[col].dropna().values for _, qq in q.groupby("group")]
            if len(gs_) == 2:
                rows.append((coh, np.median(gs_[1]) - np.median(gs_[0]), stats.mannwhitneyu(*gs_).pvalue))
        for i, (coh, d_, p_) in enumerate(rows):
            ax.scatter(d_, k * 2.2 + i, s=34, color=MU if p_ < .05 else GR)
            ax.text(.98, k * 2.2 + i, f"P = {p_:.3f}", fontsize=6, va="center", ha="right", transform=ax.get_yaxis_transform())
    ax.axvline(0, color="k", lw=.8)
    ax.set_yticks([k * 2.2 + i for k in range(3) for i in range(2)])
    ax.set_yticklabels([f"{lab.replace('clock genes, ','clock ').replace('800-gene projection','800-gene')} · {coh}" for _, lab in ests for coh in sorted(D.coh.unique())], fontsize=6)
    ax.tick_params(axis='y', pad=1)
    ax.set_xlabel("resistant − sensitive (inferred amplitude)"); ax.set_xlim(-0.95, 1.9)
    src(ax, "only the unvalidated 800-gene projection separates the groups"); panel(ax, "D")
    ax = fig.add_subplot(gs[1, 1])
    for coh, q in D.groupby("coh"):
        ax.scatter(q.amp_free, q.amp_800, s=10, alpha=.6, label=coh)
    r = stats.spearmanr(D.amp_free, D.amp_800)
    ax.set(xlabel="amplitude, free ordering", ylabel="amplitude, 800-gene projection")
    ax.legend(frameon=False, loc="upper right", title=f"ρ = {r.statistic:.2f}", title_fontsize=7)
    src(ax, "the two estimators measure different things"); panel(ax, "E")
    # F  dispersion: la prueba correcta si el mecanismo es perdida de entrainment
    ax = fig.add_subplot(gs[1, 2])
    rng2 = np.random.default_rng(1)
    RNG = lambda v: np.percentile(v, 95) - np.percentile(v, 5)
    rows = []
    for k, (col, lab) in enumerate(ests):
        for coh, q in D.groupby("coh"):
            gsx = [qq[col].dropna().values for _, qq in q.groupby("group")]
            if len(gsx) != 2: continue
            a, b = gsx[0], gsx[1]
            obs = RNG(a) / RNG(b); pool = np.concatenate([a, b]); na = len(a)
            null = np.array([RNG(p[:na]) / RNG(p[na:]) for p in (rng2.permutation(pool) for _ in range(3000))])
            pv = (np.sum(np.abs(np.log(null)) >= abs(np.log(obs))) + 1) / 3001
            rows.append((f"{lab} · {coh}", obs, pv, stats.levene(a, b, center="median").pvalue))
    for i, (lab, o, pv, lv) in enumerate(rows):
        ax.scatter(o, i, s=36, color=MU if pv < .05 else GR)
        ax.text(.98, i, f"P = {pv:.3f}  (Levene {lv:.3f})", fontsize=5.8, va="center", ha="right", transform=ax.get_yaxis_transform())
    ax.axvline(1, color="k", lw=.8)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0].replace("clock genes, ", "clock ").replace("800-gene projection", "800-gene") for r in rows], fontsize=6)
    ax.tick_params(axis="y", pad=1)
    ax.set_xlabel("5–95 range ratio, resistant / sensitive"); ax.set_xlim(.85, 3.2)
    src(ax, "scale test: the median is blind to a dispersion effect")
    panel(ax, "F")
    tidy(fig); save(fig, "ED2")

# ---------------------------------------------------------------- ED3
def ed3():
    """Methodological sensitivity."""
    fig = plt.figure(figsize=(14.5, 4.6)); gs = fig.add_gridspec(1, 3, wspace=.62)
    A = pd.read_csv(f"{OUT}/new/AUDIT_global_axis.tsv", sep="\t")
    ax = fig.add_subplot(gs[0, 0])
    ax.scatter(A.null_med, A.rho, s=34, color=GR, label="before sample centring")
    ax.scatter(A.null_med_c, A.rho_centrado, s=34, color=AD, label="after sample centring")
    for _, r in A.iterrows():
        ax.plot([r.null_med, r.null_med_c], [r.rho, r.rho_centrado], color="#cccccc", lw=.8, zorder=0)
    ax.axvline(0, color="k", lw=.8, ls=":")
    ax.set(xlabel="median ρ of a random gene with the clock score", ylabel="ρ  NR3C1 ~ clock output")
    ax.legend(frameon=False, loc="lower right")
    src(ax, "ten cohort-tissue matrices"); panel(ax, "A")
    ax = fig.add_subplot(gs[0, 1])
    S = pd.read_csv(f"{OUT}/new/msigdb_symbol_coverage.tsv", sep="\t")
    x = np.arange(len(S)); w = .38
    ax.bar(x - w / 2, S.en_msigdb_crudo, w, color=GR, label="platform symbol")
    ax.bar(x + w / 2, S.en_msigdb_mapeado, w, color=AD, label="Entrez → current symbol")
    for i, g in enumerate(S.ganancia):
        if abs(g) > 100: ax.text(i, max(S.en_msigdb_mapeado[i], S.en_msigdb_crudo[i]) * 1.02, f"+{int(g)}", ha="center", fontsize=6.5)
    ax.set_xticks(x); ax.set_xticklabels([c.replace(" adiposo", "\nadipose").replace(" músculo", "\nmuscle") for c in S.cohorte], rotation=0, fontsize=5.6)
    ax.set_ylabel("genes matched in MSigDB"); ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0, -0.16), ncol=2)
    src(ax, "MSigDB v2026.2, H + C2 + C5"); panel(ax, "B")
    ax = fig.add_subplot(gs[0, 2])
    T = pd.read_csv(f"{OUT}/new/symbol_audit_renamed_genes.tsv", sep="\t")
    cnt = T.groupby(T.columns[0]).size().sort_values()
    ax.barh(range(len(cnt)), cnt.values, color=IS)
    ax.set_yticks(range(len(cnt))); ax.set_yticklabels([i.split(" ")[0] + "\n" + " ".join(i.split(" ")[1:]) for i in cnt.index], fontsize=6); ax.tick_params(axis="y", pad=1)
    ax.set_xlabel("probes carrying a retired symbol")
    src(ax, "per platform"); panel(ax, "C")
    tidy(fig); save(fig, "ED3")

# ---------------------------------------------------------------- ED5
def ed5():
    """Within-person clamp and TZD, full diagnostic."""
    S = pd.read_csv(f"{OUT}/clamp_scores_change.tsv", sep="\t")
    fig = plt.figure(figsize=(13, 7)); gs = fig.add_gridspec(2, 2, hspace=.6, wspace=.32)
    con = list(dict.fromkeys(S.contrast))
    for k, c in enumerate(con[:4]):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        q = S[S.contrast == c].sort_values("median_delta")
        ax.barh(range(len(q)), q.median_delta, color=[MU if p < .05 else GR for p in q.p_wilcoxon])
        for i, (d_, p_, f_) in enumerate(zip(q.median_delta, q.p_wilcoxon, q.frac_down)):
            ax.text(d_ + (.02 if d_ > 0 else -.02), i, f"P = {p_:.1g}   {f_*100:.0f}% down",
                    fontsize=5.8, va="center", ha="left" if d_ > 0 else "right")
        ax.set_yticks(range(len(q))); ax.set_yticklabels(q.score, fontsize=6.5)
        ax.axvline(0, color="k", lw=.8); ax.set_xlabel("median per-person change (score units)")
        ax.set_xlim(min(q.median_delta) * 1.9 - .1, max(q.median_delta) * 1.9 + .1)
        src(ax, f"{c} · n = {int(q.n_pairs.iloc[0])} pairs"); panel(ax, "ABCD"[k])
    tidy(fig); save(fig, "ED5")

# ---------------------------------------------------------------- ED6
def ed6():
    """The empirical GRAS signature and its failure to validate."""
    T = pd.read_csv(f"{OUT}/new/GSE148640_HC_response_stats.tsv", sep="\t", index_col=0)
    G = pd.read_csv(f"{OUT}/new/GRAS_signature.tsv", sep="\t")
    fig = plt.figure(figsize=(13, 4.6)); gs = fig.add_gridspec(1, 3, wspace=.42)
    ax = fig.add_subplot(gs[0, 0])
    pcol = "p" if "p" in T.columns else "P"
    ax.scatter(T["median"], -np.log10(T[pcol].clip(lower=1e-6)), s=5, color="#dddddd")
    sel = T.loc[[g for g in G.gene if g in T.index]]
    ax.scatter(sel["median"], -np.log10(sel[pcol].clip(lower=1e-6)), s=12, color=IS)
    ax.set(xlabel="hydrocortisone − placebo (log$_2$)", ylabel="−log$_{10}$ P (Wilcoxon, n = 10)")
    src(ax, "GSE148640 · adipose · 12,641 genes; 68 selected (10/10 consistent)"); panel(ax, "A")
    ax = fig.add_subplot(gs[0, 1])
    nup = T.n_up.value_counts().sort_index()
    ax.bar(nup.index, nup.values, color=[IS if i in (0, 10) else GR for i in nup.index])
    ax.set(xlabel="participants changing in the same direction (of 10)", ylabel="genes")
    src(ax, "no gene reaches FDR < 0.10 with n = 10"); panel(ax, "B")
    ax = fig.add_subplot(gs[0, 2])
    tests = ["Cushing\nadipose", "dexamethasone\nmuscle", "dexamethasone\nadipose", "GSE13070\nIR adipose"]
    gras = [0.39, -0.06, -0.33, 0.12]; nine = [0.14, 0.48, 0.48, 0.21]
    x = np.arange(4); w = .38
    ax.bar(x - w / 2, gras, w, color=GR, label="empirical GRAS (68 genes)")
    ax.bar(x + w / 2, nine, w, color=AD, label="pre-specified set (9 genes)")
    ax.axhline(0, color="k", lw=.8); ax.set_xticks(x); ax.set_xticklabels(tests, fontsize=6.5)
    ax.set_ylabel("score difference (log$_2$)"); ax.legend(frameon=False, loc="lower left")
    src(ax, "external validation: none of these cohorts was used to build GRAS"); panel(ax, "C")
    tidy(fig); save(fig, "ED6")

# ---------------------------------------------------------------- ED7
def ed7():
    """Model robustness."""
    import forced_clock as fc
    fig = plt.figure(figsize=(13, 4.6)); gs = fig.add_gridspec(1, 3, wspace=.42)
    ax = fig.add_subplot(gs[0, 0])
    for dt, c in ((0.2, GR), (0.1, GC), (0.05, AD), (0.02, IS)):
        amps = []
        for fk in np.linspace(0, .2, 9):
            n = 60; rng = np.random.default_rng(3); taus = rng.normal(24, 1.2, n); om = 2 * np.pi / taus
            ph = rng.uniform(0, 2 * np.pi, n); rr = np.ones(n)
            for t in np.arange(0, 15 * 24, dt):
                d = 2 * np.pi * t / 24 - ph
                ph = ph + dt * (om + fk * np.sin(d) / np.maximum(rr, 1e-3))
                rr = np.clip(rr + dt * (.05 * rr * (1 - rr) + fk * np.cos(d)), 1e-3, 10)
            amps.append(fc.ensemble(ph % (2 * np.pi), rr)[0])
        ax.plot(np.linspace(0, .2, 9), amps, "o-", ms=3, color=c, label=f"dt = {dt} h")
    ax.set(xlabel="effective forcing  F · K", ylabel="ensemble amplitude")
    ax.legend(frameon=False, loc="lower right")
    src(ax, "integrator convergence"); panel(ax, "A")
    ax = fig.add_subplot(gs[0, 1])
    lv = {"morning": [], "uniform": []}; av = {"morning": [], "uniform": []}
    out = ["PER3", "DBP", "HLF", "TEF", "NR1D1", "NR1D2", "PER2", "CIART"]
    for hsd, key in ((1.5, "morning"), (7.5, "uniform")):
        for s_ in range(5):
            res = []
            for Km, seed in ((1.4, 100 + s_), (0.5, 200 + s_)):
                rng = np.random.default_rng(seed); K = np.clip(rng.normal(Km, .15, 70), .05, None)
                hour = rng.normal(8, hsd, 70) if hsd < 7 else rng.uniform(0, 24, 70)
                rows = []
                for i, k in enumerate(K):
                    phi, r = fc.simulate_population(.04, k, n_cells=40, days=12, seed=int(rng.integers(1e6)))
                    a, p, R = fc.ensemble(phi, r); rows.append(dict(subject=i, K=k, hour=hour[i], amp_ens=a, phase_ens=p, R=R))
                dd = pd.DataFrame(rows); res.append((dd, fc.synthetic_expression(dd, noise=.6, seed=int(rng.integers(1e6)))))
            lv[key].append(res[1][1].loc[out].mean().mean() - res[0][1].loc[out].mean().mean())
            av[key].append(res[1][0].amp_ens.mean() - res[0][0].amp_ens.mean())
    x = np.arange(2); w = .38
    ax.bar(x - w / 2, [np.mean(lv[k]) for k in lv], w, yerr=[np.std(lv[k]) for k in lv], capsize=4, color=IS, label="LEVEL difference")
    ax.bar(x + w / 2, [np.mean(av[k]) for k in av], w, yerr=[np.std(av[k]) for k in av], capsize=4, color=GR, label="AMPLITUDE difference")
    ax.axhline(0, color="k", lw=.8); ax.set_xticks(x); ax.set_xticklabels(["morning window", "uniform hours"], fontsize=7)
    ax.set_ylabel("low − high coupling"); ax.legend(frameon=False, loc="upper left")
    src(ax, "five seeds per condition; mean ± s.d."); panel(ax, "B")
    ax = fig.add_subplot(gs[0, 2])
    for sig, c in ((0.6, GC), (1.2, AD), (2.0, MU)):
        Rs = []
        for fk in np.linspace(0, .12, 10):
            n = 80; rng = np.random.default_rng(11); taus = rng.normal(24, sig, n); om = 2 * np.pi / taus
            ph = rng.uniform(0, 2 * np.pi, n); rr = np.ones(n); dt = .05
            for t in np.arange(0, 15 * 24, dt):
                d = 2 * np.pi * t / 24 - ph
                ph = ph + dt * (om + fk * np.sin(d) / np.maximum(rr, 1e-3))
                rr = np.clip(rr + dt * (.05 * rr * (1 - rr) + fk * np.cos(d)), 1e-3, 10)
            Rs.append(fc.ensemble(ph % (2 * np.pi), rr)[2])
        ax.plot(np.linspace(0, .12, 10), Rs, "o-", ms=3, color=c, label=f"period s.d. = {sig} h")
    ax.set(xlabel="effective forcing  F · K", ylabel="ensemble synchrony R")
    ax.legend(frameon=False, loc="lower right")
    src(ax, "sensitivity to the dispersion of intrinsic periods"); panel(ax, "C")
    tidy(fig); save(fig, "ED7")

# ---------------------------------------------------------------- ED1
def ed1():
    """GTEx: the receptor-clock axis across tissues."""
    PH = pd.read_pickle(f"{OUT}/new/GTEx_pheno.pkl")
    agemid = {'20-29': 25, '30-39': 35, '40-49': 45, '50-59': 55, '60-69': 65, '70-79': 75}
    SA = pd.read_csv("data/gtex/SampleAttributes.txt", sep="\t", low_memory=False)
    SA["donor"] = SA.SAMPID.str.split("-").str[:2].str.join("-")
    TIS = [("adipose sc", "GTEx_adipose_sc", "Adipose - Subcutaneous", AD),
           ("adipose visceral", "GTEx_adipose_visc", "Adipose - Visceral (Omentum)", IS),
           ("skeletal muscle", "GTEx_muscle", "Muscle - Skeletal", MU)]
    fig = plt.figure(figsize=(13, 4.6)); gs = fig.add_gridspec(1, 3, wspace=.42)
    rows = []
    for lab, fn, smtsd, c in TIS:
        X = pd.read_pickle(f"{OUT}/new/{fn}.pkl"); Xc = X.sub(X.mean(0), axis=1)
        sa = SA[SA.SMTSD == smtsd].drop_duplicates("donor").set_index("donor").reindex(X.columns)
        p_ = PH.reindex(X.columns)
        H = pd.get_dummies(p_.DTHHRDY.fillna(-1).astype(int), prefix="h").values.astype(float)
        isch = sa.SMTSISCH.astype(float); isch = isch.fillna(isch.median()).values
        rin = sa.SMRIN.astype(float); rin = rin.fillna(rin.median()).values
        base = [(p_.SEX == 1).astype(float).values, p_.AGE.map(agemid).values.astype(float)] + \
               [H[:, k] for k in range(H.shape[1])]
        co = sc(Xc, SETS["Clock output"]["genes"]); nr = Xc.loc["NR3C1"]
        def res(y, cv):
            D = np.column_stack([np.ones(len(y))] + cv)
            b, *_ = np.linalg.lstsq(D, np.asarray(y, float), rcond=None); return np.asarray(y, float) - D @ b
        rows.append(dict(tissue=lab, color=c, n=X.shape[1],
                         base=stats.spearmanr(res(nr, base), res(co, base)).statistic,
                         full=stats.spearmanr(res(nr, base + [isch, rin]), res(co, base + [isch, rin])).statistic,
                         isch_clock=stats.spearmanr(isch, co).statistic,
                         rin_clock=stats.spearmanr(rin, co).statistic))
    R = pd.DataFrame(rows)
    ax = fig.add_subplot(gs[0, 0]); x = np.arange(len(R)); w = .38
    ax.bar(x - w / 2, R.base, w, color=GR, label="sex, age, Hardy scale")
    ax.bar(x + w / 2, R.full, w, color=AD, label="+ ischaemic time and RIN")
    ax.set_xticks(x); ax.set_xticklabels([f"{t}\n(n = {int(n)})" for t, n in zip(R.tissue, R.n)], fontsize=6.5)
    ax.set_ylabel("ρ  NR3C1 ~ clock output")
    ax.set_ylim(0, max(R.base.max(), R.full.max()) * 1.55)
    ax.legend(frameon=False, loc="upper left", borderaxespad=.2)
    src(ax, "GTEx v11 · sample-centred expression"); panel(ax, "A")
    ax = fig.add_subplot(gs[0, 1])
    ax.bar(x - w / 2, R.isch_clock, w, color=MU, label="ischaemic time")
    ax.bar(x + w / 2, R.rin_clock, w, color=GC, label="RIN")
    ax.axhline(0, color="k", lw=.8); ax.set_xticks(x); ax.set_xticklabels(R.tissue, fontsize=6.5)
    ax.set_ylabel("ρ with clock-output score")
    lo, hi = min(R.isch_clock.min(), R.rin_clock.min()), max(R.isch_clock.max(), R.rin_clock.max())
    ax.set_ylim(lo - abs(lo) * .25 - .03, hi + (hi - lo) * .55)
    ax.legend(frameon=False, loc="upper left", borderaxespad=.2)
    src(ax, "GTEx v11 · same three tissues"); panel(ax, "B")
    ax = fig.add_subplot(gs[0, 2])
    X = pd.read_pickle(f"{OUT}/new/GTEx_adipose_visc.pkl"); Xc = X.sub(X.mean(0), axis=1)
    co = sc(Xc, SETS["Clock output"]["genes"])
    ax.scatter(Xc.loc["NR3C1"], co, s=5, alpha=.4, color=IS)
    b = np.polyfit(np.asarray(Xc.loc["NR3C1"], float), np.asarray(co, float), 1)
    xs = np.linspace(Xc.loc["NR3C1"].min(), Xc.loc["NR3C1"].max(), 40)
    ax.plot(xs, np.polyval(b, xs), "k--", lw=1.1)
    ax.set(xlabel="NR3C1 (sample-centred log$_2$)", ylabel="clock-output score (z)")
    src(ax, "GTEx visceral adipose · n = 587"); panel(ax, "C")
    tidy(fig); save(fig, "ED1")

if __name__ == "__main__":
    for f, n in ((ed1, "ED1"), (ed2, "ED2"), (ed3, "ED3"), (ed5, "ED5"), (ed6, "ED6"), (ed7, "ED7")):
        try: f(); print(n, "ok")
        except Exception as e: print(n, "FAILED:", type(e).__name__, e)
