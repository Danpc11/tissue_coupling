import os as _os, sys as _sys
_RAW = _os.environ.get("T2D_RAW", "data/raw"); _OUT = _os.environ.get("T2D_OUT", "results")   # las mismas rutas que run_all.sh
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../lib"))
import pandas as pd, numpy as np, warnings; warnings.filterwarnings("ignore", category=FutureWarning); warnings.filterwarnings("ignore", category=DeprecationWarning)
from scipy import stats
D=pd.read_csv(f"{_OUT}/new/amplitude_phase_distributions_two_cohorts.tsv",sep="\t",index_col=0)
D["coh"]=D.cohort.str.replace(" músculo basal","",regex=False).str.replace(" músculo","",regex=False).str.strip()
rng=np.random.default_rng(0); ROWS=[]
def permutation_ratio(a,b,fn,n=5000):
    """Razon de dispersiones fn(a)/fn(b) y su P por PERMUTACION de las etiquetas de grupo (no es un
    bootstrap). Si la dispersion de algun grupo es cero la razon no esta definida: devuelve nan en
    vez de un infinito, y las permutaciones degeneradas no cuentan."""
    fa,fb=fn(a),fn(b)
    if not (fa>0 and fb>0): return np.nan,np.nan
    obs=fa/fb; pool=np.concatenate([a,b]); na=len(a); null=[]
    for _ in range(n):
        p=rng.permutation(pool); x,y=fn(p[:na]),fn(p[na:])
        if x>0 and y>0: null.append(x/y)
    null=np.array(null)
    if len(null)<n/2: return obs,np.nan
    return obs,(np.sum(np.abs(np.log(null))>=abs(np.log(obs)))+1)/(len(null)+1)
def ad_p(a,b):
    """P de Anderson-Darling. scipy lo acota entre 0.001 y 0.25: en los extremos es una cota."""
    p=stats.anderson_ksamp([a,b]).pvalue
    return p, (">= 0.25" if p>=0.25 else "<= 0.001" if p<=0.001 else f"{p:.4f}")
IQR=lambda v: np.subtract(*np.percentile(v,[75,25]))
RNG=lambda v: np.percentile(v,95)-np.percentile(v,5)
for est in ("amp_800","amp_clock","amp_free"):
    print(f"\n=== {est} ===")
    for coh,q in D.groupby("coh"):
        gs={k:v[est].dropna().values for k,v in q.groupby("group")}
        if len(gs)!=2: continue
        (ka,a),(kb,b)=sorted(gs.items())
        print(f"  {coh}  ({ka} n={len(a)} vs {kb} n={len(b)})")
        print(f"     mediana   {np.median(a):.2f} vs {np.median(b):.2f}   MW P {stats.mannwhitneyu(a,b).pvalue:.3f}")
        print(f"     IQR       {IQR(a):.2f} vs {IQR(b):.2f}   razón {IQR(a)/IQR(b):.2f}   P permutación {permutation_ratio(a,b,IQR)[1]:.4f}")
        print(f"     rango 5-95 {RNG(a):.2f} vs {RNG(b):.2f}   razón {RNG(a)/RNG(b):.2f}   P {permutation_ratio(a,b,RNG)[1]:.4f}")
        print(f"     sd        {np.std(a,ddof=1):.2f} vs {np.std(b,ddof=1):.2f}   Levene P {stats.levene(a,b,center='median').pvalue:.4f}   Fligner P {stats.fligner(a,b).pvalue:.4f}")
        print(f"     cuantiles a: {np.percentile(a,[10,25,50,75,90]).round(2)}")
        print(f"     cuantiles b: {np.percentile(b,[10,25,50,75,90]).round(2)}")
        print(f"     Kolmogorov-Smirnov P {stats.ks_2samp(a,b).pvalue:.4f}   Anderson-Darling P {ad_p(a,b)[1]}")
        ROWS.append(dict(estimator=est,cohort=coh,group_a=ka,n_a=len(a),group_b=kb,n_b=len(b),median_a=np.median(a),median_b=np.median(b),
                         p_mannwhitney=stats.mannwhitneyu(a,b).pvalue,iqr_ratio=IQR(a)/IQR(b) if IQR(b)>0 else np.nan,p_iqr=permutation_ratio(a,b,IQR)[1],
                         range_5_95_ratio=RNG(a)/RNG(b) if RNG(b)>0 else np.nan,p_range=permutation_ratio(a,b,RNG)[1],
                         p_levene=stats.levene(a,b,center="median").pvalue,p_fligner=stats.fligner(a,b).pvalue,
                         p_ks=stats.ks_2samp(a,b).pvalue,p_anderson_darling=ad_p(a,b)[0],anderson_darling_note=ad_p(a,b)[1]))
pd.DataFrame(ROWS).to_csv(f"{_OUT}/new/amplitude_dispersion.tsv",sep="\t",index=False)
