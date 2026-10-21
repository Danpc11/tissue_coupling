import pandas as pd, numpy as np, warnings; warnings.filterwarnings("ignore")
from scipy import stats
D=pd.read_csv("results/new/amplitude_phase_distributions_two_cohorts.tsv",sep="\t",index_col=0)
D["coh"]=D.cohort.str.replace(" músculo basal","",regex=False).str.replace(" músculo","",regex=False).str.strip()
rng=np.random.default_rng(0)
def boot_ratio(a,b,fn,n=5000):
    obs=fn(a)/fn(b)
    null=[]
    pool=np.concatenate([a,b]); na=len(a)
    for _ in range(n):
        p=rng.permutation(pool); null.append(fn(p[:na])/fn(p[na:]))
    null=np.array(null); return obs,(np.sum(np.abs(np.log(null))>=abs(np.log(obs)))+1)/(n+1)
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
        print(f"     IQR       {IQR(a):.2f} vs {IQR(b):.2f}   razón {IQR(a)/IQR(b):.2f}   P permutación {boot_ratio(a,b,IQR)[1]:.4f}")
        print(f"     rango 5-95 {RNG(a):.2f} vs {RNG(b):.2f}   razón {RNG(a)/RNG(b):.2f}   P {boot_ratio(a,b,RNG)[1]:.4f}")
        print(f"     sd        {np.std(a,ddof=1):.2f} vs {np.std(b,ddof=1):.2f}   Levene P {stats.levene(a,b,center='median').pvalue:.4f}   Fligner P {stats.fligner(a,b).pvalue:.4f}")
        print(f"     cuantiles a: {np.percentile(a,[10,25,50,75,90]).round(2)}")
        print(f"     cuantiles b: {np.percentile(b,[10,25,50,75,90]).round(2)}")
        print(f"     Kolmogorov-Smirnov P {stats.ks_2samp(a,b).pvalue:.4f}   Anderson-Darling P {stats.anderson_ksamp([a,b]).pvalue:.4f}")
