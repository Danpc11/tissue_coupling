import pandas as pd, numpy as np, sys, warnings; warnings.filterwarnings("ignore")
from scipy import stats
sys.path.insert(0,"python/lib"); import phase as ph
X=pd.read_pickle("results/new/GSE13070_entrez.pkl"); Xc=X.sub(X.mean(0),axis=1)
S=pd.read_csv("results/clamp_samples.tsv",sep="\t")
S=S[(S.contrast=="clamp_baseline")&(S.condition=="pre")].drop_duplicates("gsm")
S=S[S.gsm.isin(Xc.columns)]
genes=[g for g in ph.CLK if g in Xc.index]
print(f"GSE13070 músculo basal: {len(S)} personas ({(S.sensitivity=='resistant').sum()} R / {(S.sensitivity=='sensitive').sum()} S), {len(genes)} genes del reloj")
rng=np.random.default_rng(0); B=80; n=max(4,int(len(genes)*0.6))
cols=list(S.gsm); A=[]
for _ in range(B):
    g=list(rng.choice(genes,n,replace=False))
    Z=Xc.loc[g,cols]; Z=((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T
    ang=ph.unsup(Z.values)
    a=np.array([np.mean(Z.loc[x].values*np.cos(ang)) for x in g])
    b=np.array([np.mean(Z.loc[x].values*np.sin(ang)) for x in g])
    v=Z.values
    th=np.arctan2(v.T@b,v.T@a)
    A.append(np.array([np.dot(v[:,i],a*np.cos(th[i])+b*np.sin(th[i])) for i in range(v.shape[1])])/np.sqrt(n))
A=np.array(A)                                   # B x personas
noise=A.std(0); mean=A.mean(0)
grp=S.set_index("gsm").sensitivity.reindex(cols).values
nr,ns=noise[grp=="resistant"],noise[grp=="sensitive"]
print(f"\nRUIDO del estimador por persona (sd entre {B} reinferencias con subconjuntos de genes):")
print(f"  resistentes mediana {np.median(nr):.3f}  |  sensibles mediana {np.median(ns):.3f}  |  razón {np.median(nr)/np.median(ns):.2f}")
print(f"  Mann-Whitney P = {stats.mannwhitneyu(nr,ns).pvalue:.3f}   Levene P = {stats.levene(nr,ns,center='median').pvalue:.3f}")
mr,ms=mean[grp=="resistant"],mean[grp=="sensitive"]
print(f"\nAMPLITUD media por persona:")
print(f"  var resistentes {mr.var(ddof=1):.3f} | sensibles {ms.var(ddof=1):.3f} | razón {mr.var(ddof=1)/ms.var(ddof=1):.2f}")
vr=mr.var(ddof=1)-np.mean(nr**2); vs=ms.var(ddof=1)-np.mean(ns**2)
print(f"  descontando el ruido medido de cada grupo: var real R {vr:.3f}, S {vs:.3f}, razón {vr/vs:.2f}" if vs>0 else "  varianza de sensibles menor que su ruido")
RNG=lambda v: np.percentile(v,95)-np.percentile(v,5)
obs=RNG(mr)/RNG(ms); pool=np.concatenate([mr,ms]); na=len(mr)
null=np.array([RNG(p[:na])/RNG(p[na:]) for p in (rng.permutation(pool) for _ in range(3000))])
print(f"  razón de rango 5-95 {obs:.2f}, P permutación {(np.sum(np.abs(np.log(null))>=abs(np.log(obs)))+1)/3001:.4f}")
pd.DataFrame(dict(gsm=cols,group=grp,amp=mean,noise=noise)).to_csv("results/new/GSE13070_amplitude_noise_by_group.tsv",sep="\t",index=False)
