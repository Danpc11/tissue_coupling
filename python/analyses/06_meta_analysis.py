"""Metaanalisis de efectos aleatorios del eje receptor-reloj.

DerSimonian-Laird sobre Fisher-z, igual que el metaanalisis de la Tabla S17 del manuscrito.
Una estimacion por cohorte y tejido, residualizada frente a las covariables de cada cohorte.
Salidas: results/new/meta_per_cohort.tsv y meta_pooled.tsv
"""
import pandas as pd, numpy as np, sys, gzip, importlib.util, os
from scipy import stats
os.environ["T2D_OUT"]="results"; sys.path.insert(0,"python/lib"); from gene_sets import SETS, resolve
def sc(X,g):
    gg=[r for r in (resolve(x,X.index) for x in g) if r]; Z=X.loc[gg]; return ((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T.mean(0)
def pres(y,cv):
    D=np.column_stack([np.ones(len(y))]+cv) if cv else np.ones((len(y),1))
    b,*_=np.linalg.lstsq(D,np.asarray(y,float),rcond=None); return np.asarray(y,float)-D@b
CL=SETS["Clock output"]["genes"]; GC=SETS["Glucocorticoid"]["genes"]
PH=pd.read_pickle("results/new/GTEx_pheno.pkl"); agemid={'20-29':25,'30-39':35,'40-49':45,'50-59':55,'60-69':65,'70-79':75}
def gtexcov(X):
    p=PH.reindex(X.columns); H=pd.get_dummies(p.DTHHRDY.fillna(-1).astype(int),prefix="h").values.astype(float)
    cv=[(p.SEX==1).astype(float).values,p.AGE.map(agemid).values.astype(float)]+[H[:,k] for k in range(H.shape[1])]
    return cv
A=pd.read_pickle("results/new/AAGMEx_adipose_log2.pkl"); M=pd.read_pickle("results/new/AAGMEx_muscle_log2.pkl")
with gzip.open("data/raw/GSE95674_series_matrix.txt.gz","rt") as f:
    for l in f:
        if l.startswith("!Sample_title"): t=[x.strip('"') for x in l.rstrip().split("\t")[1:]]
        if l.startswith("!Sample_characteristics") and "gender" in l: sx=[x.strip('"').split(": ")[1] for x in l.rstrip().split("\t")[1:]]
        if l.startswith("!Sample_characteristics") and "age" in l: ag=[float(x.strip('"').split(": ")[1]) for x in l.rstrip().split("\t")[1:]]
cov=pd.DataFrame({"sex":sx,"age":ag},index=t).reindex(A.columns); aag=[(cov.sex=="M").astype(float).values,cov.age.values]
Ea=pd.read_pickle("results/new/GSE41168_adipose.pkl"); Em=pd.read_pickle("results/new/GSE41168_muscle.pkl")
sp=importlib.util.spec_from_file_location("c5","python/analyses/05_clamp_within_person.py"); c5=importlib.util.module_from_spec(sp); sp.loader.exec_module(c5)
E,d=c5.load(); mus13=d[(d.tissue=="muscle")&(d.tzd=="pre")&(d.clamp=="pre")].drop_duplicates("subject").set_index("subject")
adi13=d[(d.tissue=="adipose")&(d.tzd=="pre")].drop_duplicates("subject").set_index("subject")
G1=pd.read_pickle("results/new/GTEx_adipose_sc.pkl"); G2=pd.read_pickle("results/new/GTEx_adipose_visc.pkl"); G3=pd.read_pickle("results/new/GTEx_muscle.pkl")
rows=[]
def add(coh,tis,X,cv,x_lab="NR3C1"):
    x=pres(X.loc[x_lab],cv); y=pres(sc(X,CL),cv); g=pres(sc(X,GC),cv)
    for lab,v in (("NR3C1~reloj",x),("GC~reloj",g)):
        r=stats.spearmanr(v,y).statistic; n=len(y)
        z=np.arctanh(np.clip(r,-.999,.999)); se=1/np.sqrt(n-3-len(cv))
        rows.append(dict(contraste=lab,cohorte=coh,tejido=tis,n=n,rho=r,z=z,se=se))
add("AAGMEx","adiposo",A,aag); add("AAGMEx","músculo",M,aag)
add("GSE41168","adiposo",Ea,[]); add("GSE41168","músculo",Em,[])
add("GSE13070","adiposo",E[list(adi13.gsm)],[]); add("GSE13070","músculo",E[list(mus13.gsm)],[])
add("GTEx","adiposo sc",G1,gtexcov(G1)); add("GTEx","adiposo visceral",G2,gtexcov(G2)); add("GTEx","músculo",G3,gtexcov(G3))
R=pd.DataFrame(rows)
def dl(sub):
    w=1/sub.se**2; m=(w*sub.z).sum()/w.sum(); Q=(w*(sub.z-m)**2).sum(); k=len(sub)
    tau2=max(0,(Q-(k-1))/(w.sum()-(w**2).sum()/w.sum())); w2=1/(sub.se**2+tau2); m2=(w2*sub.z).sum()/w2.sum()
    se2=np.sqrt(1/w2.sum()); I2=max(0,(Q-(k-1))/Q*100) if Q>0 else 0
    return np.tanh(m2),np.tanh(m2-1.96*se2),np.tanh(m2+1.96*se2),Q,k,I2,2*(1-stats.norm.cdf(abs(m2/se2)))
out=[]
for lab in R.contraste.unique():
    for tej in ("todos","adiposo","músculo"):
        sub=R[(R.contraste==lab)&((R.tejido.str.contains("adiposo") if tej=="adiposo" else R.tejido.str.contains("músculo")) if tej!="todos" else True)]
        if len(sub)<2: continue
        rr,lo,hi,Q,k,I2,p=dl(sub)
        print(f"{lab:12s} {tej:8s} k={k}: ρ combinado {rr:+.3f} (IC95 {lo:+.3f} a {hi:+.3f}), P {p:.2e}, Q={Q:.1f}, I²={I2:.0f}%")
        out.append(dict(contraste=lab,grupo=tej,k=k,rho=rr,lo=lo,hi=hi,p=p,Q=Q,I2=I2))
R.to_csv("results/new/meta_per_cohort.tsv",sep="\t",index=False); pd.DataFrame(out).to_csv("results/new/meta_pooled.tsv",sep="\t",index=False)
print("\nPor cohorte:"); print(R.round(3).to_string(index=False))
