import pandas as pd, numpy as np, sys, warnings; warnings.filterwarnings("ignore")
from scipy import stats
sys.path.insert(0,"python/lib")
exec(open("python/analyses/19_gtex_tissue_panel.py").read().split('if __name__')[0])
PH=pd.read_pickle("results/new/GTEx_pheno.pkl")
SA=pd.read_csv("data/gtex/SampleAttributes.txt",sep="\t",low_memory=False)
SA["donor"]=SA.SAMPID.str.split("-").str[:2].str.join("-")
EXTRA=[("GTEx_adipose_visc","Adipose - Visceral (Omentum)"),("GTEx_muscle","Muscle - Skeletal")]
print(f"{'tejido':26s} {'PC1 %':>6s} {'sd nulo':>8s} {'rho':>7s} {'pct':>5s} | {'tras quitar 3 PC':>16s} {'sd nulo':>8s} {'pct':>5s}")
rows=[]
items=[(f,s,None) for f,s in TIS]+[(None,s,f) for f,s in EXTRA]
for f,smtsd,pk in items:
    L=load(f) if f else pd.read_pickle(f"results/new/{pk}.pkl").pipe(lambda d: d.sub(d.mean(0),axis=1))
    sa=SA[SA.SMTSD==smtsd].drop_duplicates("donor").set_index("donor").reindex(L.columns)
    p=PH.reindex(L.columns)
    H=pd.get_dummies(p.DTHHRDY.fillna(-1).astype(int),prefix="h").values.astype(float)
    isch=sa.SMTSISCH.astype(float); isch=isch.fillna(isch.median()).values
    rin=sa.SMRIN.astype(float); rin=rin.fillna(rin.median()).values
    base=[(p.SEX==1).astype(float).values,p.AGE.map(AGE).values.astype(float)]+[H[:,k] for k in range(H.shape[1])]+[isch,rin]
    Z=L.sub(L.mean(1),axis=0).div(L.std(1)+1e-9,axis=0)
    U,S,Vt=np.linalg.svd(Z.values-Z.values.mean(1,keepdims=True),full_matrices=False)
    var=S**2/ (S**2).sum()
    pcs=[Vt[k] for k in range(3)]
    def run(cv,tag):
        D=np.column_stack([np.ones(L.shape[1])]+cv); Q,_=np.linalg.qr(D)
        r=lambda y: np.asarray(y,float)-Q@(Q.T@np.asarray(y,float))
        co,_=sc(L,SETS["Clock output"]["genes"]); cor=r(co)
        rho=stats.spearmanr(r(L.loc["NR3C1"]),cor).statistic
        rng=np.random.default_rng(0)
        null=np.array([stats.spearmanr(r(L.loc[g]),cor).statistic for g in rng.choice(L.index,800,replace=False)])
        return rho,null.std(),np.mean(null<rho)*100
    a=run(base,"base"); b=run(base+pcs,"pc")
    print(f"{smtsd[:26]:26s} {var[0]*100:5.1f}% {a[1]:8.3f} {a[0]:+7.3f} {a[2]:4.0f}% | {b[0]:+16.3f} {b[1]:8.3f} {b[2]:4.0f}%")
    rows.append(dict(tissue=smtsd,pc1_var=var[0]*100,sd_null=a[1],rho=a[0],pct=a[2],
                     rho_pc=b[0],sd_null_pc=b[1],pct_pc=b[2]))
R=pd.DataFrame(rows); R.to_csv("results/new/GTEx_panel_pc_adjusted.tsv",sep="\t",index=False)
print(f"\ncorrelación entre varianza del PC1 y anchura del nulo: r = {stats.pearsonr(R.pc1_var,R.sd_null)[0]:+.2f}")
