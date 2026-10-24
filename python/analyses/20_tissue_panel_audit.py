import pandas as pd, numpy as np, sys, warnings; warnings.filterwarnings("ignore")
from scipy import stats
sys.path.insert(0,"python/analyses"); sys.path.insert(0,"python/lib")
from importlib import import_module
m=import_module("19_gtex_tissue_panel") if False else None
exec(open("python/analyses/19_gtex_tissue_panel.py").read().split('if __name__')[0])
PH=pd.read_pickle("results/new/GTEx_pheno.pkl")
SA=pd.read_csv("data/gtex/SampleAttributes.txt",sep="\t",low_memory=False)
SA["donor"]=SA.SAMPID.str.split("-").str[:2].str.join("-")
print(f"{'tejido':26s} {'rho':>7s} {'pct':>6s} {'nulo sd':>8s} {'|r|>rho':>8s} {'CD68':>7s} {'PTPRC':>7s} {'rho|comp':>9s}")
rows=[]
for f,smtsd in TIS:
    L=load(f)
    sa=SA[SA.SMTSD==smtsd].drop_duplicates("donor").set_index("donor").reindex(L.columns)
    p=PH.reindex(L.columns)
    H=pd.get_dummies(p.DTHHRDY.fillna(-1).astype(int),prefix="h").values.astype(float)
    isch=sa.SMTSISCH.astype(float); isch=isch.fillna(isch.median()).values
    rin=sa.SMRIN.astype(float); rin=rin.fillna(rin.median()).values
    base=[(p.SEX==1).astype(float).values,p.AGE.map(AGE).values.astype(float)]+[H[:,k] for k in range(H.shape[1])]+[isch,rin]
    def resid(y,cv):
        D=np.column_stack([np.ones(len(y))]+cv); Q,_=np.linalg.qr(D); y=np.asarray(y,float); return y-Q@(Q.T@y)
    co,_=sc(L,SETS["Clock output"]["genes"]); cor=resid(co,base)
    rho=stats.spearmanr(resid(L.loc["NR3C1"],base),cor).statistic
    rng=np.random.default_rng(0)
    null=np.array([stats.spearmanr(resid(L.loc[g],base),cor).statistic for g in rng.choice(L.index,1500,replace=False)])
    # ajuste adicional por composicion: marcador inmune y de musculo liso donde exista
    comp=[L.loc[g].values for g in ("CD68","PTPRC") if g in L.index]
    rho_c=stats.spearmanr(resid(L.loc["NR3C1"],base+comp),resid(co,base+comp)).statistic if comp else np.nan
    rcd=stats.spearmanr(L.loc["CD68"],co).statistic if "CD68" in L.index else np.nan
    rpt=stats.spearmanr(L.loc["PTPRC"],co).statistic if "PTPRC" in L.index else np.nan
    print(f"{smtsd[:26]:26s} {rho:+7.3f} {np.mean(null<rho)*100:5.0f}% {null.std():8.3f} {np.mean(np.abs(null)>abs(rho))*100:7.1f}% {rcd:+7.2f} {rpt:+7.2f} {rho_c:+9.3f}")
    rows.append(dict(tissue=smtsd,rho=rho,pct=np.mean(null<rho)*100,null_sd=null.std(),
                     frac_abs_gt=np.mean(np.abs(null)>abs(rho))*100,rho_comp_adj=rho_c,
                     clock_vs_CD68=rcd,clock_vs_PTPRC=rpt))
pd.DataFrame(rows).to_csv("results/new/GTEx_tissue_panel_audit.tsv",sep="\t",index=False)
