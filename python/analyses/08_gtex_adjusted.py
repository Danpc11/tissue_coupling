"""GTEx ajustado por tiempo isquemico y RIN, sobre matrices centradas por muestra.

El tiempo isquemico correlaciona con el score de reloj (rho 0,21-0,29) y con NR3C1, asi que hay
que descartarlo como explicacion. Requiere GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt
"""
import pandas as pd, numpy as np, sys, warnings; warnings.filterwarnings("ignore")
from scipy import stats
sys.path.insert(0,"python/lib"); from gene_sets import SETS, resolve
SA=pd.read_csv("data/gtex/SampleAttributes.txt",sep="\t",low_memory=False)
SA["donor"]=SA.SAMPID.str.split("-").str[:2].str.join("-")
PH=pd.read_pickle("results/new/GTEx_pheno.pkl"); agemid={'20-29':25,'30-39':35,'40-49':45,'50-59':55,'60-69':65,'70-79':75}
def sc(X,g):
    gg=[r for r in (resolve(x,X.index) for x in g) if r]; Z=X.loc[gg]; return ((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T.mean(0)
def pres(y,cv):
    D=np.column_stack([np.ones(len(y))]+cv); b,*_=np.linalg.lstsq(D,np.asarray(y,float),rcond=None); return np.asarray(y,float)-D@b
TIS={"adiposo sc":("GTEx_adipose_sc","Adipose - Subcutaneous"),"visceral":("GTEx_adipose_visc","Adipose - Visceral (Omentum)"),"músculo":("GTEx_muscle","Muscle - Skeletal")}
for lab,(fn,smtsd) in TIS.items():
    X=pd.read_pickle(f"results/new/{fn}.pkl")
    sa=SA[SA.SMTSD==smtsd].drop_duplicates("donor").set_index("donor").reindex(X.columns)
    p=PH.reindex(X.columns)
    H=pd.get_dummies(p.DTHHRDY.fillna(-1).astype(int),prefix="h").values.astype(float)
    base=[(p.SEX==1).astype(float).values,p.AGE.map(agemid).values.astype(float)]+[H[:,k] for k in range(H.shape[1])]
    isch=sa.SMTSISCH.astype(float).fillna(sa.SMTSISCH.astype(float).median()).values
    rin=sa.SMRIN.astype(float).fillna(sa.SMRIN.astype(float).median()).values
    Xc=X.sub(X.mean(0),axis=1)   # centrado por muestra (corrección §30)
    co=sc(Xc,SETS["Clock output"]["genes"]); gc=sc(Xc,SETS["Glucocorticoid"]["genes"]); nr=Xc.loc["NR3C1"]
    print(f"\n{lab} (n={X.shape[1]}, RIN {np.nanmedian(rin):.1f}, isquemia {np.nanmedian(isch):.0f} min):")
    print(f"   isquemia vs reloj ρ {stats.spearmanr(isch,co).statistic:+.3f} | vs NR3C1 ρ {stats.spearmanr(isch,nr).statistic:+.3f}")
    print(f"   RIN      vs reloj ρ {stats.spearmanr(rin,co).statistic:+.3f} | vs NR3C1 ρ {stats.spearmanr(rin,nr).statistic:+.3f}")
    for name,cv in (("sexo, edad, Hardy",base),("+ isquemia y RIN",base+[isch,rin])):
        r=stats.spearmanr(pres(nr,cv),pres(co,cv)); r2=stats.spearmanr(pres(gc,cv),pres(co,cv))
        print(f"   {name:22s}: NR3C1~reloj {r.statistic:+.3f} (P {r.pvalue:.1e}) | GC~reloj {r2.statistic:+.3f}")
