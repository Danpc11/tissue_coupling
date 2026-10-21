#!/usr/bin/env python3
"""El eje receptor-reloj en el panel ampliado de tejidos de GTEx v11.

Pregunta: el acoplamiento entre el nivel del receptor y la salida del reloj, ¿es propio de los
tejidos metabolicos o una propiedad general de los tejidos perifericos humanos?

Mismo procedimiento en todos: CPM log2, genes con >5 cuentas en la mitad de las muestras, centrado
por muestra, y residuales frente a sexo, edad, clase de Hardy, tiempo isquemico y RIN. El nulo de
cada tejido se estima con 400 genes al azar para situar a NR3C1 en su distribucion.
"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(HERE,"../lib"))
from gene_sets import SETS, resolve
OUT=os.environ.get("T2D_OUT","results"); GTEX=os.environ.get("T2D_GTEX","data/gtex")
AGE={'20-29':25,'30-39':35,'40-49':45,'50-59':55,'60-69':65,'70-79':75}
TIS=[("liver","Liver"),("whole_blood","Whole Blood"),("stomach","Stomach"),
     ("small_intestine_terminal_ileum","Small Intestine - Terminal Ileum"),
     ("colon_transverse","Colon - Transverse"),("colon_sigmoid","Colon - Sigmoid"),
     ("nerve_tibial","Nerve - Tibial"),("artery_tibial","Artery - Tibial"),
     ("adrenal_gland","Adrenal Gland")]

def load(f):
    X=pd.read_csv(f"{GTEX}/{f}.gct.gz",sep="\t",skiprows=2)
    X=X.drop(columns=["Name"]).rename(columns={"Description":"gene"}).groupby("gene").sum()
    X.columns=["-".join(c.split("-")[:2]) for c in X.columns]
    X=X.loc[:,~X.columns.duplicated()]
    X=X[(X>5).sum(1)>=X.shape[1]*.5]
    L=np.log2(X/X.sum()*1e6+1)
    return L.sub(L.mean(0),axis=1)

def sc(X,g):
    gg=[r for r in (resolve(x,X.index) for x in g) if r]
    Z=X.loc[gg]; return ((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T.mean(0), len(gg)

if __name__=="__main__":
    PH=pd.read_pickle(f"{OUT}/new/GTEx_pheno.pkl")
    SA=pd.read_csv(f"{GTEX}/SampleAttributes.txt",sep="\t",low_memory=False)
    SA["donor"]=SA.SAMPID.str.split("-").str[:2].str.join("-")
    rows=[]
    for f,smtsd in TIS:
        L=load(f)
        sa=SA[SA.SMTSD==smtsd].drop_duplicates("donor").set_index("donor").reindex(L.columns)
        p=PH.reindex(L.columns)
        H=pd.get_dummies(p.DTHHRDY.fillna(-1).astype(int),prefix="h").values.astype(float)
        isch=sa.SMTSISCH.astype(float); isch=isch.fillna(isch.median()).values
        rin=sa.SMRIN.astype(float); rin=rin.fillna(rin.median()).values
        cv=[(p.SEX==1).astype(float).values,p.AGE.map(AGE).values.astype(float)]+ \
           [H[:,k] for k in range(H.shape[1])]+[isch,rin]
        D=np.column_stack([np.ones(L.shape[1])]+cv); Q,_=np.linalg.qr(D)
        r=lambda y: np.asarray(y,float)-Q@(Q.T@np.asarray(y,float))
        co,n1=sc(L,SETS["Clock output"]["genes"]); gc,n2=sc(L,SETS["Glucocorticoid"]["genes"])
        cor=r(co)
        rho=stats.spearmanr(r(L.loc["NR3C1"]),cor)
        rgc=stats.spearmanr(r(gc),cor)
        rng=np.random.default_rng(0)
        null=np.array([stats.spearmanr(r(L.loc[g]),cor).statistic
                       for g in rng.choice(L.index,400,replace=False)])
        rows.append(dict(tissue=smtsd,n=L.shape[1],n_clock=n1,
                         rho_NR3C1=rho.statistic,P=rho.pvalue,
                         rho_GC=rgc.statistic,P_GC=rgc.pvalue,
                         null_median=np.median(null),pct=np.mean(null<rho.statistic)*100))
        print(f"{smtsd:34s} n={L.shape[1]:4d}  NR3C1~clock {rho.statistic:+.3f} (P {rho.pvalue:.1e})"
              f"  GC~clock {rgc.statistic:+.3f}  nulo {np.median(null):+.3f}  pct {np.mean(null<rho.statistic)*100:.0f}%")
    pd.DataFrame(rows).to_csv(f"{OUT}/new/GTEx_tissue_panel.tsv",sep="\t",index=False)
