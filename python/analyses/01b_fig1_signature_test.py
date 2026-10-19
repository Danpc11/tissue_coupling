#!/usr/bin/env python3
"""Figura 1: la expresión diferencial convencional no recupera la firma.

Pregunta de la figura: ¿puede el análisis convencional —comparar grupos dentro de un tejido— dar
la firma glucocorticoide y de salida del reloj sobre la que descansa el resto del artículo?
Respuesta: no. Ni gen a gen (ningún miembro de los conjuntos pasa FDR 0,1 en ninguna cohorte) ni,
en la mayoría de cohortes, a nivel de conjunto. La firma solo aparece cuando se estima por muestra
y se correlaciona dentro del individuo, que es lo que hacen las figuras 2 a 6.

Cohortes de la figura, las mismas cinco del manuscrito:
  músculo  GSE182120 (T2D vs NGT), GSE25462 (historia familiar de T2D), GSE159984 (transversal)
  adiposo  GSE78721, GSE26637
GSE13070 NO entra en esta figura: su contraste adiposo es 35 vs 6 y domina el panel con ~3.000
genes significativos que son consecuencia del desbalance, no del efecto. GSE13070 se usa en las
figuras 4 y 6, donde el diseño es pareado dentro de la persona y el desbalance no aplica.
Este script corre las que estén en data/raw y avisa de las que falten.

Método: moderación de varianza estilo limma (ajuste por momentos de la distribución a priori sobre
log s², contracción de cada varianza hacia ella), Benjamini-Hochberg, y enriquecimiento de conjunto
por z frente a 2.000 conjuntos del mismo tamaño muestreados del universo medido.

Salidas: results/fig1_de_per_cohort.tsv, fig1_signature_genes.tsv, fig1_set_enrichment.tsv
"""
import os, sys, gzip, glob, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(HERE,"../lib"))
from gene_sets import SETS, resolve
RAW=os.environ.get("T2D_RAW","data/raw"); OUT=os.environ.get("T2D_OUT","results")

def moderated_t(X, grp):
    """t moderada estilo limma. grp: booleano, True = caso."""
    a=X.values[:,grp]; b=X.values[:,~grp]; na,nb=a.shape[1],b.shape[1]
    df=na+nb-2
    s2=((a.var(1,ddof=1)*(na-1))+(b.var(1,ddof=1)*(nb-1)))/df
    s2=np.clip(s2,1e-12,None)
    # a priori por momentos sobre log(s2)  (Smyth 2004, aproximación)
    z=np.log(s2); e=z-stats.psi(df/2)+np.log(df/2) if hasattr(stats,"psi") else z-(np.log(df/2)-1/df)
    from scipy.special import digamma, polygamma
    e=np.log(s2)-digamma(df/2)+np.log(df/2)
    ebar=e.mean(); var_e=e.var(ddof=1)-polygamma(1,df/2)
    if var_e<=0: d0,s02=np.inf,np.exp(ebar)
    else:
        # resolver polygamma(1,d0/2)=var_e
        lo,hi=1e-4,1e4
        for _ in range(80):
            mid=np.sqrt(lo*hi)
            if polygamma(1,mid/2)>var_e: lo=mid
            else: hi=mid
        d0=np.sqrt(lo*hi); s02=np.exp(ebar+digamma(d0/2)-np.log(d0/2))
    s2_post=(d0*s02+df*s2)/(d0+df) if np.isfinite(d0) else np.full_like(s2,s02)
    t=(a.mean(1)-b.mean(1))/np.sqrt(s2_post*(1/na+1/nb))
    dft=df+(d0 if np.isfinite(d0) else 1e6)
    p=2*stats.t.sf(np.abs(t),dft)
    return pd.DataFrame({"logFC":a.mean(1)-b.mean(1),"t":t,"P":p,
                         "FDR":stats.false_discovery_control(p)},index=X.index)

def set_z(res, genes, n=2000, seed=0):
    rng=np.random.default_rng(seed)
    g=[r for r in (resolve(x,res.index) for x in genes) if r]
    if len(g)<3: return np.nan,0
    obs=res.t[g].mean(); null=np.array([res.t.sample(len(g)).mean() for _ in range(n)])
    return (obs-null.mean())/null.std(), len(g)

if __name__=="__main__":
    import importlib.util
    spec=importlib.util.spec_from_file_location("c5",os.path.join(HERE,"05_clamp_within_person.py"))
    c5=importlib.util.module_from_spec(spec); spec.loader.exec_module(c5)
    COH={}
    # --- GSE25462 (GPL570): músculo, diabetes (10) frente a historia familiar negativa (15).
    #     Trae SI por persona (FSIGT), que se usa como fenotipo continuo en el panel de distribuciones.
    if os.path.exists(f"{OUT}/new/GSE25462_log2.pkl"):
        X=pd.read_pickle(f"{OUT}/new/GSE25462_log2.pkl"); md=pd.read_pickle(f"{OUT}/new/GSE25462_meta.pkl").set_index("gsm")
        keep=md.index[md.fh.isin(["DM","Family history negative"])]
        COH["GSE25462 · músculo"]=(X[list(keep)],(md.fh.reindex(keep)=="DM").values,"músculo")
    # --- GSE26637 (GPL570): adiposo, ayuno, 5 resistentes / 5 sensibles (el brazo de clamp no entra)
    if os.path.exists(f"{OUT}/new/GSE26637_log2.pkl"):
        X=pd.read_pickle(f"{OUT}/new/GSE26637_log2.pkl"); md=pd.read_pickle(f"{OUT}/new/GSE26637_meta.pkl").set_index("gsm")
        keep=md.index[md.stim=="fasting"]
        COH["GSE26637 · adiposo"]=(X[list(keep)],(md.status.reindex(keep)=="resistant").values,"adiposo")
    # --- GSE78721 (PrimeView): adiposo de extremidad, 68 T2D / 62 normales
    if os.path.exists(f"{OUT}/new/GSE78721_log2.pkl"):
        X=pd.read_pickle(f"{OUT}/new/GSE78721_log2.pkl"); md=pd.read_pickle(f"{OUT}/new/GSE78721_meta.pkl").set_index("gsm")
        COH["GSE78721 · adiposo"]=(X,(md.dx.reindex(X.columns)!="Normal").values,"adiposo")
    # --- GSE182120 (HTA 2.0): músculo 25 T2D / 24 NGT, centrado por laboratorio
    if os.path.exists(f"{OUT}/new/GSE182120_log2.pkl"):
        X=pd.read_pickle(f"{OUT}/new/GSE182120_log2.pkl"); md=pd.read_pickle(f"{OUT}/new/GSE182120_meta.pkl").set_index("gsm")
        Xc=X.copy()
        for lab,q in md.groupby("lab"): Xc[list(q.index)]=X[list(q.index)].sub(X[list(q.index)].mean(1),axis=0)
        COH["GSE182120 · músculo"]=(Xc,(md.disease.reindex(Xc.columns)=="T2D").values,"músculo")
    FALTAN=[g for g in ("GSE159984",) if not glob.glob(f"{RAW}/{g}*")]
    for g in FALTAN: print(f"AVISO: falta {g} en {RAW}; la figura completa requiere correrlo en el clúster")
    if FALTAN: print(f"   -> esta corrida usa {len(COH)} de 5 cohortes; los z de conjunto publicados\n"
                     "      incluyen clock output enriquecido en GSE25462 (z 4,3) y GSE78721 (z 2,5)")
    rows=[],; rows=[]; sg=[]; se=[]
    for name,(X,grp,tis) in COH.items():
        res=moderated_t(X,grp); res.to_csv(f"{OUT}/fig1_de_{name.split(' ')[0]}_{tis}.tsv",sep="\t")
        n01=int((res.FDR<0.1).sum())
        rows.append(dict(cohorte=name,tejido=tis,casos=int(grp.sum()),controles=int((~grp).sum()),
                         genes_FDR01=n01,t_mediana=float(res.t.median())))
        for setname in ("Clock output","Glucocorticoid","Insulin core","Oxidative","Core oscillator"):
            z,k=set_z(res,SETS[setname]["genes"]); se.append(dict(cohorte=name,tejido=tis,conjunto=setname,n_genes=k,z=z))
        for setname in ("Clock output","Glucocorticoid"):
            for x in SETS[setname]["genes"]:
                r=resolve(x,res.index)
                if r: sg.append(dict(cohorte=name,tejido=tis,conjunto=setname,gen=x,
                                     logFC=res.logFC[r],t=res.t[r],P=res.P[r],FDR=res.FDR[r]))
    D=pd.DataFrame(rows); S=pd.DataFrame(sg); Z=pd.DataFrame(se)
    D.to_csv(f"{OUT}/fig1_de_per_cohort.tsv",sep="\t",index=False)
    S.to_csv(f"{OUT}/fig1_signature_genes.tsv",sep="\t",index=False)
    Z.to_csv(f"{OUT}/fig1_set_enrichment.tsv",sep="\t",index=False)
    pd.set_option("display.width",200)
    print("\n== genes a FDR<0.1 por cohorte =="); print(D.to_string(index=False))
    print("\n== los 18 genes de las dos firmas: ¿alguno pasa FDR 0,1? ==")
    print(f"   total de pruebas: {len(S)} | con FDR<0.1: {int((S.FDR<0.1).sum())} | con P<0.05 sin corregir: {int((S.P<0.05).sum())}")
    print(S.sort_values("P").head(10).round(3).to_string(index=False))
    print("\n== enriquecimiento de conjunto (z frente a conjuntos del mismo tamaño) ==")
    print(Z.pivot_table(index=["cohorte","tejido"],columns="conjunto",values="z").round(2).to_string())
