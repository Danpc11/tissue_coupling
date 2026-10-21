import pandas as pd, numpy as np, sys, warnings; warnings.filterwarnings("ignore")
from scipy import stats
print("¿La mayor dispersión en resistentes es ruido de medida? (A6 -> 43.4)\n")
L=pd.read_csv("results/new/GSE129843_loo_validation.tsv",sep="\t")
# ruido del estimador DENTRO de la misma serie (misma persona, misma condicion): sd de las inferencias
w=L.groupby(["ind","cond"]).amp_inf.std()
print(f"1) Ruido del estimador dentro de una misma serie (GSE129843, 6 biopsias/persona):")
print(f"   sd mediana dentro de serie = {w.median():.3f}   rango {w.min():.3f}-{w.max():.3f}")
b=L.groupby(["ind","cond"]).amp_inf.mean()
print(f"   sd ENTRE series (señal + ruido)          = {b.std():.3f}")
icc=max(0,(b.var()-w.pow(2).mean()/6))/b.var()
print(f"   -> fraccion de la varianza entre personas que NO es ruido: {icc*100:.0f}%")
print(f"   -> con UNA sola muestra, la varianza observada entre personas se infla en")
print(f"      {w.pow(2).mean():.3f} unidades (varianza de medida) sobre la varianza real\n")
D=pd.read_csv("results/new/amplitude_phase_distributions_two_cohorts.tsv",sep="\t",index_col=0)
D["coh"]=D.cohort.str.replace(" músculo basal","",regex=False).str.replace(" músculo","",regex=False).str.strip()
print("2) ¿La diferencia de dispersión sobrevive al descontar varianza de medida constante?")
print("   (si el ruido es igual en ambos grupos, restarlo no cambia la RAZON de varianzas,")
print("    pero sí la acerca a 1: var_real = var_obs - var_ruido)")
q=D[D.coh=="GSE13070"]
for col in ("amp_free","amp_clock"):
    gs={k:v[col].dropna().values for k,v in q.groupby("group")}
    a,b2=gs["resistant"],gs["sensitive"]
    va,vb=a.var(ddof=1),b2.var(ddof=1)
    print(f"   {col}: var resistentes {va:.3f}, sensibles {vb:.3f}, razon {va/vb:.2f}")
    for vn in (0.0,0.25*vb,0.5*vb,0.75*vb):
        if vb-vn>0: print(f"      si la varianza de medida fuese {vn/vb*100:3.0f}% de la de sensibles -> razon real {(va-vn)/(vb-vn):.2f}")
print("\n3) Prueba directa: ¿difiere el RUIDO del estimador entre grupos?")
print("   Requiere varias muestras por persona en GSE13070; no las hay (una biopsia basal).")
print("   -> no se puede descartar con los datos de esa cohorte.")
