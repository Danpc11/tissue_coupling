"""Metaanalisis de efectos aleatorios del eje receptor-reloj.

DerSimonian-Laird sobre Fisher-z. Una estimacion por cohorte y tejido, residualizada frente a las
covariables de cada cohorte. El modelo que respeta que varias estimaciones comparten donantes es
el de tres niveles (13_meta_three_level.py), que lee la tabla que escribe este paso.

LAS TRECE ESTIMACIONES SON OBLIGATORIAS. La lista esta en MANIFEST. Si falta alguna, porque falta
su matriz, el paso lo dice y termina con error en vez de agrupar las que haya: un resultado con
once estimaciones no es el del articulo. T2D_ALLOW_PARTIAL=1 permite seguir, y la tabla lo anota.

ERROR ESTANDAR. `se` es la aproximacion de Fisher, 1/sqrt(n - 3 - covariables), que es la de una
correlacion de Pearson. Aqui se aplica a una de Spearman sobre residuos, asi que la tabla lleva
tambien `se_boot`, el error estandar por bootstrap de individuos (se remuestrean personas y se
repite la residualizacion), para comprobar que el resultado no depende de esa aproximacion.

GRUPO DIAGNOSTICO. Las estimaciones mezclan casos y controles, porque la pregunta es si receptor y
salida del reloj varian juntos entre personas. `rho_adjusted_for_group` es la misma correlacion con
el grupo (diagnostico, o brazo del estudio) como covariable, donde la cohorte lo tiene: dice si la
asociacion existe DENTRO de los grupos o es solo la diferencia entre ellos.

CENTRADO. Las estimaciones usan expresion centrada por muestra. `rho_uncentred` es la misma
correlacion sin centrar, para que la sensibilidad a esa decision quede escrita en la tabla.

ESPECIFICIDAD. La misma correlacion parcial se calcula para todos los genes de cada matriz y se
guarda (meta_all_genes_rho.pkl); `percentile_among_all_genes` dice donde queda NR3C1 entre ellos, y
`rho_given_the_other` es la correlacion de cada predictor con la salida del reloj una vez quitado
lo que explica el otro (NR3C1 dada la senal glucocorticoide, y al reves).

Salidas: results/new/meta_per_cohort.tsv, meta_pooled.tsv y meta_all_genes_rho.pkl
"""
import os as _os, sys as _sys
_RAW = _os.environ.get("T2D_RAW", "data/raw"); _OUT = _os.environ.get("T2D_OUT", "results")   # las mismas rutas que run_all.sh
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../lib"))
import pandas as pd, numpy as np, sys, gzip, importlib.util, os
from scipy import stats
OUT=os.environ.get("T2D_OUT","results"); RAW=os.environ.get("T2D_RAW","data/raw"); from gene_sets import SETS, resolve; import islets
def sc(X,g):
    gg=[r for r in (resolve(x,X.index) for x in g) if r]; Z=X.loc[gg]; return ((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T.mean(0)
def design(n,cv):
    """Matriz de covariables con intercepto, comprobada: sin valores no finitos y con muestras de sobra."""
    D=np.column_stack([np.ones(n)]+cv) if cv else np.ones((n,1))
    if not np.isfinite(D).all(): raise ValueError("hay covariables con valores ausentes o no finitos")
    if n-np.linalg.matrix_rank(D)<10: raise ValueError(f"quedan menos de 10 grados de libertad ({n} muestras, rango {np.linalg.matrix_rank(D)})")
    return D
def pres(y,cv,D=None):
    D=design(len(y),cv) if D is None else D
    y=np.asarray(y,float)
    if not np.isfinite(y).all(): raise ValueError("la variable a residualizar tiene valores no finitos")
    b,*_=np.linalg.lstsq(D,y,rcond=None); return y-D@b
NBOOT=int(os.environ.get("T2D_BOOT_META",500))
MANIFEST=[("AAGMEx","adiposo"),("AAGMEx","músculo"),("GSE41168","adiposo"),("GSE41168","músculo"),("GSE13070","adiposo"),("GSE13070","músculo"),
          ("GTEx","adiposo sc"),("GTEx","adiposo visceral"),("GTEx","músculo"),("GSE78721","adiposo"),("GSE25462","músculo"),
          ("GSE159984","islotes (donantes)"),("GSE159984","islotes (cultivos)")]
CL=SETS["Clock output"]["genes"]; GC=SETS["Glucocorticoid"]["genes"]
PH=pd.read_pickle(f"{_OUT}/new/GTEx_pheno.pkl"); agemid={'20-29':25,'30-39':35,'40-49':45,'50-59':55,'60-69':65,'70-79':75}
def gtexcov(X):
    p=PH.reindex(X.columns); H=pd.get_dummies(p.DTHHRDY.fillna(-1).astype(int),prefix="h",drop_first=True).values.astype(float)
    cv=[(p.SEX==1).astype(float).values,p.AGE.map(agemid).values.astype(float)]+[H[:,k] for k in range(H.shape[1])]
    return cv
A=pd.read_pickle(f"{_OUT}/new/AAGMEx_adipose_log2.pkl"); M=pd.read_pickle(f"{_OUT}/new/AAGMEx_muscle_log2.pkl")
with gzip.open(f"{_RAW}/GSE95674_series_matrix.txt.gz","rt") as f:
    for l in f:
        if l.startswith("!Sample_title"): t=[x.strip('"') for x in l.rstrip().split("\t")[1:]]
        if l.startswith("!Sample_characteristics") and "gender" in l: sx=[x.strip('"').split(": ")[1] for x in l.rstrip().split("\t")[1:]]
        if l.startswith("!Sample_characteristics") and "age" in l: ag=[float(x.strip('"').split(": ")[1]) for x in l.rstrip().split("\t")[1:]]
cov=pd.DataFrame({"sex":sx,"age":ag},index=t).reindex(A.columns); aag=[(cov.sex=="M").astype(float).values,cov.age.values]
Ea=pd.read_pickle(f"{_OUT}/new/GSE41168_adipose.pkl"); Em=pd.read_pickle(f"{_OUT}/new/GSE41168_muscle.pkl")
sp=importlib.util.spec_from_file_location("c5",os.path.join(os.path.dirname(os.path.abspath(__file__)),"05_clamp_within_person.py")); c5=importlib.util.module_from_spec(sp); sp.loader.exec_module(c5)
E_sym,d=c5.load()
# GSE13070 es GPL570, como GSE25462 y GSE41168: se anota igual que ellas, por Entrez y con la sonda de
# mayor expresion media (python/lib/gpl570.py). Hasta v2.7.0 entraba anotada por el simbolo de la
# plataforma y promediando sondas, que es la regla del analisis del clamp, y sus dos estimaciones
# salian 0.1 mas altas. T2D_GSE13070_SYMBOLS=1 recupera aquella anotacion.
_f13=f"{OUT}/new/GSE13070_entrez.pkl"
E=E_sym if (os.environ.get("T2D_GSE13070_SYMBOLS")=="1" or not os.path.exists(_f13)) else pd.read_pickle(_f13)
mus13=d[(d.tissue=="muscle")&(d.tzd=="pre")&(d.clamp=="pre")].drop_duplicates("subject").set_index("subject")
adi13=d[(d.tissue=="adipose")&(d.tzd=="pre")].drop_duplicates("subject").set_index("subject")
G1=pd.read_pickle(f"{_OUT}/new/GTEx_adipose_sc.pkl"); G2=pd.read_pickle(f"{_OUT}/new/GTEx_adipose_visc.pkl"); G3=pd.read_pickle(f"{_OUT}/new/GTEx_muscle.pkl")
rows=[]; ALLG={}
def add(coh,tis,X0,cv,x_lab="NR3C1",grp=None):
    X=X0.sub(X0.mean(0),axis=1) if os.environ.get("T2D_CENTER","1")!="0" else X0   # centrado por muestra (Methods, "Sample centring"); T2D_CENTER=0 lo desactiva
    n=X.shape[1]; D=design(n,cv)
    raw={"NR3C1~reloj":np.asarray(X.loc[x_lab],float),"GC~reloj":np.asarray(sc(X,GC),float)}; yraw=np.asarray(sc(X,CL),float)
    unc={"NR3C1~reloj":np.asarray(X0.loc[x_lab],float),"GC~reloj":np.asarray(sc(X0,GC),float)}; yunc=np.asarray(sc(X0,CL),float)
    y=pres(yraw,cv,D); rng=np.random.default_rng(0)
    # la misma correlacion parcial para TODOS los genes: es el fondo contra el que se mide NR3C1
    V=X.values.T.astype(float); V=V-D@np.linalg.lstsq(D,V,rcond=None)[0]
    rk=np.apply_along_axis(stats.rankdata,0,V); sd=rk.std(0); sd[sd==0]=np.nan
    ry=stats.rankdata(y); allrho=pd.Series(((rk-rk.mean(0))/sd*((ry-ry.mean())/ry.std())[:,None]).mean(0),index=X.index).dropna()
    ALLG[f"{coh} | {tis}"]=allrho.astype("float32"); del V,rk
    inscore={r for r in (resolve(x,X.index) for x in CL) if r}
    bg=allrho.drop([g for g in allrho.index if g in inscore or g==x_lab])
    pct_all=float((bg<allrho[x_lab]).mean()*100)
    def partial(a,b,c):                                   # rho de a con b quitando a los dos lo que explica c
        Dc=np.column_stack([D,c]); return stats.spearmanr(pres(a,cv,Dc),pres(b,cv,Dc)).statistic
    part={"NR3C1~reloj":partial(raw["NR3C1~reloj"],yraw,raw["GC~reloj"]),"GC~reloj":partial(raw["GC~reloj"],yraw,raw["NR3C1~reloj"])}
    idx=[rng.integers(0,n,n) for _ in range(NBOOT)]          # las mismas remuestras para los dos contrastes
    for lab in raw:
        r=stats.spearmanr(pres(raw[lab],cv,D),y).statistic
        z=np.arctanh(np.clip(r,-.999,.999)); se=1/np.sqrt(n-3-len(cv))
        zb=[]
        for i in idx:                                         # bootstrap de individuos: se repite la residualizacion
            Di=D[i]
            if np.linalg.matrix_rank(Di)<Di.shape[1] and Di.shape[1]>1: continue
            rb=stats.spearmanr(pres(raw[lab][i],cv,Di),pres(yraw[i],cv,Di)).statistic
            zb.append(np.arctanh(np.clip(rb,-.999,.999)))
        ru=stats.spearmanr(pres(unc[lab],cv,D),pres(yunc,cv,D)).statistic
        # con el grupo diagnostico o experimental como covariable: es la asociacion solo cosa de casos frente a controles?
        rg=np.nan
        if grp is not None:
            Gd=pd.get_dummies(pd.Series(np.asarray(grp)).astype(str),drop_first=True).values.astype(float)
            Dg=np.column_stack([D,Gd])
            if n-np.linalg.matrix_rank(Dg)>=10: rg=stats.spearmanr(pres(raw[lab],cv,Dg),pres(yraw,cv,Dg)).statistic
        rows.append(dict(contraste=lab,cohorte=coh,tejido=tis,n=n,rho=r,z=z,se=se,
                         se_boot=float(np.std(zb,ddof=1)) if len(zb)>50 else np.nan,rho_uncentred=ru,
                         rho_adjusted_for_group=rg,group_levels=(pd.Series(np.asarray(grp)).nunique() if grp is not None else 0),
                         rho_given_the_other=part[lab],percentile_among_all_genes=pct_all if lab=="NR3C1~reloj" else np.nan))
add("AAGMEx","adiposo",A,aag); add("AAGMEx","músculo",M,aag)
_p41=pd.read_pickle(f"{OUT}/new/GSE41168_meta.pkl"); _p41=_p41[_p41.treatment=="before"]
_g41=lambda X,key: _p41[_p41.tissue.str.contains(key)].drop_duplicates("subject").set_index("subject").group.reindex(X.columns).values
add("GSE41168","adiposo",Ea,[],grp=_g41(Ea,"adipose")); add("GSE41168","músculo",Em,[],grp=_g41(Em,"muscle"))
add("GSE13070","adiposo",E[list(adi13.gsm)],[],grp=adi13.group.values); add("GSE13070","músculo",E[list(mus13.gsm)],[],grp=mus13.group.values)
add("GTEx","adiposo sc",G1,gtexcov(G1)); add("GTEx","adiposo visceral",G2,gtexcov(G2)); add("GTEx","músculo",G3,gtexcov(G3))
# --- las cuatro estimaciones restantes de la Fig. 2a (13 en total) ---
if os.path.exists(f"{OUT}/new/GSE78721_log2.pkl"):
    X=pd.read_pickle(f"{OUT}/new/GSE78721_log2.pkl"); md=pd.read_pickle(f"{OUT}/new/GSE78721_meta.pkl").set_index("gsm").reindex(X.columns)
    keep=list(md.index[md.one_per_person]); X=X[keep]; md=md.reindex(keep)   # una muestra por persona (95): extremidad o subcutaneo
    add("GSE78721","adiposo",X,[(md.sex=="M").astype(float).values,(md.depot=="Subcutaneous").astype(float).values],grp=md.dx.values)
if os.path.exists(f"{OUT}/new/GSE25462_entrez.pkl"):
    X=pd.read_pickle(f"{OUT}/new/GSE25462_entrez.pkl"); md=pd.read_pickle(f"{OUT}/new/GSE25462_meta.pkl").set_index("gsm").reindex(X.columns)
    add("GSE25462","músculo",X,[(md.sex=="male").astype(float).values,md.age.fillna(md.age.median()).values.astype(float)],grp=md.fh.values)
if os.path.exists(f"{RAW}/GSE159984_raw_counts.tsv.gz"):
    don,cul=islets.groups(RAW)
    _dx=islets.diagnosis(RAW)                                # diabetes tipo 2 o no, de los donantes
    for lab,cols in (("islotes (donantes)",don),("islotes (cultivos)",cul)):
        L=islets.log2cpm_symbols(f"{RAW}/GSE159984_raw_counts.tsv.gz",cols,f"{RAW}/Human_GRCh38_p13_annot_tsv.gz")
        if L.shape[1]<10: print(f"  GSE159984 {lab}: {L.shape[1]} muestras en el fichero de cuentas; se omite"); continue
        add("GSE159984",lab,L,[],grp=(_dx.reindex(L.columns).values if lab.endswith("(donantes)") and _dx.reindex(L.columns).notna().all() else None))
R=pd.DataFrame(rows)
pd.DataFrame(ALLG).to_pickle(f"{OUT}/new/meta_all_genes_rho.pkl")    # gen x estimacion: lo lee 30_receptor_specificity.py
have=set(zip(R.cohorte,R.tejido)); missing=[m for m in MANIFEST if m not in have]
R["complete"]=not missing
if missing:
    print(f"\nFALTAN {len(missing)} de las {len(MANIFEST)} estimaciones del articulo: "+"; ".join(f"{c} {t}" for c,t in missing))
    print("El agrupado de abajo NO es el del articulo."+("" if os.environ.get("T2D_ALLOW_PARTIAL")=="1" else " T2D_ALLOW_PARTIAL=1 permite continuar."))
    if os.environ.get("T2D_ALLOW_PARTIAL")!="1":
        R.to_csv(f"{OUT}/new/meta_per_cohort.tsv",sep="\t",index=False); sys.exit(1)
def dl(sub):
    w=1/sub.se**2; m=(w*sub.z).sum()/w.sum(); Q=(w*(sub.z-m)**2).sum(); k=len(sub)
    tau2=max(0,(Q-(k-1))/(w.sum()-(w**2).sum()/w.sum())); w2=1/(sub.se**2+tau2); m2=(w2*sub.z).sum()/w2.sum()
    se2=np.sqrt(1/w2.sum()); I2=max(0,(Q-(k-1))/Q*100) if Q>0 else 0
    return np.tanh(m2),np.tanh(m2-1.96*se2),np.tanh(m2+1.96*se2),Q,k,I2,2*stats.norm.sf(abs(m2/se2))
out=[]
for lab in R.contraste.unique():
    for tej in ("todos","adiposo","músculo","islotes"):
        sub=R[(R.contraste==lab)&(R.tejido.str.contains(tej) if tej!="todos" else True)]
        if len(sub)<2: continue
        rr,lo,hi,Q,k,I2,p=dl(sub)
        print(f"{lab:12s} {tej:8s} k={k}: ρ combinado {rr:+.3f} (IC95 {lo:+.3f} a {hi:+.3f}), P {p:.1e}, Q={Q:.1f}, I²={I2:.0f}%")
        out.append(dict(contraste=lab,grupo=tej,k=k,rho=rr,lo=lo,hi=hi,p=p,Q=Q,I2=I2))
R.to_csv(f"{OUT}/new/meta_per_cohort.tsv",sep="\t",index=False); pd.DataFrame(out).to_csv(f"{OUT}/new/meta_pooled.tsv",sep="\t",index=False)
print("\nPor cohorte:"); print(R.round(3).to_string(index=False))
