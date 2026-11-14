"""Control externo: un retraso de 10 h del periodo de sueno, en laboratorio, con horas registradas
(GSE107537, celulas mononucleares de sangre, 8 participantes; Kervezee et al. 2018).

DOS HERRAMIENTAS, CADA UNA DONDE CORRESPONDE. La amplitud de un ritmo se puede medir de dos maneras
y el articulo usa las dos, segun lo que el diseno permite:

  - COSINOR, cuando hay una serie con horas registradas en cada persona. Es la medida directa: se
    ajusta una sinusoide a cada gen sobre las horas reales. Se usa en este control, y para obtener
    la amplitud medida contra la que se valida el otro estimador (miotubos y musculo in vivo).
  - El ESTIMADOR DE MUESTRA UNICA, cuando hay una sola biopsia por persona y ninguna hora. Infiere
    la amplitud de la posicion de la muestra en el ciclo que forman los genes del reloj en la
    cohorte. Se usa en las cohortes transversales, y solo en tejidos donde esos genes forman un
    ciclo reconocible (control de antifase BMAL1-DBP).

Este control es del primer tipo: tiene horas, y la pregunta es si el conjunto de genes de salida del
reloj lleva el efecto del desajuste. La respuesta sale del cosinor. El estimador de muestra unica se
aplica ademas al mismo diseno para dejar escrito donde esta su limite: la sangre no supera el
control de antifase, de modo que alli no es la herramienta adecuada, y no reproduce el efecto.

Salida: results/new/GSE107537_shift_response.tsv, una fila por participante con las dos medidas.
"""

import os as _os, sys as _sys
_RAW = _os.environ.get("T2D_RAW", "data/raw"); _OUT = _os.environ.get("T2D_OUT", "results")   # las mismas rutas que run_all.sh
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../lib"))
import gzip, io, re, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings("ignore", category=FutureWarning); warnings.filterwarnings("ignore", category=DeprecationWarning)
from scipy import stats
from gene_sets import SETS
# anotacion: TC id -> simbolo, del campo gene_assignment
ann={}
with open(f"{_RAW}/GPL24324-69524.txt", errors="ignore") as fh:
    for l in fh:
        if l.startswith("#") or l.startswith("ID\t"): continue
        p=l.rstrip("\r\n").split("\t")
        if len(p)<10 or not p[0]: continue
        m=re.search(r'\(([A-Za-z0-9orf\-]{2,15})\),', p[9]) if len(p)>9 else None
        if m: ann[p[0]]=m.group(1)
print("sondas anotadas:",len(ann))
# matriz
with gzip.open(f"{_RAW}/GSE107537_series_matrix.txt.gz","rt",errors="ignore") as fh: L=fh.read().split("\n")
i=[k for k,x in enumerate(L) if x.startswith("!series_matrix_table_begin")][0]
j=[k for k,x in enumerate(L) if x.startswith("!series_matrix_table_end")][0]
X=pd.read_csv(io.StringIO("\n".join(L[i+1:j])),sep="\t",index_col=0)
X.columns=[c.strip('"') for c in X.columns]; X.index=[str(k).strip('"') for k in X.index]
print("solapamiento con la anotacion:", len(set(X.index)&set(ann)))
X=X.loc[[k for k in X.index if k in ann]]; X.index=[ann[k] for k in X.index]
X=X.groupby(level=0).mean()
meta={}
for l in L:
    if l.startswith("!Sample_geo_accession"): meta["gsm"]=[x.strip('"') for x in l.rstrip().split("\t")[1:]]
    elif l.startswith("!Sample_characteristics"): meta.setdefault("ch",[]).append([x.strip('"') for x in l.rstrip().split("\t")[1:]])
M=pd.DataFrame({"gsm":meta["gsm"]})
for k,c in enumerate(meta["ch"]): M[f"c{k}"]=c
M["subject"]=M.c0.str.replace("subject: ",""); M["cond"]=M.c1.str.replace("condition: ","")
M["t"]=M.c2.str.replace("relclocktime: ","").astype(float)
M=M.set_index("gsm").reindex(X.columns)
clk=[g for g in SETS["Clock output"]["genes"] if g in X.index]
print(f"\ngenes del reloj medidos: {clk}")
def cosinor(y,t,per=24.0):
    w=2*np.pi/per; D=np.column_stack([np.ones(len(t)),np.cos(w*t),np.sin(w*t)])
    b,*_=np.linalg.lstsq(D,y,rcond=None); return np.hypot(b[1],b[2]), (np.arctan2(b[2],b[1])%(2*np.pi))*24/(2*np.pi)
print(f"\n{'sujeto':8s} {'dia: amp':>9s} {'noche: amp':>11s} {'cambio %':>9s} | {'fase dia':>9s} {'fase noche':>11s}")
rows=[]
for s,q in M.groupby("subject"):
    r={}
    for c in ("DayOrientedSchedule","NightOrientedSchedule"):
        qq=q[q.cond==c]
        if len(qq)<5: continue
        a=[cosinor(np.asarray(X.loc[g,qq.index],float), qq.t.values) for g in clk]
        r[c]=(np.mean([x[0] for x in a]), np.angle(np.mean([np.exp(1j*x[1]*2*np.pi/24) for x in a]))*24/(2*np.pi)%24)
    if len(r)==2:
        d,n=r["DayOrientedSchedule"],r["NightOrientedSchedule"]
        ch=(n[0]-d[0])/d[0]*100
        print(f"{s:8s} {d[0]:9.3f} {n[0]:11.3f} {ch:+8.1f}% | {d[1]:9.1f} {n[1]:11.1f}")
        rows.append(dict(subject=s,amp_day=d[0],amp_night=n[0],pct=ch,phase_day=d[1],phase_night=n[1]))
R=pd.DataFrame(rows)

# --- el estimador de muestra unica bajo el mismo diseno -----------------------------------------
import phase as ph
Xc=X.sub(X.mean(0),axis=1)                                   # centrado por muestra, como en el resto
ref=ph.fit_reference(Xc)                                     # ordenamiento libre sobre el panel del reloj
M["amp_inferred"]=[ph.infer_from_reference(Xc[g],ref)[1] for g in M.index]
Z=Xc.loc[ref["genes"]]; Z=((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T
clock_ok,sep=ph.is_clock(Xc,ph.unsup(Z.values),genes=ref["genes"])
inf=M.groupby(["subject","cond"]).amp_inferred.mean().unstack()
R["amp_inferred_day"]=R.subject.map(inf["DayOrientedSchedule"]); R["amp_inferred_night"]=R.subject.map(inf["NightOrientedSchedule"])
R["pct_inferred"]=(R.amp_inferred_night-R.amp_inferred_day)/R.amp_inferred_day.abs()*100
print(f"\nn = {len(R)} sujetos")
print(f"  amplitud: mediana del cambio {R.pct.median():+.1f}%, {(R.pct<0).sum()}/{len(R)} bajan, Wilcoxon P = {stats.wilcoxon(R.amp_night,R.amp_day).pvalue:.4f}")
d=np.angle(np.exp(1j*(R.phase_night-R.phase_day)*2*np.pi/24))*24/(2*np.pi)
print(f"  fase: desplazamiento mediano {np.median(d):+.1f} h, rango {d.min():+.1f} a {d.max():+.1f} h")
pw=stats.wilcoxon(R.amp_inferred_night,R.amp_inferred_day).pvalue
rr=stats.spearmanr(R.pct,R.pct_inferred)
print(f"\nestimador de muestra unica, mismo diseno ({len(ref['genes'])} genes del panel):")
print("  cambio por participante (%): "+", ".join(f"{v:+.0f}" for v in R.pct_inferred.sort_values()))
print(f"  mediana del cambio {R.pct_inferred.median():+.1f}%, {(R.pct_inferred<0).sum()}/{len(R)} bajan, Wilcoxon P = {pw:.4f}")
print(f"  concordancia con el cambio por cosinor entre participantes: rho {rr.statistic:+.2f} (P {rr.pvalue:.2f})")
print(f"  control de antifase BMAL1-DBP en este tejido: {'superado' if clock_ok else 'NO superado'} (separacion {sep:.1f} h; se exige 8 h o mas)")
if not clock_ok:
    print("  -> en este tejido los genes del reloj no forman un ciclo reconocible, asi que el estimador de muestra\n"
          "     unica no es la herramienta: la medida que corresponde a este diseno, con horas, es el cosinor.")
elif pw>=0.05 or (R.pct_inferred<0).mean()<0.75:
    print("  -> el estimador de muestra unica no reproduce aqui el efecto que mide el cosinor.")
R["antiphase_control_passed"]=bool(clock_ok); R["bmal1_dbp_separation_h"]=sep
R.to_csv(f"{_OUT}/new/GSE107537_shift_response.tsv",sep="\t",index=False)
