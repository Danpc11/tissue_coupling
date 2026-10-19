"""Figura de mecanismo del envio anterior; conservada como diagnostico.
La version del manuscrito es paper_figures.py (F3 y F4).
"""
import pandas as pd, numpy as np, sys, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
sys.path.insert(0,"python/lib"); from gene_sets import SETS, resolve
BLUE,ORANGE,GREEN,GREY,PURPLE="#1f77b4","#d95f02","#2a9d8f","#999999","#6a3d9a"
R=pd.read_csv("results/new/meta_per_cohort.tsv",sep="\t")
def dl(sub):
    w=1/sub.se**2; m=(w*sub.z).sum()/w.sum(); Q=(w*(sub.z-m)**2).sum(); k=len(sub)
    tau2=max(0,(Q-(k-1))/(w.sum()-(w**2).sum()/w.sum())); w2=1/(sub.se**2+tau2)
    m2=(w2*sub.z).sum()/w2.sum(); se2=np.sqrt(1/w2.sum()); I2=max(0,(Q-(k-1))/Q*100) if Q>0 else 0
    return np.tanh(m2),np.tanh(m2-1.96*se2),np.tanh(m2+1.96*se2),I2,2*(1-stats.norm.cdf(abs(m2/se2)))
fig=plt.figure(figsize=(16.5,9.8))
gs=fig.add_gridspec(2,3,height_ratios=[1.15,1],hspace=.55,wspace=.42)
# --- a) forest NR3C1
ax=fig.add_subplot(gs[0,0])
sub=R[R.contraste=="NR3C1~reloj"].copy()
sub["lo"]=np.tanh(sub.z-1.96*sub.se); sub["hi"]=np.tanh(sub.z+1.96*sub.se)
order=["adiposo","adiposo sc","adiposo visceral","músculo"]
sub["ord"]=sub.tejido.apply(lambda t: 0 if "adiposo" in t else 1)
sub=sub.sort_values(["ord","n"],ascending=[True,False]).reset_index(drop=True)
y=np.arange(len(sub))
for i,r in sub.iterrows():
    c=GREEN if "adiposo" in r.tejido else ORANGE
    ax.plot([r.lo,r.hi],[i,i],color=c,lw=1.4); ax.scatter(r.rho,i,s=18+r.n/12,color=c,zorder=3)
ad=dl(sub[sub["ord"]==0]); mu=dl(sub[sub["ord"]==1]); allr=dl(sub)
for j,(lab,v,c) in enumerate(((f"adiposo (k=5)",ad,GREEN),(f"músculo (k=4)",mu,ORANGE),("global (k=9)",allr,"k"))):
    i=len(sub)+j+.6; ax.plot([v[1],v[2]],[i,i],color=c,lw=2.5); ax.scatter(v[0],i,marker="D",s=55,color=c,zorder=3)
ax.set_yticks(list(y)+[len(sub)+.6,len(sub)+1.6,len(sub)+2.6])
ax.set_yticklabels([f"{r.cohorte} · {r.tejido} (n={r.n})" for _,r in sub.iterrows()]+
                   [f"adiposo: {ad[0]:.2f} (I²={ad[3]:.0f}%)",f"músculo: {mu[0]:.2f} (I²={mu[3]:.0f}%)",f"global: {allr[0]:.2f}, P={allr[4]:.0e}"],fontsize=8)
ax.invert_yaxis(); ax.axvline(0,color="k",lw=.8); ax.axvline(allr[0],color="k",ls=":",lw=.8)
ax.set(xlabel="ρ receptor (NR3C1) ~ reloj",xlim=(-0.2,0.95),title="a  El reloj sigue al receptor\nen nueve cohortes-tejido")
# --- b) forest GC, vivos vs post mortem
ax=fig.add_subplot(gs[0,1])
s2=R[R.contraste=="GC~reloj"].copy(); s2["lo"]=np.tanh(s2.z-1.96*s2.se); s2["hi"]=np.tanh(s2.z+1.96*s2.se)
s2["pm"]=(s2.cohorte=="GTEx").astype(int); s2=s2.sort_values(["pm","n"],ascending=[True,False]).reset_index(drop=True)
for i,r in s2.iterrows():
    c=PURPLE if r.pm==0 else GREY
    ax.plot([r.lo,r.hi],[i,i],color=c,lw=1.4); ax.scatter(r.rho,i,s=18+r.n/12,color=c,zorder=3)
vv=dl(s2[s2.pm==0]); pm=dl(s2[s2.pm==1]); vad=dl(s2[(s2.pm==0)&(s2.tejido.str.contains("adiposo"))])
for j,(lab,v,c) in enumerate(((f"vivos (k=6)",vv,PURPLE),(f"vivos, adiposo (k=3)",vad,"#4a1486"),(f"GTEx post mortem (k=3)",pm,GREY))):
    i=len(s2)+j+.6; ax.plot([v[1],v[2]],[i,i],color=c,lw=2.5); ax.scatter(v[0],i,marker="D",s=55,color=c,zorder=3)
ax.set_yticks(list(np.arange(len(s2)))+[len(s2)+.6,len(s2)+1.6,len(s2)+2.6])
ax.set_yticklabels([f"{r.cohorte} · {r.tejido} (n={r.n})" for _,r in s2.iterrows()]+
                   [f"vivos: {vv[0]:.2f}, P={vv[4]:.0e}",f"vivos adiposo: {vad[0]:.2f} (I²={vad[3]:.0f}%)",f"post mortem: {pm[0]:.2f}, P={pm[4]:.2f}"],fontsize=8)
ax.invert_yaxis(); ax.axvline(0,color="k",lw=.8)
ax.set(xlabel="ρ señalización glucocorticoide ~ reloj",xlim=(-0.55,0.85),title="b  La señalización se mide en vivos,\nno en tejido post mortem")
# --- c) Hallmark
ax=fig.add_subplot(gs[0,2])
T=pd.read_csv("results/new/AAGMEx_NR3C1_enrichment_extra.tsv",sep="\t"); H=T[T.lib=="h_all"].copy()
H["lab"]=H.set.str.replace("HALLMARK_","").str.replace("_"," ").str.capitalize()
H=pd.concat([H.nlargest(4,"z"),H.nsmallest(14,"z")]).sort_values("z")
ax.barh(range(len(H)),H.z,color=[BLUE if v>0 else ORANGE for v in H.z])
ax.set_yticks(range(len(H))); ax.set_yticklabels(H.lab,fontsize=7); ax.axvline(0,color="k",lw=.8)
ax.set(xlabel="z frente a conjuntos del mismo tamaño",title="c  El eje del receptor es antiinflamatorio y va\ncontra OXPHOS, adipogénesis y mTORC1")
# --- d) TZD no repara la firma
ax=fig.add_subplot(gs[1,0])
sig=pd.read_csv("/mnt/user-data/outputs/new_datasets/signatures_log2_common_genes.tsv",sep="\t",index_col=0)
genes=["PER3","DBP","HLF","TEF","CIART","NR1D1","NR1D2","PER2","BMAL1"]
cols={"IR adiposo":"IR adipose (GSE13070, 35 vs 6)","Cushing":"Cushing adipose (GSE66446)","hidrocortisona":"hydrocortisone adipose (GSE148640)","TZD 3 meses":"TZD 3mo adipose (GSE13070)"}
Tm=sig.loc[[g for g in genes if g in sig.index],list(cols.values())]; Tm.columns=list(cols.keys())
x=np.arange(len(Tm)); w=.2
for k,(c,col) in enumerate(zip(Tm.columns,[ORANGE,"#6a3d9a",BLUE,GREEN])):
    ax.bar(x+(k-1.5)*w,Tm[c],w,color=col,label=c)
ax.axhline(0,color="k",lw=.8); ax.set_xticks(x); ax.set_xticklabels(Tm.index,rotation=45,fontsize=8,ha="right")
r1=stats.spearmanr(Tm["IR adiposo"],Tm["Cushing"]); r2=stats.spearmanr(Tm["IR adiposo"],Tm["TZD 3 meses"])
ax.set(ylabel="log2",title=f"d  La firma del adiposo resistente es la del cortisol;\nla TZD no la revierte\nρ con Cushing {r1.statistic:.2f} · con TZD {r2.statistic:.2f}")
ax.legend(fontsize=7,frameon=False)
# --- e) clamp = mañana
ax=fig.add_subplot(gs[1,1])
C=pd.read_csv("/mnt/user-data/outputs/new_datasets/time_calibration_vs_clamp.tsv",sep="\t",index_col=0)
gl=["NR1D1","DBP","CIART","PER1","PER3","TEF"]
x=np.arange(len(gl)); w=.2
for k,(c,col) in enumerate(zip(C.columns,[GREY,GREEN,"#999955",ORANGE])):
    ax.bar(x+(k-1.5)*w,C.loc[gl,c],w,color=col,label=c.replace("GSE129843 ","129843 ").replace("GSE108539 ","108539 ").replace("GSE13070 ","13070 "))
ax.axhline(0,color="k",lw=.8); ax.set_xticks(x); ax.set_xticklabels(gl,rotation=45,fontsize=8)
ax.set(ylabel="cambio (log2)",title="e  Tres horas de clamp = tres horas de mañana"); ax.legend(fontsize=6.5,frameon=False)
# --- f) modelo
ax=fig.add_subplot(gs[1,2]); ax.axis("off")
ax.set_title("f  Modelo",fontsize=10,loc="left")
txt=[(0.00,0.93,"TONO GLUCOCORTICOIDE DEL INDIVIDUO","k",10,"bold"),
     (0.03,0.84,"el receptor correlaciona entre tejidos del mismo donante (ρ 0,66)","k",8,"normal"),
     (0.00,0.73,"↓","k",13,"bold"),
     (0.00,0.63,"ADIPOSO","#2a9d8f",9.5,"bold"),
     (0.00,0.55,"reloj apagado","k",8.5,"normal"),
     (0.00,0.48,"PER3, DBP, HLF, TEF,","k",7.5,"normal"),
     (0.00,0.42,"NR1D1/2 ↓   BMAL1 ↑","k",7.5,"normal"),
     (0.00,0.35,"= firma de Cushing","k",7.5,"italic"),
     (0.00,0.29,"  e hidrocortisona","k",7.5,"italic"),
     (0.55,0.63,"MÚSCULO","#d95f02",9.5,"bold"),
     (0.55,0.55,"ayuno encendido","k",8.5,"normal"),
     (0.55,0.48,"TXNIP, PDK4,","k",7.5,"normal"),
     (0.55,0.42,"ANGPTL4, RRAD,","k",7.5,"normal"),
     (0.55,0.35,"G0S2 ↑","k",7.5,"normal"),
     (0.00,0.17,"INSULINA AGUDA: sin efecto sobre el reloj","#777777",9,"bold"),
     (0.03,0.10,"el clamp reproduce el paso del tiempo, no una señal","k",8,"normal"),
     (0.00,0.00,"La firma mide señal glucocorticoide efectiva, no resistencia a la insulina.","k",8.5,"bold")]
for xx,yy,s,c,fs,w_ in txt: ax.text(xx,yy,s,color=c,fontsize=fs,fontweight=("bold" if w_=="bold" else "normal"),style=("italic" if w_=="italic" else "normal"),transform=ax.transAxes)
for a in fig.axes:
    a.spines[["top","right"]].set_visible(False); a.title.set_fontsize(9.5)
fig.savefig("results/figs/F6_mecanismo.png",dpi=160,bbox_inches="tight"); fig.savefig("results/figs/F6_mecanismo.pdf",bbox_inches="tight")
print("ok")
