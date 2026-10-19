#!/usr/bin/env python3
"""Las cinco figuras del manuscrito, en la distribución fusionada.

F1  El problema y la herramienta        F2  El eje
F3  Cortisol sí, insulina no            F4  Regulación y especificidad
F5  El modelo
Requiere: 01b (Fig 1), 05 (clamp), 06 (metaanálisis), 07 (TF), model/make_*.py (Fig 5).
Los paneles de diagrama (esquema de diseño, línea de tiempo, modelo) se dejan como marco vacío
rotulado: se dibujan a mano, porque las imágenes generadas por IA no las aceptan los editores.
"""
import os,sys,warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,os.path.join(HERE,"../lib")); sys.path.insert(0,os.path.join(HERE,"../model"))
from gene_sets import SETS, resolve
OUT=os.environ.get("T2D_OUT","results"); RAW=os.environ.get("T2D_RAW","data/raw")
AD,MU,IS,GC,GR="#2a9d8f","#d95f02","#6a3d9a","#1f77b4","#999999"
def sc(X,g):
    gg=[r for r in (resolve(x,X.index) for x in g) if r]; Z=X.loc[gg]; return ((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T.mean(0)
def center(X): return X.sub(X.mean(0),axis=1)
def placeholder(ax,txt):
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.02,0.05),0.96,0.88,fill=False,ls="--",lw=1,ec="#aaaaaa",transform=ax.transAxes))
    ax.text(.5,.5,txt,ha="center",va="center",fontsize=8,color="#777777",transform=ax.transAxes,style="italic")
def tidy(fig):
    for a in fig.axes:
        a.spines[["top","right"]].set_visible(False); a.title.set_fontsize(9)
# ---------------------------------------------------------------- F1
def fig1():
    fig=plt.figure(figsize=(16.5,12.5)); gs=fig.add_gridspec(3,4,hspace=.6,wspace=.4)
    ax=fig.add_subplot(gs[0,0]); placeholder(ax,"a  Diseño\n\ncinco cohortes, dos tejidos,\nqué puede responder cada una\n\n(dibujo a mano)")
    S18=pd.read_csv(f"{OUT}/fig1_signature_genes.tsv",sep="\t")
    D=pd.read_csv(f"{OUT}/fig1_de_per_cohort.tsv",sep="\t")
    ax=fig.add_subplot(gs[0,1:3])
    cohs=list(D.sort_values("tejido").cohorte); genes=[g for g in SETS["Clock output"]["genes"]+SETS["Glucocorticoid"]["genes"]]
    genes=[g for g in genes if g in set(S18.gen)]
    for i,c in enumerate(cohs):
        q=S18[S18.cohorte==c].set_index("gen")
        for j,g in enumerate(genes):
            if g in q.index:
                t=q.loc[g,"t"]; pv=q.loc[g,"P"]
                ax.scatter(j,i,s=12+70*min(-np.log10(max(pv,1e-6))/3,1),c=[t],cmap="RdBu_r",vmin=-4,vmax=4,
                           edgecolors="k",linewidths=.3)
    ax.axvline(len(SETS["Clock output"]["genes"])-.5,color="k",lw=.8,ls="--")
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes,rotation=90,fontsize=7)
    ax.set_yticks(range(len(cohs))); ax.set_yticklabels([f"{c}\n{int(D.set_index(chr(99)+chr(111)+chr(104)+chr(111)+chr(114)+chr(116)+chr(101)).genes_FDR01[c])} genes FDR<0,1" for c in cohs],fontsize=7)
    sm=plt.cm.ScalarMappable(cmap="RdBu_r",norm=plt.Normalize(-4,4)); plt.colorbar(sm,ax=ax,label="t moderada (caso − control)",shrink=.8)
    ax.text(1.5,len(cohs)-.3,"salida del reloj",fontsize=7,ha="center"); ax.text(len(genes)-4,len(cohs)-.3,"glucocorticoide",fontsize=7,ha="center")
    ax.set(title="b  Ningún gen de las dos firmas pasa FDR 0,1 en ninguna cohorte\n(tamaño del punto: −log10 P sin corregir)",ylim=(-.6,len(cohs)-.2))
    Z=pd.read_csv(f"{OUT}/fig1_set_enrichment.tsv",sep="\t")
    p=Z.pivot_table(index="cohorte",columns="conjunto",values="z").reindex(columns=["Clock output","Glucocorticoid","Insulin core","Oxidative"]).astype(float)
    ax=fig.add_subplot(gs[0,3]); im=ax.imshow(p.values,cmap="RdBu_r",vmin=-7,vmax=7,aspect="auto")
    ax.set_xticks(range(4)); ax.set_xticklabels(["salida\nreloj","gluco-\ncorticoide","insulina","oxidativo"],fontsize=7)
    ax.set_yticks(range(len(p))); ax.set_yticklabels(p.index,fontsize=7)
    for i in range(p.shape[0]):
        for j in range(p.shape[1]):
            v=p.values[i,j]; ax.text(j,i,"—" if np.isnan(v) else f"{v:.1f}",ha="center",va="center",fontsize=7,color="w" if (not np.isnan(v) and abs(v)>4) else "k")
    plt.colorbar(im,ax=ax,label="z",shrink=.75)
    ax.set(title="c  El reloj ya sale alterado;\nel eje glucocorticoide no")
    ax=fig.add_subplot(gs[1,0]); placeholder(ax,"d  Estimador\n\nuna biopsia → fase y amplitud,\nvalidado contra la serie temporal\nde la misma persona\n\n(dibujo a mano)")
    # e amplitud inferida vs medida
    ax=fig.add_subplot(gs[1,1])
    M=pd.read_csv(f"{OUT}/new/GSE182117_inferred_free_vs_measured.tsv",sep="\t")
    for trt,c in (("CTL",GC),("HGI",MU)):
        q=M[M.trt==trt]; r=stats.spearmanr(q.amp_free,q.amp_mean_log2)
        ax.scatter(q.amp_mean_log2,q.amp_free,s=26,color=c,label=f"{trt}  ρ={r.statistic:.2f}")
    ax.set(xlabel="amplitud medida por cosinor (log2)",ylabel="amplitud inferida",
           title="e  Miotubos sincronizados:\nla amplitud inferida es la medida")
    ax.legend(fontsize=7,frameon=False)
    ax=fig.add_subplot(gs[1,2])
    G=pd.read_csv(f"{OUT}/new/GSE129843_amplitude_per_person.tsv",sep="\t")
    r=stats.spearmanr(G.amp_measured,G.amp_free)
    for cond,c in (("EXF",GC),("TRF",AD)):
        q=G[G.cond==cond]; ax.scatter(q.amp_measured,q.amp_free,s=30,color=c,label=cond)
    ax.set(xlabel="amplitud medida (cosinor, 6 puntos)",ylabel="amplitud inferida",
           title=f"f  Músculo in vivo: ρ = {r.statistic:.2f}\n(P = {r.pvalue:.1e}, 22 series)")
    ax.legend(fontsize=7,frameon=False)
    # g distribuciones por persona
    ax=fig.add_subplot(gs[1,3]); rng=np.random.default_rng(0)
    COH=[]
    X=pd.read_pickle(f"{OUT}/new/GSE25462_entrez.pkl"); md=pd.read_pickle(f"{OUT}/new/GSE25462_meta.pkl").set_index("gsm").reindex(X.columns)
    k=md.index[md.fh.isin(["DM","Family history negative"])]; COH.append(("GSE25462",X[list(k)],(md.fh.reindex(k)=="DM").values,MU))
    X=pd.read_pickle(f"{OUT}/new/GSE26637_entrez.pkl"); md=pd.read_pickle(f"{OUT}/new/GSE26637_meta.pkl").set_index("gsm").reindex(X.columns)
    k=md.index[md.stim=="fasting"]; COH.append(("GSE26637",X[list(k)],(md.status.reindex(k)=="resistant").values,AD))
    X=pd.read_pickle(f"{OUT}/new/GSE78721_log2.pkl"); md=pd.read_pickle(f"{OUT}/new/GSE78721_meta.pkl").set_index("gsm").reindex(X.columns)
    COH.append(("GSE78721",X,(md.dx!="Normal").values,AD))
    for i,(n,X,grp,c) in enumerate(COH):
        s=np.asarray(sc(center(X),SETS["Clock output"]["genes"]),float)
        for k2,m in enumerate((~grp,grp)):
            v=s[m]; xx=i*2.4+k2*.85
            ax.scatter(np.full(v.shape,xx)+rng.uniform(-.18,.18,len(v)),v,s=9,alpha=.6,color=c,marker="o" if k2 else "s")
            ax.hlines(np.median(v),xx-.28,xx+.28,color="k",lw=1.5)
        ax.text(i*2.4+.42,2.5,f"P={stats.mannwhitneyu(s[grp],s[~grp]).pvalue:.3f}",fontsize=7,ha="center")
    ax.set_xticks([i*2.4+.42 for i in range(3)]); ax.set_xticklabels([c[0] for c in COH],fontsize=7.5)
    ax.set(ylabel="score de salida del reloj (z)",ylim=(-2.7,2.9),
           title="g  Por persona: control (cuadrado)\nfrente a caso (círculo)")
    # h rho entre personas
    ax=fig.add_subplot(gs[2,0])
    R=pd.read_csv(f"{OUT}/new/meta_per_cohort.tsv",sep="\t")
    s2=R[(R.contraste=="NR3C1~reloj")&(R.cohorte.isin(["GSE25462","GSE182120","GSE13070"]))].sort_values("rho",ascending=False)
    lo=np.tanh(s2.z-1.96*s2.se); hi=np.tanh(s2.z+1.96*s2.se)
    ax.bar(range(len(s2)),s2.rho,color=[AD if "adiposo" in t else MU for t in s2.tejido],
           yerr=[s2.rho-lo,hi-s2.rho],capsize=4)
    ax.set_xticks(range(len(s2))); ax.set_xticklabels([f"{r.cohorte}\n{r.tejido}" for _,r in s2.iterrows()],fontsize=7)
    ax.axhline(0,color="k",lw=.8)
    ax.set(ylabel="ρ receptor ~ salida del reloj",title="h  La señal que falta está\nentre personas, no entre grupos")
    ax=fig.add_subplot(gs[2,2])
    for c,col in zip(["GSE25462_músculo","GSE78721_adiposo","GSE26637_adiposo"],[MU,AD,IS]):
        f=f"{OUT}/fig1_de_{c}.tsv"
        if os.path.exists(f):
            r=pd.read_csv(f,sep="\t",index_col=0)
            ax.hist(r.P.dropna(),bins=40,histtype="step",lw=1.4,color=col,density=True,label=c.replace("_"," "))
    ax.axhline(1,color="k",lw=.8,ls=":")
    ax.set(xlabel="P sin corregir, todos los genes",ylabel="densidad",
           title="i  La distribución de P es casi uniforme:\nno hay señal que la corrección esté ocultando")
    ax.legend(fontsize=7,frameon=False)
    ax=fig.add_subplot(gs[2,3]); placeholder(ax,"j  De comparar grupos\na correlacionar dentro del individuo\n\n(dibujo a mano)")
    tidy(fig); fig.savefig(f"{OUT}/figs/F1.png",dpi=160,bbox_inches="tight"); fig.savefig(f"{OUT}/figs/F1.pdf",bbox_inches="tight")

# ---------------------------------------------------------------- F2 (el eje)
def fig2():
    import subprocess
    subprocess.run([sys.executable,os.path.join(HERE,"fig3_eje.py")],check=True)
    for ext in ("png","pdf"):
        os.replace(f"{OUT}/figs/F3_eje.{ext}",f"{OUT}/figs/F2.{ext}")

# ---------------------------------------------------------------- F3 (cortisol sí, insulina no)
def fig3():
    fig=plt.figure(figsize=(16.5,9)); gs=fig.add_gridspec(2,4,hspace=.55,wspace=.38)
    ax=fig.add_subplot(gs[0,0]); placeholder(ax,"a  Las dos perturbaciones\n\nclamp: 0 h y 3 h, misma mañana\nhidrocortisona: 7:00 vs 7:00,\ncruzado con placebo\n\n(dibujo a mano)")
    C=pd.read_csv(f"{OUT}/new/time_calibration_vs_clamp.tsv",sep="\t",index_col=0)
    gl=["NR1D1","DBP","CIART","PER1","PER3","TEF"]; x=np.arange(len(gl)); w=.2
    ax=fig.add_subplot(gs[0,1])
    for k,(c,col) in enumerate(zip(C.columns,[GR,AD,"#999955",MU])):
        ax.bar(x+(k-1.5)*w,C.loc[gl,c],w,color=col,label=c.split(" ")[0].replace("GSE",""))
    ax.axhline(0,color="k",lw=.8); ax.set_xticks(x); ax.set_xticklabels(gl,rotation=45,fontsize=7.5)
    ax.set(ylabel="cambio (log2)",title="b  Tres horas de clamp =\ntres horas de mañana"); ax.legend(fontsize=6.5,frameon=False)
    # c clamp adiposo pareado (GSE26637)
    ax=fig.add_subplot(gs[0,2])
    D=pd.read_csv(f"{OUT}/new/GSE26637_clamp_adipose_per_person.tsv",sep="\t",index_col=0)
    gl2=["PDK4","TXNIP","DDIT4","PER1","NR1D2","CIART","PER3","HK2"]
    for i,g in enumerate(gl2):
        if g in D.index:
            v=D.loc[g].values
            ax.scatter(np.full(v.shape,i)+np.random.default_rng(i).uniform(-.16,.16,len(v)),v,s=14,color=AD,alpha=.7)
            ax.hlines(np.median(v),i-.3,i+.3,color="k",lw=1.6)
    ax.axhline(0,color="k",lw=.8); ax.set_xticks(range(len(gl2))); ax.set_xticklabels(gl2,rotation=45,fontsize=7.5)
    ax.set(ylabel="hiperinsulinemia − ayuno (log2)",title="c  Adiposo, clamp pareado (n=10):\nel brazo de salida cae en 10/10")
    # d PER1: sube con cortisol, baja con insulina
    ax=fig.add_subplot(gs[0,3])
    H=pd.read_csv(f"{OUT}/new/GSE148640_HC_response_stats.tsv",sep="\t",index_col=0)
    HC=pd.read_pickle(f"{OUT}/GSE148640_diff_genes.pkl")
    rng=np.random.default_rng(3)
    for i,(v,c,lab) in enumerate(((HC.loc["PER1"].values,IS,"hidrocortisona\nvs placebo"),
                                  (D.loc["PER1"].values,AD,"clamp\nvs ayuno"))):
        ax.scatter(np.full(v.shape,i)+rng.uniform(-.14,.14,len(v)),v,s=26,color=c,alpha=.8)
        ax.hlines(np.median(v),i-.28,i+.28,color="k",lw=2)
        ax.text(i,max(v)+.18,f"{int((v>0).sum())}/{len(v)} suben",ha="center",fontsize=7)
    ax.set_xticks([0,1]); ax.set_xticklabels(["hidrocortisona\nvs placebo","clamp\nvs ayuno"],fontsize=8)
    ax.axhline(0,color="k",lw=.8); ax.set(ylabel="PER1, cambio por persona (log2)",
       title="d  Dos perturbaciones opuestas\nsobre el mismo gen")
    # e hidrocortisona gen a gen
    ax=fig.add_subplot(gs[1,0])
    gl3=["PER1","PER3","BMAL1","CIART","NR1D1","DBP","TEF","ZBTB16","FKBP5","KLF15"]
    gl3=[g for g in gl3 if g in H.index]
    v=[H.loc[g,"median"] for g in gl3]; nup=[int(H.loc[g,"n_up"]) for g in gl3]
    ax.bar(range(len(gl3)),v,color=[IS if x>0 else MU for x in v])
    for i,(vv,n) in enumerate(zip(v,nup)): ax.text(i,vv+(0.08 if vv>0 else -0.16),f"{n}/10",ha="center",fontsize=6.5)
    ax.axhline(0,color="k",lw=.8); ax.set_xticks(range(len(gl3))); ax.set_xticklabels(gl3,rotation=45,fontsize=7.5)
    ax.set(ylabel="hidrocortisona − placebo (log2)",title="e  Misma hora del día:\nacción por genes concretos")
    # f dexametasona
    ax=fig.add_subplot(gs[1,1])
    L=pd.read_csv(f"{OUT}/new/GSE17992_dex_log2ratio_muscle.tsv",sep="\t",index_col=0)
    gl4=[g for g in ["KLF15","DBP","CRY2","CLOCK","TEF","HLF","PER3","BMAL1","ZBTB16"] if g in L.index]
    for i,g in enumerate(gl4):
        v=L.loc[g].dropna().values
        ax.scatter(np.full(v.shape,i)+np.random.default_rng(i).uniform(-.16,.16,len(v)),v,s=14,color=MU,alpha=.7)
        ax.hlines(np.median(v),i-.3,i+.3,color="k",lw=1.6)
    ax.axhline(0,color="k",lw=.8); ax.set_xticks(range(len(gl4))); ax.set_xticklabels(gl4,rotation=45,fontsize=7.5)
    ax.set(ylabel="dexametasona 4 d (log2)",title="f  Dexametasona en músculo (n=10)")
    # g Cushing
    ax=fig.add_subplot(gs[1,2])
    Xc_=pd.read_pickle(f"{OUT}/new/GSE66446_log2cpm.pkl"); Xc_.index=[{"ARNTL":"BMAL1"}.get(i,i) for i in Xc_.index]
    mdc=pd.read_pickle(f"{OUT}/new/GSE66446_meta.pkl").set_index("gsm").reindex(Xc_.columns)
    cu=mdc.c0.str.contains("Cushing").values; Xc2=center(Xc_)
    gl5=[g for g in ["PER3","DBP","PER2","TEF","HLF","NR1D1","CIART","BMAL1"] if g in Xc2.index]
    rng=np.random.default_rng(4)
    for i,g in enumerate(gl5):
        for k,(m,mk) in enumerate(((~cu,"s"),(cu,"o"))):
            v=Xc2.loc[g,m].values-Xc2.loc[g,~cu].median()
            xx=i*1.0+k*0.32-0.16
            ax.scatter(np.full(v.shape,xx)+rng.uniform(-.07,.07,len(v)),v,s=13,
                       color=IS if k else GR,marker=mk,alpha=.8)
            ax.hlines(np.median(v),xx-.13,xx+.13,color="k",lw=1.3)
    ax.axhline(0,color="k",lw=.8); ax.set_xticks(range(len(gl5))); ax.set_xticklabels(gl5,rotation=45,fontsize=7.5)
    ax.set(ylabel="log2 (centrado en el control)",title="g  Cushing adiposo: control (gris)\nfrente a Cushing (morado), 11 vs 5")
    # h firma de enfermedad frente a las tres
    ax=fig.add_subplot(gs[1,3])
    S=pd.read_csv(f"{OUT}/new/signatures_log2_common_genes.tsv",sep="\t",index_col=0)
    cols={"IR adiposo":"IR adipose (GSE13070, 35 vs 6)","Cushing":"Cushing adipose (GSE66446)",
          "hidrocortisona":"hydrocortisone adipose (GSE148640)","TZD":"TZD 3mo adipose (GSE13070)"}
    gl6=[g for g in ["PER3","DBP","HLF","TEF","CIART","NR1D1","NR1D2","PER2","BMAL1"] if g in S.index]
    T=S.loc[gl6,list(cols.values())]; T.columns=list(cols.keys())
    x=np.arange(len(gl6)); w=.2
    for k,(c,col) in enumerate(zip(T.columns,[MU,IS,GC,AD])):
        ax.bar(x+(k-1.5)*w,T[c],w,color=col,label=c)
    ax.axhline(0,color="k",lw=.8); ax.set_xticks(x); ax.set_xticklabels(gl6,rotation=45,fontsize=7.5)
    r1=stats.spearmanr(T["IR adiposo"],T["Cushing"]).statistic
    ax.set(ylabel="log2",title=f"h  El adiposo resistente lleva la firma\ndel cortisol (ρ={r1:.2f}); la TZD no la revierte")
    ax.legend(fontsize=6.5,frameon=False)
    tidy(fig); fig.savefig(f"{OUT}/figs/F3.png",dpi=160,bbox_inches="tight"); fig.savefig(f"{OUT}/figs/F3.pdf",bbox_inches="tight")


# ---------------------------------------------------------------- F4 (regulación y especificidad)
def fig4():
    fig=plt.figure(figsize=(16,12)); gs=fig.add_gridspec(3,3,hspace=.55,wspace=.62)
    T=pd.read_csv(f"{OUT}/new/AAGMEx_TF_activity_vs_clock.tsv",sep="\t")
    ax=fig.add_subplot(gs[0,0])
    top=T.nlargest(8,"rho"); bot=T.nsmallest(8,"rho"); P=pd.concat([bot,top])
    sz=20+120*(P.n_dianas-P.n_dianas.min())/(P.n_dianas.max()-P.n_dianas.min()+1e-9)
    ax.scatter(P.rho,range(len(P)),s=sz,c=P.rho,cmap="RdBu_r",vmin=-.45,vmax=.45,edgecolors="k",linewidths=.4)
    ax.set_yticks(range(len(P))); ax.set_yticklabels(P.TF,fontsize=7.5); ax.axvline(0,color="k",lw=.8)
    for n,yy in ((P.n_dianas.min(),-1.3),(P.n_dianas.max(),-2.3)):
        ax.scatter([-0.42],[yy],s=20+120*((n-P.n_dianas.min())/(P.n_dianas.max()-P.n_dianas.min()+1e-9)),
                   color="#cccccc",edgecolors="k",linewidths=.4,clip_on=False)
        ax.text(-0.38,yy,f"{int(n)} dianas",fontsize=6.5,va="center")
    ax.set(xlabel="ρ actividad del TF ~ salida del reloj",
           title=f"a  FOXO3 encabeza 134 factores;\nel eje NF-κB está al otro extremo")
    H=pd.read_csv(f"{OUT}/new/AAGMEx_NR3C1_enrichment_extra.tsv",sep="\t"); Hh=H[H.lib=="h_all"].copy()
    Hh["lab"]=Hh.set.str.replace("HALLMARK_","").str.replace("_"," ").str.capitalize()
    Hh=pd.concat([Hh.nlargest(3,"z"),Hh.nsmallest(11,"z")]).sort_values("z")
    ax=fig.add_subplot(gs[0,1])
    sz=20+130*(Hh.n-Hh.n.min())/(Hh.n.max()-Hh.n.min()+1e-9)
    sc_=ax.scatter(Hh.z,range(len(Hh)),s=sz,c=Hh.z,cmap="RdBu_r",vmin=-9,vmax=9,edgecolors="k",linewidths=.4)
    ax.set_yticks(range(len(Hh))); ax.set_yticklabels(Hh.lab,fontsize=7); ax.axvline(0,color="k",lw=.8)
    plt.colorbar(sc_,ax=ax,label="z",shrink=.75)
    ax.set(xlabel="z frente a conjuntos del mismo tamaño",
           title="b  El eje del receptor es antiinflamatorio\ny va contra OXPHOS y adipogénesis")
    E=pd.read_csv(f"{OUT}/new/AAGMEx_NR3C1_enrichment.tsv",sep="\t")
    hp=E[E.lib=="c5_hpo"].set_index("set")
    want={"HP_TRUNCAL_OBESITY":"obesidad troncal","HP_INCREASED_URINARY_CORTISOL_LEVEL":"cortisol urinario alto",
          "HP_INSULIN_RESISTANCE":"resistencia a la insulina","HP_HYPERGLYCEMIA":"hiperglucemia",
          "HP_TYPE_II_DIABETES_MELLITUS":"diabetes tipo 2"}
    lab=[];val=[];nn=[]
    for k,v in want.items():
        if k in hp.index: lab.append(v); val.append(float(hp.loc[k,"z"])); nn.append(float(hp.loc[k,"n"]))
    ax=fig.add_subplot(gs[0,2])
    sz=20+130*(np.array(nn)-min(nn))/(max(nn)-min(nn)+1e-9)
    sc2=ax.scatter(val,range(len(lab)),s=sz,c=val,cmap="RdBu_r",vmin=-4,vmax=4,edgecolors="k",linewidths=.4)
    ax.set_yticks(range(len(lab))); ax.set_yticklabels(lab,fontsize=8); ax.axvline(0,color="k",lw=.8)
    plt.colorbar(sc2,ax=ax,label="z",shrink=.75)
    ax.set(xlabel="z del término HPO en el eje del receptor",
           title="c  Fenotipo glucocorticoide a favor,\nfenotipo diabético en contra")
    # d TZD no revierte
    ax=fig.add_subplot(gs[1,0])
    S=pd.read_csv(f"{OUT}/new/signatures_log2_common_genes.tsv",sep="\t",index_col=0)
    gl=[g for g in ["PER3","DBP","HLF","TEF","CIART","NR1D1","NR1D2","BMAL1"] if g in S.index]
    ir=S.loc[gl,"IR adipose (GSE13070, 35 vs 6)"]; tz=S.loc[gl,"TZD 3mo adipose (GSE13070)"]
    ax.scatter(ir,tz,s=40,color=AD)
    for g in gl: ax.annotate(g,(ir[g],tz[g]),fontsize=6.5,xytext=(3,3),textcoords="offset points")
    ax.axhline(0,color="k",lw=.8); ax.axvline(0,color="k",lw=.8)
    r=stats.spearmanr(ir,tz)
    ax.set(xlabel="resistente − sensible (log2)",ylabel="TZD 3 meses (log2)",
           title=f"d  La TZD mejora la sensibilidad y no\nrevierte la firma (ρ = {r.statistic:.2f})")
    # e el receptor frente a la adiposidad
    ax=fig.add_subplot(gs[1,1])
    G=pd.read_csv(f"{OUT}/new/AAGMEx_GCreceptor_vs_adiposity.tsv",sep="\t")
    x=np.arange(len(G)); w=.38
    ax.bar(x-w/2,G.rho_vs_LEP,w,color=GR,label="ρ con LEP (adiposidad)")
    ax.bar(x+w/2,G.rho_vs_clock,w,color=AD,label="ρ con la salida del reloj")
    ax.axhline(0,color="k",lw=.8); ax.set_xticks(x); ax.set_xticklabels(G.gene,rotation=45,fontsize=7.5)
    ax.set(ylabel="ρ",title="e  Con la adiposidad sube la enzima\ny baja el receptor"); ax.legend(fontsize=7,frameon=False)
    # f esquema
    ax=fig.add_subplot(gs[1,2])
    R=pd.read_csv(f"{OUT}/new/AAGMEx_rank_vs_NR3C1.tsv",sep="\t",index_col=0)["rho"]
    ax.hist(R.dropna(),bins=60,color="#dddddd",edgecolor="none")
    for g,c in (("FOXO3",AD),("NFKB1",MU),("HSD11B1",IS)):
        if g in R.index: ax.axvline(R[g],color=c,lw=1.6); ax.text(R[g],ax.get_ylim()[1]*.88,g,rotation=90,fontsize=7,color=c,ha="right")
    ax.set_xlim(-0.6,0.6)
    ax.set(xlabel="ρ con NR3C1 (ajustado)",ylabel="genes",
           title="f  El eje del receptor en el\ntranscriptoma completo")
    ax=fig.add_subplot(gs[1,2] if False else gs[2,0]); placeholder(ax,"g  Arquitectura reguladora\n\nFOXO3 y el receptor sobre el brazo\nde salida; NF-κB en contra\n\n(dibujo a mano)")
    E2=pd.read_csv(f"{OUT}/new/AAGMEx_NR3C1_enrichment.tsv",sep="\t")
    for k,(lib,ttl) in enumerate((("c2_cp_reactome","h  Reactome"),("c5_go_bp","i  GO procesos biológicos"))):
        q=E2[E2.lib==lib].copy(); q=pd.concat([q.nlargest(4,"z"),q.nsmallest(10,"z")]).sort_values("z")
        q["lab"]=q.set.str.split("_").str[1:].str.join(" ").str.capitalize().str.slice(0,40)
        a=fig.add_subplot(gs[2,k+1])
        sz=20+130*(q.n-q.n.min())/(q.n.max()-q.n.min()+1e-9)
        s3=a.scatter(q.z,range(len(q)),s=sz,c=q.z,cmap="RdBu_r",vmin=-9,vmax=9,edgecolors="k",linewidths=.4)
        a.set_yticks(range(len(q))); a.set_yticklabels(q.lab,fontsize=6); a.axvline(0,color="k",lw=.8)
        plt.colorbar(s3,ax=a,label="z",shrink=.75)
        a.set(xlabel="z frente a conjuntos del mismo tamaño",title=ttl+" en el eje del receptor\n(tamaño: nº de genes)")
    tidy(fig); fig.savefig(f"{OUT}/figs/F4.png",dpi=160,bbox_inches="tight"); fig.savefig(f"{OUT}/figs/F4.pdf",bbox_inches="tight")

def fig5():
    import subprocess
    subprocess.run([sys.executable,os.path.join(HERE,"fig7_modelo.py")],check=True)
    for ext in ("png","pdf"):
        os.replace(f"{OUT}/figs/F7_modelo.{ext}",f"{OUT}/figs/F5.{ext}")

if __name__=="__main__":
    fig1(); print("F1 ok")
    fig2(); print("F2 ok")
    fig3(); print("F3 ok")
    fig4(); print("F4 ok")
    fig5(); print("F5 ok")
