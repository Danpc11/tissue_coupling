#!/usr/bin/env python3
"""Figura 7. Un oscilador forzado por el ritmo glucocorticoide reproduce lo observado.
Requiere python/model/make_fig7_data.py (precalcula las simulaciones) y /tmp/tongue (barrido)."""
import os,sys,warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
OUT=os.environ.get("T2D_OUT","results")
BLUE,ORANGE,GREEN,PURPLE,GREY="#1f77b4","#d95f02","#2a9d8f","#6a3d9a","#999999"
D=np.load(f"{OUT}/new/fig7_data.npz",allow_pickle=True)
fig=plt.figure(figsize=(16,8.8)); gs=fig.add_gridspec(2,3,hspace=.5,wspace=.35)
ax=fig.add_subplot(gs[0,0])
R=np.load(f"{OUT}/new/tongue_R.npy"); FK,TAU=np.load(f"{OUT}/new/tongue_axes.npy")
im=ax.pcolormesh(TAU,FK,R,cmap="viridis",shading="auto",vmin=0,vmax=1)
ax.contour(TAU,FK,R,levels=[0.6],colors="w",linewidths=1.6)
plt.colorbar(im,ax=ax,label="sincronía del ensamble (R)")
ax.scatter([24],[0.045],marker="*",s=150,color="w",edgecolors="k",zorder=4)
ax.scatter([24],[0.012],marker="X",s=90,color="w",edgecolors="k",zorder=4)
ax.text(24.3,0.047,"acoplamiento\nconservado",fontsize=7,color="w",va="center")
ax.text(24.3,0.014,"acoplamiento\nreducido",fontsize=7,color="w",va="center")
ax.set(xlabel="período intrínseco del tejido (h)",ylabel="forzamiento efectivo  F · K",
       title="a  Lengua de entrainment: el forzamiento es el\nproducto de ritmo de cortisol y receptor")
ax=fig.add_subplot(gs[0,1])
ax.plot(D["fks"],D["amps"],"o-",color=BLUE)
ax.set(xlabel="forzamiento efectivo  F · K",ylabel="amplitud del ensamble",
       title="b  Sin Zeitgeber las células siguen oscilando\npero el ensamble se apaga")
ax2=ax.twinx(); ax2.plot(D["fks"],D["Rs"],"s--",color=ORANGE); ax2.set_ylabel("sincronía R",color=ORANGE)
ax2.spines[["top"]].set_visible(False)
ax=fig.add_subplot(gs[0,2])
sim=D["sig_sim"]; genes=list(D["sig_genes"])
OBS=pd.read_csv(f"{OUT}/new/observed_signature_for_model.tsv",sep="\t",index_col=0)["IR_minus_sensitive_log2"]
obs=np.array([OBS.get(g,np.nan) for g in genes])
x=np.arange(len(sim))
ax.bar(x-.2,sim,.4,color=PURPLE,label="modelo (acoplamiento bajo − alto)")
ax.bar(x+.2,obs,.4,color=GREEN,label="observado (resistente − sensible)")
ax.axhline(0,color="k",lw=.8); ax.set_xticks(x); ax.set_xticklabels(genes,rotation=45,fontsize=7.5,ha="right")
ok=~np.isnan(obs); conc=int(np.sum(np.sign(sim[ok])==np.sign(obs[ok])))
ratio=np.abs(sim[ok]).mean()/np.abs(obs[ok]).mean()
ax.set(ylabel="cambio (log2 / unidades del modelo)",
       title=f"c  Sin ajustar ningún parámetro: {conc}/{int(ok.sum())} signos coinciden\ny la magnitud media está dentro del {abs(1-ratio)*100:.0f} %")
ax.legend(fontsize=7,frameon=False)
ax=fig.add_subplot(gs[1,0])
W=D["window"]; x=np.arange(3)
ax.bar(x-.2,W[:,0],.4,color=PURPLE,label="diferencia de NIVEL")
ax.bar(x+.2,W[:,1],.4,color=GREY,label="diferencia de AMPLITUD")
ax.axhline(0,color="k",lw=.8); ax.set_xticks(x)
ax.set_xticklabels(["matinal\n(8:00 ± 1,5 h)","amplia\n(± 4 h)","uniforme\n(24 h)"],fontsize=7.5)
ax.set(ylabel="bajo − alto acoplamiento",
       title="d  Predicción falsable: el cambio de nivel\nsolo existe con biopsia a hora fija")
ax.legend(fontsize=7,frameon=False)
ax=fig.add_subplot(gs[1,1])
E=D["est"]
ax.plot(E[:,0],E[:,1],"o-",color=BLUE,label="acoplamiento verdadero (K)")
ax.axhspan(0.09,0.22,color=GREEN,alpha=.25)
ax.text(0.62,0.26,"rango observado con NR3C1\n(0,09–0,22)",fontsize=7.5,color="#1d6b5f")
ax.set(xlabel="ruido de medida",ylabel="ρ (acoplamiento, amplitud inferida)",ylim=(0,0.85),
       title="e  Nuestro estimador sobre datos del modelo:\nel mRNA del receptor es un proxy atenuado")
ax.legend(fontsize=7,frameon=False)
ax=fig.add_subplot(gs[1,2]); ax.axis("off"); ax.set_title("f  El mecanismo, en una línea",fontsize=9.5,loc="left")
for xx,yy,s,c,fs,w in [(0,0.90,"forzamiento efectivo  =  F · K","k",10.5,"bold"),
     (0,0.81,"F = amplitud del ritmo de cortisol","k",8.5,"normal"),
     (0,0.74,"K = disponibilidad de receptor","k",8.5,"normal"),
     (0,0.60,"tres situaciones, un solo efecto:","k",9,"bold"),
     (0.02,0.51,"obesidad y resistencia  →  F baja (ritmo aplanado)","#d95f02",8,"normal"),
     (0.02,0.44,"Cushing, hidrocortisona  →  F baja (exposición constante)","#6a3d9a",8,"normal"),
     (0.02,0.37,"adiposidad, exposición   →  K baja (menos receptor)","#2a9d8f",8,"normal"),
     (0,0.24,"↓ forzamiento → el ensamble se desincroniza","k",9,"bold"),
     (0.02,0.16,"PER3, DBP, HLF, TEF, NR1D1/2 ↓     BMAL1 ↑","k",8,"normal"),
     (0,0.04,"la insulina aguda no entra en el producto","#777777",8.5,"italic")]:
    ax.text(xx,yy,s,color=c,fontsize=fs,fontweight="bold" if w=="bold" else "normal",
            style="italic" if w=="italic" else "normal",transform=ax.transAxes)
for a in fig.axes:
    a.spines[["top","right"]].set_visible(False); a.title.set_fontsize(9.5)
ax2.spines[["right"]].set_visible(True)
fig.savefig(f"{OUT}/figs/F7_modelo.png",dpi=160,bbox_inches="tight"); fig.savefig(f"{OUT}/figs/F7_modelo.pdf",bbox_inches="tight")
print("ok")
