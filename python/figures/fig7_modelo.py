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
ax=fig.add_subplot(gs[0,0]); ax.text(-0.22,1.07,"A",transform=ax.transAxes,fontsize=11,fontweight="bold",va="top")
R=np.load(f"{OUT}/new/tongue_R.npy"); FK,TAU=np.load(f"{OUT}/new/tongue_axes.npy")
im=ax.pcolormesh(TAU,FK,R,cmap="viridis",shading="auto",vmin=0,vmax=1)
ax.contour(TAU,FK,R,levels=[0.6],colors="w",linewidths=1.6)
plt.colorbar(im,ax=ax,label="ensemble synchrony (R)")
ax.scatter([24],[0.045],marker="*",s=150,color="w",edgecolors="k",zorder=4)
ax.scatter([24],[0.012],marker="X",s=90,color="w",edgecolors="k",zorder=4)
ax.text(24.3,0.047,"coupling\npreserved",fontsize=7,color="w",va="center")
ax.text(24.3,0.014,"coupling\nreduced",fontsize=7,color="w",va="center")
ax.set(xlabel="intrinsic tissue period (h)",ylabel="effective forcing  F · K",
       )
ax=fig.add_subplot(gs[0,1]); ax.text(-0.22,1.07,"B",transform=ax.transAxes,fontsize=11,fontweight="bold",va="top")
ax.plot(D["fks"],D["amps"],"o-",color=BLUE)
ax.set(xlabel="effective forcing  F · K",ylabel="ensemble amplitude",
       )
ax2=ax.twinx(); ax2.plot(D["fks"],D["Rs"],"s--",color=ORANGE); ax2.set_ylabel("synchrony R",color=ORANGE)
ax2.spines[["top"]].set_visible(False)
ax=fig.add_subplot(gs[0,2]); ax.text(-0.22,1.07,"C",transform=ax.transAxes,fontsize=11,fontweight="bold",va="top")
sim=D["sig_sim"]; genes=list(D["sig_genes"])
OBS=pd.read_csv(f"{OUT}/new/observed_signature_for_model.tsv",sep="\t",index_col=0)["IR_minus_sensitive_log2"]
obs=np.array([OBS.get(g,np.nan) for g in genes])
x=np.arange(len(sim))
ax.bar(x-.2,sim,.4,color=PURPLE,label="model (low − high coupling)")
ax.bar(x+.2,obs,.4,color=GREEN,label="observed (resistant − sensitive)")
ax.axhline(0,color="k",lw=.8); ax.set_xticks(x); ax.set_xticklabels(genes,rotation=45,fontsize=7.5,ha="right")
ok=~np.isnan(obs); conc=int(np.sum(np.sign(sim[ok])==np.sign(obs[ok])))
ratio=np.abs(sim[ok]).mean()/np.abs(obs[ok]).mean()
ax.set(ylabel="change (log$_2$ / model units)",
       )
ax.legend(fontsize=7,frameon=False)
ax=fig.add_subplot(gs[1,0]); ax.text(-0.22,1.07,"D",transform=ax.transAxes,fontsize=11,fontweight="bold",va="top")
W=D["window"]; x=np.arange(3)
ax.bar(x-.2,W[:,0],.4,color=PURPLE,label="LEVEL difference")
ax.bar(x+.2,W[:,1],.4,color=GREY,label="AMPLITUDE difference")
ax.axhline(0,color="k",lw=.8); ax.set_xticks(x)
ax.set_xticklabels(["morning\n(08:00 ± 1.5 h)","wide\n(± 4 h)","uniform\n(24 h)"],fontsize=7.5)
ax.set(ylabel="low − high coupling",
       )
ax.legend(fontsize=7,frameon=False)
ax=fig.add_subplot(gs[1,1]); ax.text(-0.22,1.07,"E",transform=ax.transAxes,fontsize=11,fontweight="bold",va="top")
E=D["est"]
ax.plot(E[:,0],E[:,1],"o-",color=BLUE,label="true coupling (K)")
ax.axhspan(0.09,0.22,color=GREEN,alpha=.25)
ax.text(0.62,0.26,"range observed with NR3C1\n(0.09–0.22)",fontsize=7.5,color="#1d6b5f")
ax.set(xlabel="measurement noise",ylabel="ρ (coupling, inferred amplitude)",ylim=(0,0.85),
       )
ax.legend(fontsize=7,frameon=False)
ax=fig.add_subplot(gs[1,2]); ax.text(-0.22,1.07,"F",transform=ax.transAxes,fontsize=11,fontweight="bold",va="top"); ax.axis("off"); 
for xx,yy,s,c,fs,w in [(0,0.90,"effective forcing  =  F · K","k",10.5,"bold"),
     (0,0.81,"F = amplitude of the cortisol rhythm","k",8.5,"normal"),
     (0,0.74,"K = receptor availability","k",8.5,"normal"),
     (0,0.60,"three situations, one single effect:","k",9,"bold"),
     (0.02,0.51,"obesity, insulin resistance  →  F falls (flattened rhythm)","#d95f02",8,"normal"),
     (0.02,0.44,"Cushing, hydrocortisone  →  F falls (constant exposure)","#6a3d9a",8,"normal"),
     (0.02,0.37,"adiposity, chronic exposure  →  K falls (less receptor)","#2a9d8f",8,"normal"),
     (0,0.24,"↓ forcing → the ensemble desynchronises","k",9,"bold"),
     (0.02,0.16,"PER3, DBP, HLF, TEF, NR1D1/2 ↓     BMAL1 ↑","k",8,"normal"),
     (0,0.04,"acute insulin does not enter the product","#777777",8.5,"italic")]:
    ax.text(xx,yy,s,color=c,fontsize=fs,fontweight="bold" if w=="bold" else "normal",
            style="italic" if w=="italic" else "normal",transform=ax.transAxes)
for a in fig.axes:
    a.spines[["top","right"]].set_visible(False); a.title.set_fontsize(8.5); a.title.set_ha("left"); a.title.set_position([0,1])
ax2.spines[["right"]].set_visible(True)
fig.savefig(f"{OUT}/figs/F7_modelo.png",dpi=400,bbox_inches="tight"); fig.savefig(f"{OUT}/figs/F7_modelo.pdf",bbox_inches="tight")
print("ok")
