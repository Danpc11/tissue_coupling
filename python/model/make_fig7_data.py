#!/usr/bin/env python3
"""Precalcula los datos de la figura 7 (lento: ~3 min). Guarda en results/new/fig7_data.npz"""
import os,sys,warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE); sys.path.insert(0,os.path.join(HERE,"../lib"))
import forced_clock as fc
OUT=os.environ.get("T2D_OUT","results")
NC=60; DAYS=18
def cohort(Km,seed,hour_sd=1.5,n=200,F=0.04,noise=0.6,Ksd=0.15,hour_mu=8.0):
    rng=np.random.default_rng(seed); K=np.clip(rng.normal(Km,Ksd,n),0.05,None)
    hour=rng.normal(hour_mu,hour_sd,n) if hour_sd<7 else rng.uniform(0,24,n)
    rows=[]
    for i,k in enumerate(K):
        phi,r=fc.simulate_population(F,k,n_cells=NC,days=DAYS,seed=int(rng.integers(1e6)))
        a,p,Rr=fc.ensemble(phi,r); rows.append(dict(subject=i,K=k,hour=hour[i],amp_ens=a,phase_ens=p,R=Rr))
    d=pd.DataFrame(rows); return d,fc.synthetic_expression(d,noise=noise,seed=int(rng.integers(1e6)))
D={}
# b) curva de forzamiento
fks=np.linspace(0,0.25,13); amps=[];Rs=[]
for fk in fks:
    phi,r=fc.simulate_population(fk,1.0,n_cells=150,days=DAYS,seed=1); a,_,Rr=fc.ensemble(phi,r); amps.append(a); Rs.append(Rr)
D["fks"]=fks; D["amps"]=amps; D["Rs"]=Rs
# c) firma
out=["PER3","DBP","HLF","TEF","NR1D1","NR1D2","PER2","CIART"]
dh,Eh=cohort(1.4,11); dl,El=cohort(0.5,12)
D["sig_sim"]=np.array([El.loc[g].mean()-Eh.loc[g].mean() for g in out]+[El.loc["BMAL1"].mean()-Eh.loc["BMAL1"].mean()])
D["sig_genes"]=np.array(out+["BMAL1"])
# d) ventana horaria
lv=[]
for sd in (1.5,4.0,7.5):
    a,Ea=cohort(1.4,21,hour_sd=sd); b,Eb=cohort(0.5,22,hour_sd=sd)
    lv.append([Eb.loc[out].mean().mean()-Ea.loc[out].mean().mean(), b.amp_ens.mean()-a.amp_ens.mean()])
D["window"]=np.array(lv)
# e) estimador
from scipy import stats
ev=[]
for noise in (0.3,0.6,1.0):
    d,E=cohort(1.0,31,hour_sd=7.5,n=200,noise=noise,Ksd=0.35)
    ang,amp=fc.run_estimator(E); ev.append([noise,stats.spearmanr(d.K,amp).statistic])
D["est"]=np.array(ev)
np.savez(f"{OUT}/new/fig7_data.npz",**D)
print("guardado", {k:(np.shape(v)) for k,v in D.items()})
