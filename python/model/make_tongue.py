#!/usr/bin/env python3
"""Barrido de la lengua de entrainment: forzamiento efectivo F*K frente a período intrínseco."""
import os,sys,warnings; warnings.filterwarnings("ignore")
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import forced_clock as fc
OUT=os.environ.get("T2D_OUT","results")
FK=np.linspace(0,0.06,24); TAU=np.linspace(21.5,26.5,24); SIG=1.5
R=np.zeros((len(FK),len(TAU))); A=np.zeros_like(R)
rng0=np.random.default_rng(7)
for i,fk in enumerate(FK):
    for j,tau in enumerate(TAU):
        n=80; rng=np.random.default_rng(100+j)
        taus=rng.normal(tau,SIG,n); om=2*np.pi/taus
        ph=rng.uniform(0,2*np.pi,n); rr=np.ones(n); dt=0.05
        for t in np.arange(0,20*24,dt):
            thz=2*np.pi*t/24.0; d=thz-ph
            ph=ph+dt*(om+fk*np.sin(d)/np.maximum(rr,1e-3))
            rr=np.clip(rr+dt*(0.05*rr*(1.0-rr)+fk*np.cos(d)),1e-3,10)
        amp,_,Rs=fc.ensemble(ph%(2*np.pi),rr); R[i,j]=Rs; A[i,j]=amp
np.save(f"{OUT}/new/tongue_R.npy",R); np.save(f"{OUT}/new/tongue_A.npy",A)
np.save(f"{OUT}/new/tongue_axes.npy",np.array([FK,TAU]))
print("R rango",R.min().round(2),R.max().round(2))
