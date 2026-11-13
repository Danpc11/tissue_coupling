#!/usr/bin/env python3
"""Barrido de la lengua de entrainment: forzamiento efectivo F*K frente a período intrínseco."""
import os,sys,warnings; warnings.filterwarnings("ignore")
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import forced_clock as fc
OUT=os.environ.get("T2D_OUT","results")
# Las mismas ecuaciones que el resto del modelo (fc.simulate_population), no una copia. Dos parametros
# difieren de los demas paneles y se dicen: la dispersion de periodos entre celulas es 1.5 h aqui y
# 1.2 h en el resto, y se usan 80 celulas y 20 dias.
FK=np.linspace(0,0.06,24); TAU=np.linspace(21.5,26.5,24); SIG=1.5
R=np.zeros((len(FK),len(TAU))); A=np.zeros_like(R)
for i,fk in enumerate(FK):
    for j,tau in enumerate(TAU):
        ph,rr=fc.simulate_population(1.0,float(fk),n_cells=80,sigma_tau=SIG,days=20,seed=100+j,tau_mean=float(tau))
        amp,_,Rs=fc.ensemble(ph,rr); R[i,j]=Rs; A[i,j]=amp
np.save(f"{OUT}/new/tongue_R.npy",R); np.save(f"{OUT}/new/tongue_A.npy",A)
np.save(f"{OUT}/new/tongue_axes.npy",np.array([FK,TAU]))
print("R rango",R.min().round(2),R.max().round(2))
