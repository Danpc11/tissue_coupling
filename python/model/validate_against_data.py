"""Comprueba si el modelo de oscilador forzado reproduce lo observado en biopsias unicas:
atenuacion del proxy NR3C1 frente al acoplamiento real y firma predicha con ventana horaria.
"""
import sys, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy import stats
sys.path.insert(0,"python/model"); sys.path.insert(0,"python/lib")
from forced_clock import synthetic_cohort, synthetic_expression, run_estimator
# El modelo da rho 0,45-0,73; lo observado es 0,09-0,22. ¿Qué falta? Dos cosas reales:
#  (a) NR3C1 medido no es K: es un proxy ruidoso (mRNA del receptor, no actividad)
#  (b) hay varianza de amplitud que no viene del acoplamiento (composición celular, calidad)
print("Atenuación: NR3C1 medido = K + ruido de medida. ¿Cuánta atenuación hace falta?\n")
d=synthetic_cohort(n_subjects=232,F=0.04,seed=3,n_cells=150)
E=synthetic_expression(d,noise=0.6,seed=4); ang,amp=run_estimator(E)
rng=np.random.default_rng(0)
print(f"{'rel. señal/ruido del proxy':>28s} {'ρ(proxy, amplitud)':>20s}")
for snr in (1.0,0.6,0.4,0.3,0.2,0.15,0.1):
    proxy=stats.zscore(d.K)*snr+rng.normal(0,1,len(d))*np.sqrt(1-snr**2)
    print(f"{snr:28.2f} {stats.spearmanr(proxy,amp).statistic:+20.3f}")
print("\n-> la ρ observada (0,09-0,22) corresponde a un proxy con correlación 0,15-0,35 con K real.")
print("   Es decir: el modelo es COMPATIBLE con lo observado si el mRNA de NR3C1 explica")
print("   entre el 2 % y el 12 % de la varianza de la ganancia real del acoplamiento.")
# predicción distintiva: la firma PER3/DBP/HLF abajo con BMAL1 arriba al bajar F*K
print("\n¿Reproduce el modelo la FIRMA (salida abajo, BMAL1 arriba)?")
hi=synthetic_cohort(n_subjects=400,K_mean=1.4,K_sd=0.15,F=0.04,seed=7,n_cells=150)
lo=synthetic_cohort(n_subjects=400,K_mean=0.5,K_sd=0.15,F=0.04,seed=8,n_cells=150)
Eh=synthetic_expression(hi,noise=0.6,seed=9); El=synthetic_expression(lo,noise=0.6,seed=10)
zh=((Eh.T-Eh.mean(1))/(Eh.std(1)+1e-9)).T; zl=((El.T-El.mean(1))/(El.std(1)+1e-9)).T
out=["PER3","DBP","HLF","TEF","NR1D1","NR1D2","PER1","PER2","CIART"]
print(f"   score de salida: K alto {zh.loc[out].mean().mean():+.3f} vs K bajo {zl.loc[out].mean().mean():+.3f}")
print(f"   BMAL1:           K alto {zh.loc['BMAL1'].mean():+.3f} vs K bajo {zl.loc['BMAL1'].mean():+.3f}")
print("   (con hora de biopsia uniforme, el modelo NO predice un cambio de nivel medio:")
print("    predice cambio de AMPLITUD. El desplazamiento de nivel que observamos requiere")
print("    además que las biopsias se tomen en una ventana horaria estrecha.)")
# con ventana horaria estrecha, como en las cohortes reales (biopsia matinal)
hi2=synthetic_cohort(n_subjects=400,K_mean=1.4,K_sd=0.15,F=0.04,seed=7,n_cells=150); hi2["hour"]=np.random.default_rng(1).normal(8,1.5,400)
lo2=synthetic_cohort(n_subjects=400,K_mean=0.5,K_sd=0.15,F=0.04,seed=8,n_cells=150); lo2["hour"]=np.random.default_rng(2).normal(8,1.5,400)
Eh=synthetic_expression(hi2,noise=0.6,seed=9); El=synthetic_expression(lo2,noise=0.6,seed=10)
zh=((Eh.T-Eh.mean(1))/(Eh.std(1)+1e-9)).T; zl=((El.T-El.mean(1))/(El.std(1)+1e-9)).T
print(f"\n   con biopsia matinal (8:00 +- 1,5 h), que es el diseño real:")
print(f"   score de salida: K alto {zh.loc[out].mean().mean():+.3f} vs K bajo {zl.loc[out].mean().mean():+.3f}")
print(f"   BMAL1:           K alto {zh.loc['BMAL1'].mean():+.3f} vs K bajo {zl.loc['BMAL1'].mean():+.3f}")
