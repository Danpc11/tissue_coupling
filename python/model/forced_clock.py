#!/usr/bin/env python3
"""Oscilador periférico forzado por el ritmo glucocorticoide: simulación y cohortes sintéticas.

Objetivo: comprobar si un modelo mínimo de entrainment reproduce lo que observamos en biopsias
únicas de seis cohortes humanas, sin ajustar nada a esos datos.

Estructura (Poincaré amplitude-phase forzado, la formulación estándar en cronobiología):

    dr/dt   = lambda * r * (A - r)                     relajación hacia el ciclo límite
    dphi/dt = omega + (F*K/r) * sin(theta_z - phi)     forzamiento en fase (término de acoplamiento)
    dr/dt  += F*K * cos(theta_z - phi)                 forzamiento en amplitud

con theta_z = 2*pi*t/24 la fase del Zeitgeber glucocorticoide. F es la amplitud RÍTMICA del cortisol
y K la ganancia del receptor: el forzamiento efectivo es su producto, que es la afirmación central
del modelo y la que une aplanamiento del ritmo (F baja), exposición sostenida (F baja aunque suba la
media) y caída del receptor (K baja).

Salidas del modelo por célula: fase y amplitud. Un tejido es una población de N células con períodos
intrínsecos dispersos (sigma_tau); la biopsia mide la MEDIA de la población en un instante, de modo
que la desincronización se lee como pérdida de amplitud aunque cada célula siga oscilando.

Brazos transcripcionales, con la misma lógica de §13:
  salida  (PER3, DBP, HLF, TEF, NR1D1/2): proporcional a la oscilación del ensamble
  directo (PER1):   término de GRE, sigue al nivel instantáneo de cortisol, no al ensamble
  BMAL1:            en antifase con el brazo de salida

Uso:  python python/model/forced_clock.py            (barrido + cohortes sintéticas)
"""
import numpy as np, pandas as pd

TAU_Z = 24.0

def simulate_population(F, K, n_cells=200, sigma_tau=1.2, lam=0.05, A=1.0,
                        days=30, dt=0.05, seed=0):
    """Integra N osciladores forzados y devuelve fase y radio de cada célula al final."""
    rng = np.random.default_rng(seed)
    tau = rng.normal(24.0, sigma_tau, n_cells)
    omega = 2 * np.pi / tau
    phi = rng.uniform(0, 2 * np.pi, n_cells)
    r = np.full(n_cells, A)
    FK = F * K
    for t in np.arange(0, days * 24, dt):
        thz = 2 * np.pi * t / TAU_Z
        d = thz - phi
        phi = phi + dt * (omega + FK * np.sin(d) / np.maximum(r, 1e-3))
        r = r + dt * (lam * r * (A - r) + FK * np.cos(d))
        r = np.clip(r, 1e-3, 10)
    return phi % (2 * np.pi), r

def ensemble(phi, r):
    """Lo que ve una biopsia: vector medio de la población."""
    z = np.mean(r * np.exp(1j * phi))
    return np.abs(z), np.angle(z) % (2 * np.pi), np.abs(np.mean(np.exp(1j * phi)))

def synthetic_cohort(n_subjects=200, K_mean=1.0, K_sd=0.35, F=0.06, noise=0.25,
                     n_cells=200, seed=0, **kw):
    """Cohorte de biopsias únicas: cada sujeto con su K, su hora de biopsia y ruido de medida."""
    rng = np.random.default_rng(seed)
    K = np.clip(rng.normal(K_mean, K_sd, n_subjects), 0.05, None)
    hour = rng.uniform(0, 24, n_subjects)          # hora de biopsia desconocida, como en las cohortes
    rows = []
    for i, k in enumerate(K):
        phi, r = simulate_population(F, k, n_cells=n_cells, seed=int(rng.integers(1e6)), **kw)
        amp, ph, R = ensemble(phi, r)
        th = 2 * np.pi * hour[i] / TAU_Z
        out = amp * np.cos(ph - th)                # brazo de salida, muestreado a la hora de biopsia
        rows.append(dict(subject=i, K=k, hour=hour[i], amp_ens=amp, phase_ens=ph, R=R,
                         output=out + rng.normal(0, noise),
                         bmal1=-out + rng.normal(0, noise),
                         per1_direct=np.cos(th) + rng.normal(0, noise)))
    return pd.DataFrame(rows)

if __name__ == "__main__":
    from scipy import stats
    print("== 1. BARRIDO DEL FORZAMIENTO EFECTIVO F*K ==")
    print(f"{'F*K':>7s} {'amplitud ensamble':>18s} {'R (sincronía)':>14s}")
    for fk in (0.0, 0.01, 0.02, 0.04, 0.06, 0.10, 0.20, 0.40):
        phi, r = simulate_population(fk, 1.0, seed=1)
        amp, ph, R = ensemble(phi, r)
        print(f"{fk:7.2f} {amp:18.3f} {R:14.3f}")
    print("\n== 2. COHORTE SINTÉTICA DE BIOPSIAS ÚNICAS (n=200, K variable entre personas) ==")
    for F in (0.03, 0.06, 0.12):
        d = synthetic_cohort(F=F, seed=2)
        r_amp = stats.spearmanr(d.K, d.amp_ens)
        q = pd.qcut(d.K, 5, labels=False)
        R1 = np.abs(np.mean(np.exp(1j * d.phase_ens[q == 0])))
        R5 = np.abs(np.mean(np.exp(1j * d.phase_ens[q == 4])))
        print(f"  F={F:.2f}: ρ(K, amplitud) = {r_amp.statistic:+.3f} (P {r_amp.pvalue:.1e}) | "
              f"R de fase quintil bajo {R1:.2f} -> alto {R5:.2f}")

# ---------------------------------------------------------------------------
def synthetic_expression(d, genes_phase=None, noise=0.6, seed=0):
    """Convierte la cohorte simulada en una matriz gen x sujeto, como la que daría un array.

    Cada gen del reloj es una sinusoide con su propia fase y su propio peso; la amplitud de todos
    ellos escala con la amplitud del ENSAMBLE de ese sujeto, y la fase se desplaza por la hora de
    biopsia, que es desconocida. El ruido domina la señal, como en los datos reales.
    """
    rng = np.random.default_rng(seed)
    if genes_phase is None:
        genes_phase = {"PER1": 8.0, "PER2": 10.0, "PER3": 10.5, "NR1D1": 4.0, "NR1D2": 5.0,
                       "DBP": 10.0, "TEF": 10.5, "HLF": 11.0, "CIART": 7.0,
                       "BMAL1": 22.0, "NPAS2": 22.5, "CLOCK": 21.0, "CRY1": 14.0, "CRY2": 12.0}
    M = {}
    for g, ph_h in genes_phase.items():
        a = 2 * np.pi * ph_h / TAU_Z
        th = 2 * np.pi * d.hour.values / TAU_Z
        M[g] = d.amp_ens.values * np.cos(th - d.phase_ens.values - a) + rng.normal(0, noise, len(d))
    return pd.DataFrame(M, index=d.subject).T

def run_estimator(E):
    """Aplica NUESTRO estimador (ordenamiento libre + amplitud por proyección) a la matriz."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../lib"))
    import phase as ph
    Z = E.copy(); Z = ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T
    ang = ph.unsup(Z.values)
    amp = ph.amplitude(E, ang, genes=list(E.index))
    return ang, amp
