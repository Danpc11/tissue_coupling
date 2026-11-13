#!/usr/bin/env python3
"""La prediccion del modelo, sin unidades y con fases medidas.

El modelo dice que, si baja el forzamiento, baja la amplitud del ensamble, y que una biopsia tomada
a una hora fija ve entonces cada gen desplazado hacia su media: los que a esa hora estan por encima
bajan y los que estan por debajo suben, en proporcion a su amplitud. La figura del modelo compara
esa firma, en unidades del modelo, con la observada en log2, y las dos escalas no son comparables:
la magnitud del modelo depende de una amplitud fijada en 1 y del contraste de ganancia elegido.

Aqui se comprueba lo que si es comparable, el PATRON entre genes, y sin usar las fases que el
modelo lleva escritas a mano sino las medidas en una cohorte independiente con horas registradas
(la plantilla de GSE108539):

    efecto predicho del gen g  =  - amplitud_g * cos( 2*pi*(hora de biopsia - fase_g) / 24 )

y se compara con la diferencia observada entre resistentes y sensibles: cuantos signos coinciden y
la correlacion de rangos entre genes. La hora de biopsia es la del modelo, las 08:00, fijada antes;
se informa tambien el resultado entre las 06:00 y las 11:00 para que se vea cuanto depende de ella.
El nulo asigna a los mismos valores observados la fase y la amplitud de genes de la plantilla
elegidos al azar.

Entradas: results/template_2harm.tsv, results/new/observed_signature_for_model.tsv
Salidas:  results/new/model_signature_from_phases.tsv, model_signature_from_phases_genes.tsv
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "../lib"))
from gene_sets import resolve
OUT = os.environ.get("T2D_OUT", "results")
HOUR = 8.0                      # la hora de biopsia del modelo (make_fig7_data.py)
NNULL, SEED = 5000, 0


def predicted(amp, phase_h, hour):
    return -amp * np.cos(2 * np.pi * (hour - phase_h) / 24)


if __name__ == "__main__":
    T = pd.read_csv(f"{OUT}/template_2harm.tsv", sep="\t", index_col=0)
    obs = pd.read_csv(f"{OUT}/new/observed_signature_for_model.tsv", sep="\t", index_col=0).iloc[:, 0]
    phase = (np.arctan2(T.b1, T.a1) % (2 * np.pi)) * 24 / (2 * np.pi); amp = np.hypot(T.a1, T.b1)
    rows = []
    for g in obs.index:
        r = resolve(g, T.index)
        if r is not None and np.isfinite(obs[g]):
            rows.append(dict(gene=g, observed_log2=obs[g], measured_phase_h=phase[r], measured_amplitude=amp[r]))
    D = pd.DataFrame(rows).set_index("gene")
    D["predicted_at_08h"] = predicted(D.measured_amplitude, D.measured_phase_h, HOUR)
    D["sign_agrees"] = np.sign(D.predicted_at_08h) == np.sign(D.observed_log2)
    D.to_csv(f"{OUT}/new/model_signature_from_phases_genes.tsv", sep="\t")
    print(f"{len(D)} genes con fase medida y diferencia observada: {', '.join(D.index)}")
    print(D.round(2).to_string())

    rng = np.random.default_rng(SEED); out = []
    P, A = phase.values, amp.values
    for hour in (6.0, 7.0, 8.0, 9.0, 10.0, 11.0):
        pred = predicted(D.measured_amplitude, D.measured_phase_h, hour)
        r = stats.spearmanr(pred, D.observed_log2).statistic
        null = np.empty(NNULL)
        for i in range(NNULL):
            k = rng.choice(len(T), len(D), replace=False)
            null[i] = stats.spearmanr(predicted(A[k], P[k], hour), D.observed_log2.values).statistic
        pval = (np.sum(null >= r) + 1) / (NNULL + 1)
        agree = int((np.sign(pred) == np.sign(D.observed_log2)).sum())
        out.append(dict(biopsy_hour=hour, prespecified=hour == HOUR, n_genes=len(D), signs_agree=agree, rho=r,
                        null_median=np.median(null), null_p95=np.percentile(null, 95), p_vs_random_genes=pval))
        print(f"  biopsia a las {hour:4.1f} h{' (la del modelo)' if hour == HOUR else '':16s}: signos {agree}/{len(D)}; rho {r:+.2f}; "
              f"genes al azar: mediana {np.median(null):+.2f}, percentil 95 {np.percentile(null, 95):+.2f}, P = {pval:.4f}")
    pd.DataFrame(out).to_csv(f"{OUT}/new/model_signature_from_phases.tsv", sep="\t", index=False)
