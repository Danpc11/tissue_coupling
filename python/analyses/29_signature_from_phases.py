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

    # Dos nulos, y en los dos se sacan primero los conjuntos de genes y luego se evaluan a todas las horas,
    # para poder preguntar tambien por el maximo entre horas:
    #   any template gene        genes de la plantilla tomados al azar, con su fase y su amplitud medidas.
    #   matched on amplitude     cada gen se sustituye por uno de la plantilla de su mismo decil de amplitud:
    #                            conserva cuanto oscilan los genes usados y deja libre solo la fase.
    # Ninguno conserva la coexpresion entre los ocho genes, que la plantilla no contiene: siete son del brazo
    # de salida y comparten fase, de modo que la coincidencia equivale a dos hechos (ese brazo esta cerca de
    # su maximo a la hora de la biopsia y mas bajo en los resistentes; BMAL1, en antifase, va al reves).
    rng = np.random.default_rng(SEED); out = []
    P, A = phase.values, amp.values
    HOURS = (6.0, 7.0, 8.0, 9.0, 10.0, 11.0)
    decile = np.minimum((stats.rankdata(A) - 1) * 10 // len(A), 9).astype(int)
    pools = [np.flatnonzero(decile == decile[list(T.index).index(resolve(g, T.index))]) for g in D.index]
    draws = {"any template gene": np.array([rng.choice(len(T), len(D), replace=False) for _ in range(NNULL)]),
             "matched on amplitude": np.array([[rng.choice(p) for p in pools] for _ in range(NNULL)])}
    obs_v = D.observed_log2.values
    rank_obs = stats.rankdata(obs_v); rank_obs = (rank_obs - rank_obs.mean()) / rank_obs.std()
    def rho_rows(M):                                         # Spearman de cada fila de M con lo observado
        R = stats.rankdata(M, axis=1); R = (R - R.mean(1, keepdims=True)) / R.std(1, keepdims=True)
        return (R * rank_obs).mean(1)
    NULL = {k: np.array([rho_rows(predicted(A[d], P[d], h)) for h in HOURS]) for k, d in draws.items()}          # horas x nulos
    AGREE = {k: np.array([(np.sign(predicted(A[d], P[d], h)) == np.sign(obs_v)).sum(1) for h in HOURS]) for k, d in draws.items()}
    OBS = np.array([stats.spearmanr(predicted(D.measured_amplitude, D.measured_phase_h, h), obs_v).statistic for h in HOURS])
    for j, hour in enumerate(HOURS):
        pred = predicted(D.measured_amplitude, D.measured_phase_h, hour); r = OBS[j]
        agree = int((np.sign(pred) == np.sign(D.observed_log2)).sum())
        row = dict(biopsy_hour=hour, prespecified=hour == HOUR, n_genes=len(D), signs_agree=agree, rho=r)
        for k in NULL:
            tag = k.replace(" ", "_")
            row[f"p_rho_vs_{tag}"] = (np.sum(NULL[k][j] >= r) + 1) / (NNULL + 1)
            row[f"p_signs_vs_{tag}"] = (np.sum(AGREE[k][j] >= agree) + 1) / (NNULL + 1)
            row[f"signs_null_median_{tag}"] = float(np.median(AGREE[k][j]))
        row["null_median"] = float(np.median(NULL["any template gene"][j])); row["null_p95"] = float(np.percentile(NULL["any template gene"][j], 95))
        out.append(row)
        print(f"  biopsia a las {hour:4.1f} h{' (la del modelo)' if hour == HOUR else '':16s}: signos {agree}/{len(D)}; rho {r:+.2f}; "
              f"P frente a genes al azar {row['p_rho_vs_any_template_gene']:.4f}, de la misma amplitud {row['p_rho_vs_matched_on_amplitude']:.4f}"
              f" | signos: P {row['p_signs_vs_any_template_gene']:.4f} y {row['p_signs_vs_matched_on_amplitude']:.4f}")
    # La hora mejor, corregida por haberla buscado: el maximo entre las seis horas se compara con el maximo
    # entre las seis horas de cada conjunto nulo. La hora preespecificada es las 08:00; esto es lo que vale
    # la mejor de las exploradas.
    jbest = int(np.argmax(OBS)); best = dict(biopsy_hour=HOURS[jbest], prespecified=False, n_genes=len(D), rho=OBS[jbest],
                                             signs_agree=int((np.sign(predicted(D.measured_amplitude, D.measured_phase_h, HOURS[jbest])) == np.sign(obs_v)).sum()))
    for k in NULL:
        best[f"p_rho_vs_{k.replace(' ', '_')}"] = (np.sum(NULL[k].max(0) >= OBS[jbest]) + 1) / (NNULL + 1)
    best["selected_as_the_best_of_six_hours"] = True
    print(f"  la mejor de las seis horas ({HOURS[jbest]:.0f} h, rho {OBS[jbest]:+.2f}), corregida por la busqueda: "
          f"P = {best['p_rho_vs_any_template_gene']:.4f} frente a genes al azar, {best['p_rho_vs_matched_on_amplitude']:.4f} de la misma amplitud")
    R = pd.DataFrame(out); R["selected_as_the_best_of_six_hours"] = False
    pd.concat([R, pd.DataFrame([best])], ignore_index=True).to_csv(f"{OUT}/new/model_signature_from_phases.tsv", sep="\t", index=False)
