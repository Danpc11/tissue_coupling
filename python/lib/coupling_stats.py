"""Distribucion del acoplamiento individual y contrastes sobre ella.

La medida de `within_person_coupling` da una cifra por persona, de modo que a partir de aqui se
puede aplicar estadistica ordinaria: describir su distribucion, descomponer su varianza, comparar
grupos y calcular potencia. Tres cautelas gobiernan todo lo que sigue.

1. EL CONTRASTE NO ES CONTRA CERO. Dos personas sin relacion comparten la misma plantilla, y eso
   solo ya produce un fondo. El nulo es emparejar el tejido de una persona con el de otra, y es
   el que devuelve `coupling` en la columna `null`.

2. LA VARIANZA OBSERVADA ESTA INFLADA por el error de estimacion. La fiabilidad por mitades de
   genes dice que fraccion es real; la desviacion verdadera es la observada por la raiz de la
   fiabilidad. Sin esa correccion se sobreestima la variabilidad entre personas.

3. LA MEDIDA NO ES INDEPENDIENTE DE LA DIFERENCIA DE FASE (rho = -0.57 en GTEx). Al asociarla con
   un fenotipo conviene comprobar si la asociacion sobrevive al ajuste por |dphase|, porque si no,
   lo que se esta midiendo es la diferencia de fase bajo otro nombre.
"""
import numpy as np
import pandas as pd
from scipy import stats


def describe(d, reliability=None):
    """Resumen de la distribucion, con la desviacion corregida por fiabilidad si se aporta."""
    r = d.coupling.values
    nl = d.null.values
    out = dict(n=len(r), mean=r.mean(), median=float(np.median(r)), sd=r.std(ddof=1),
               q05=np.percentile(r, 5), q25=np.percentile(r, 25),
               q75=np.percentile(r, 75), q95=np.percentile(r, 95),
               skew=float(stats.skew(r)), kurtosis=float(stats.kurtosis(r)),
               null_median=float(np.median(nl)), null_sd=nl.std(ddof=1),
               above_null_p95=float(np.mean(r > np.percentile(nl, 95))))
    if reliability is not None:
        out["reliability"] = reliability
        out["sd_true"] = r.std(ddof=1) * np.sqrt(max(reliability, 0))
        out["var_fraction_real"] = reliability
    return out


def n_components(d, kmax=3, seed=0):
    """Numero de componentes por BIC en una mezcla de gaussianas.

    Una mezcla con una componente negativa suele indicar lecturas de fase a medio ciclo, no
    subpoblaciones: comprobar `dphase` antes de interpretarla como biologia.
    """
    from sklearn.mixture import GaussianMixture
    X = d.coupling.values.reshape(-1, 1)
    bic = [GaussianMixture(k, random_state=seed).fit(X).bic(X) for k in range(1, kmax + 1)]
    k = int(np.argmin(bic)) + 1
    g = GaussianMixture(k, random_state=seed).fit(X)
    o = np.argsort(g.means_.ravel())
    comp = [dict(mean=float(g.means_.ravel()[i]), sd=float(np.sqrt(g.covariances_.ravel()[i])),
                 weight=float(g.weights_[i])) for i in o]
    return dict(k=k, bic=bic, components=comp)


def compare_groups(d, groups, adjust_dphase=True):
    """Compara el acoplamiento entre grupos, con y sin ajuste por la diferencia de fase.

    Si la diferencia desaparece al ajustar, lo que distingue a los grupos es la fase y no el
    acoplamiento.
    """
    g = pd.Series(groups, index=d.donor.values) if not isinstance(groups, pd.Series) else groups
    d = d.assign(grp=[g.get(x, np.nan) for x in d.donor]).dropna(subset=["grp"])
    lv = sorted(d.grp.unique())
    if len(lv) != 2:
        raise ValueError(f"se esperan dos grupos, hay {len(lv)}")
    a = d[d.grp == lv[0]].coupling.values
    b = d[d.grp == lv[1]].coupling.values
    out = dict(group_a=lv[0], n_a=len(a), median_a=float(np.median(a)),
               group_b=lv[1], n_b=len(b), median_b=float(np.median(b)),
               difference=float(np.median(a) - np.median(b)),
               p_mannwhitney=float(stats.mannwhitneyu(a, b).pvalue),
               cohens_d=float((a.mean() - b.mean()) /
                              np.sqrt(((len(a) - 1) * a.var(ddof=1) +
                                       (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))))
    if adjust_dphase and "dphase" in d.columns:
        X = np.column_stack([np.ones(len(d)), (d.grp == lv[0]).astype(float),
                             np.abs(d.dphase.values)])
        beta, *_ = np.linalg.lstsq(X, d.coupling.values, rcond=None)
        res = d.coupling.values - X @ beta
        se = np.sqrt((res ** 2).sum() / (len(d) - 3) * np.linalg.pinv(X.T @ X)[1, 1])
        out["beta_group_adjusted"] = float(beta[1])
        out["p_adjusted_for_dphase"] = float(2 * stats.t.sf(abs(beta[1] / se), len(d) - 3))
    return out


def power(sd_true, difference, n_per_group, alpha=0.05, nsim=2000, seed=0):
    """Potencia para detectar una diferencia de medias dada la dispersion REAL entre personas.

    sd_true debe ser la desviacion corregida por fiabilidad, no la observada: usar la observada
    sobreestima la dispersion y por tanto subestima la potencia.
    """
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(nsim):
        a = rng.normal(0, sd_true, n_per_group)
        b = rng.normal(difference, sd_true, n_per_group)
        if stats.mannwhitneyu(a, b).pvalue < alpha:
            hits += 1
    return hits / nsim


def n_required(sd_true, difference, target=0.8, alpha=0.05, grid=(10, 20, 40, 80, 150, 300, 600)):
    """Tamano por grupo necesario para alcanzar la potencia objetivo."""
    for n in grid:
        if power(sd_true, difference, n, alpha, nsim=800) >= target:
            return n
    return None


# ---------------------------------------------------------------------------------------------
# CONTRASTES DE DISTRIBUCION
#
# El acoplamiento varia entre personas y esa variacion es el objeto de estudio, no ruido: depende
# de las senales de arrastre -horario de comidas, sueno, actividad, luz- que difieren entre
# individuos. Esperar un valor unico de la especie seria lo extrano.
#
# Eso cambia como se contrasta una hipotesis de enfermedad. Comparar medianas asume que la
# condicion desplaza TODA la distribucion. Son igual de plausibles, y distinguibles, otras dos
# formas: que enriquezca la cola baja sin mover el grueso, o que ensanche la distribucion. Una
# prueba de medianas no ve ninguna de las dos. Las funciones siguientes contrastan la
# distribucion completa.
# ---------------------------------------------------------------------------------------------


def distribution_tests(x, y, labels=("a", "b"), nperm=5000, seed=0):
    """Bateria de contrastes sobre dos distribuciones de acoplamiento.

    Devuelve, ademas del desplazamiento central, pruebas sensibles a la FORMA: Kolmogorov-Smirnov
    sobre toda la distribucion, Anderson-Darling (mas sensible en las colas, que es donde se
    esperaria el enriquecimiento de un subgrupo mal acoplado), y una comparacion de dispersiones.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    out = dict(label_a=labels[0], n_a=len(x), label_b=labels[1], n_b=len(y),
               median_a=float(np.median(x)), median_b=float(np.median(y)),
               sd_a=float(x.std(ddof=1)), sd_b=float(y.std(ddof=1)))
    out["p_median_shift"] = float(stats.mannwhitneyu(x, y).pvalue)
    out["p_ks"] = float(stats.ks_2samp(x, y).pvalue)
    try:
        out["p_anderson"] = float(stats.anderson_ksamp([x, y]).significance_level)
    except Exception:
        out["p_anderson"] = np.nan
    # dispersion: Levene sobre desviaciones absolutas a la mediana, robusto a no normalidad
    out["p_dispersion"] = float(stats.levene(x, y, center="median").pvalue)
    return out


def tail_enrichment(x, y, threshold=None, labels=("a", "b"), nperm=5000, seed=0):
    """Proporcion de individuos por debajo de un umbral en cada grupo.

    Contrasta la hipotesis de que la condicion ENRIQUECE la cola de acoplamiento bajo sin
    desplazar el grueso de la distribucion. Si no se da umbral se usa el cuartil inferior del
    grupo de referencia (el primero), de modo que el umbral no se elige mirando el resultado.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    thr = float(np.percentile(x, 25)) if threshold is None else float(threshold)
    px, py = float(np.mean(x <= thr)), float(np.mean(y <= thr))
    tab = [[int((x <= thr).sum()), int((x > thr).sum())],
           [int((y <= thr).sum()), int((y > thr).sum())]]
    p = float(stats.fisher_exact(tab).pvalue)
    return dict(threshold=thr, label_a=labels[0], frac_a=px, label_b=labels[1], frac_b=py,
                odds_ratio=float(stats.fisher_exact(tab).statistic), p_fisher=p)


def mixture_shift(x, y, labels=("a", "b"), seed=0):
    """Ajusta una mezcla de dos componentes al conjunto y compara el PESO del componente bajo.

    Formaliza la hipotesis de subgrupo: las dos distribuciones comparten los mismos componentes y
    lo que cambia entre condiciones es la proporcion de personas en el de acoplamiento bajo. Si
    los componentes resultan muy solapados, la descomposicion no es informativa y conviene
    quedarse con tail_enrichment.
    """
    from sklearn.mixture import GaussianMixture
    x = np.asarray(x, float); y = np.asarray(y, float)
    g = GaussianMixture(2, random_state=seed).fit(np.concatenate([x, y]).reshape(-1, 1))
    lo = int(np.argmin(g.means_.ravel()))
    wx = float(g.predict_proba(x.reshape(-1, 1))[:, lo].mean())
    wy = float(g.predict_proba(y.reshape(-1, 1))[:, lo].mean())
    sep = float(abs(np.diff(g.means_.ravel())[0]) / np.sqrt(g.covariances_.ravel().mean()))
    return dict(mean_low=float(g.means_.ravel()[lo]),
                mean_high=float(g.means_.ravel()[1 - lo]),
                separation_sd=sep,
                label_a=labels[0], weight_low_a=wx,
                label_b=labels[1], weight_low_b=wy,
                difference=wy - wx,
                interpretable=bool(sep > 2))


def compare_distributions(d, groups, threshold=None, seed=0):
    """Contraste completo entre dos grupos: desplazamiento, forma, cola y mezcla."""
    g = pd.Series(groups, index=d.donor.values) if not isinstance(groups, pd.Series) else groups
    dd = d.assign(grp=[g.get(x, np.nan) for x in d.donor]).dropna(subset=["grp"])
    lv = sorted(dd.grp.unique())
    if len(lv) != 2:
        raise ValueError(f"se esperan dos grupos, hay {len(lv)}")
    x = dd[dd.grp == lv[0]].coupling.values
    y = dd[dd.grp == lv[1]].coupling.values
    return dict(tests=distribution_tests(x, y, lv, seed=seed),
                tail=tail_enrichment(x, y, threshold, lv, seed=seed),
                mixture=mixture_shift(x, y, lv, seed=seed))


# ---------------------------------------------------------------------------------------------
# DECONVOLUCION DEL ERROR DE MEDIDA
#
# La distribucion observada es la verdadera convolucionada con el error de estimacion. Eso NO
# sesga la media -la convolucion con un error centrado deja la media intacta- pero sí ensancha la
# distribucion y borra su estructura: una subpoblacion mal acoplada puede quedar emborronada hasta
# desaparecer. Cuanto menor la fiabilidad, mas severo el efecto (0.72 en GTEx, 0.51 en GSE13070).
#
# La deconvolucion no anade informacion y es un problema mal condicionado: amplifica el ruido de
# alta frecuencia. Por eso se regulariza y se reporta siempre junto a la observada, nunca en su
# lugar.
# ---------------------------------------------------------------------------------------------


def error_sd(d, reliability):
    """Desviacion del error de medida, a partir de la fiabilidad y la varianza observada."""
    v = d.coupling.var(ddof=1)
    return float(np.sqrt(max(v * (1 - reliability), 0)))


def deconvolve(x, sd_err, grid=None, n_grid=200, iters=400, damping=0.02, seed=0):
    """Recupera la distribucion verdadera por Richardson-Lucy sobre la rejilla.

    x       : valores observados
    sd_err  : desviacion del error de medida (de error_sd)
    damping : suavizado por iteracion; sin el, la solucion se fragmenta en picos espurios

    Devuelve la rejilla y la densidad estimada. Comprobar siempre que la media se conserva: si se
    desplaza, la regularizacion esta distorsionando y no solo afinando.
    """
    x = np.asarray(x, float)
    if grid is None:
        pad = 3 * sd_err + 0.05
        grid = np.linspace(x.min() - pad, x.max() + pad, n_grid)
    h = grid[1] - grid[0]
    obs, _ = np.histogram(x, bins=np.r_[grid - h / 2, grid[-1] + h / 2])
    obs = obs.astype(float) / max(obs.sum(), 1)
    if sd_err <= 0:
        return grid, obs
    K = np.exp(-0.5 * ((grid[:, None] - grid[None, :]) / sd_err) ** 2)
    K /= K.sum(0, keepdims=True)
    f = np.ones_like(grid) / len(grid)
    for _ in range(iters):
        c = K @ f
        f = f * (K.T @ (obs / np.maximum(c, 1e-12)))
        if damping:                      # suavizado gaussiano ligero entre iteraciones
            w = np.exp(-0.5 * ((grid[:, None] - grid[None, :]) / (damping + h)) ** 2)
            w /= w.sum(1, keepdims=True)
            f = w @ f
        f = np.maximum(f, 0); f /= max(f.sum(), 1e-12)
    return grid, f


def deconvolved_summary(d, reliability, **kw):
    """Resumen de la distribucion antes y despues de deconvolucionar.

    La comparacion util es la DE y la forma; la media debe mantenerse, y si no lo hace es senal de
    que la regularizacion esta distorsionando.
    """
    x = d.coupling.values
    se = error_sd(d, reliability)
    g, f = deconvolve(x, se, **kw)
    mu = float((g * f).sum())
    sd = float(np.sqrt(((g - mu) ** 2 * f).sum()))
    cdf = np.cumsum(f)
    qs = {q: float(g[np.searchsorted(cdf, q)]) for q in (0.05, 0.25, 0.5, 0.75, 0.95)}
    return dict(sd_error=se,
                observed_mean=float(x.mean()), observed_sd=float(x.std(ddof=1)),
                deconvolved_mean=mu, deconvolved_sd=sd,
                deconvolved_quantiles=qs, grid=g, density=f)
