"""Inferencia de fase y amplitud circadianas a partir de una sola muestra.

El problema: ninguna cohorte multiorgano registra la hora de la biopsia. La DIFERENCIA de fase
entre dos tejidos de una persona, en cambio, no necesita la hora absoluta.

El ordenamiento define la fase SALVO rotacion y reflexion, de modo que las horas absolutas no son
interpretables: es una pseudocuantificacion. Todas las comparaciones usan estadisticos invariantes
a rotacion, y los valores no se comparan entre cohortes con ventanas de muestreo distintas.

Controles obligatorios:
  is_clock()     el ciclo recuperado debe SER el reloj (antifase BMAL1-DBP, 8-16 h). Sangre,
                 hipotalamo y corteza fallan; sin este filtro la sangre saldria como el tejido
                 mejor acoplado, que es lo contrario de lo que es.
  align()        antes de comparar dos ordenaciones hay que alinearlas; sin alinear se mide la
                 indeterminacion de rotacion y reflexion, no el error, y se infla unas cinco veces.
  double_angle() delata ambiguedad de periodo residual.
  coupling()     P por permutacion DENTRO de la cohorte.
"""
import numpy as np, pandas as pd
from gene_sets import panel as _panel
CLK=_panel("Phase ordering")          # registro unico: python/lib/gene_sets.py
def kap(a, kmax=1e4):
    """Concentracion de von Mises y parametro de orden R.

    Con R proximo a 1 la estimacion de kappa diverge, lo que ocurre cuando dos series estan casi
    perfectamente alineadas. Se satura en kmax en lugar de devolver nan: una dispersion por debajo
    de la resolucion del metodo no es informativa y no debe propagarse como valor no numerico.
    """
    from scipy.special import i0,i1
    R=np.abs(np.mean(np.exp(1j*a)))
    if R>1-1e-8: return kmax,R
    if R<1e-8:   return 1e-6,R
    k=R*(2-R**2)/(1-R**2)
    for _ in range(60):
        A=i1(k)/i0(k)
        den=1-A**2-A/k
        if not np.isfinite(den) or abs(den)<1e-12: break
        k=min(max(k-(A-R)/den,1e-6),kmax)
    return k,R
def unsup(Z,iters=60):
    U,S,Vt=np.linalg.svd(Z-Z.mean(1,keepdims=True),full_matrices=False)
    ang=np.arctan2(Vt[1],Vt[0])
    for _ in range(iters):
        D=np.column_stack([np.cos(ang),np.sin(ang)])
        B,*_=np.linalg.lstsq(D,Z.T,rcond=None)
        grid=np.linspace(0,2*np.pi,721,endpoint=False)
        E=np.column_stack([np.cos(grid),np.sin(grid)])@B
        E=(E-E.mean(1,keepdims=True))/(E.std(1,keepdims=True)+1e-9)
        O=(Z-Z.mean(0))/(Z.std(0)+1e-9)
        new=grid[np.argmax(E@O,axis=0)]
        if np.mean(np.abs(np.angle(np.exp(1j*(new-ang)))))<1e-4: return new
        ang=new
    return ang
def align(a,b):
    best=None
    for s in (1,-1):
        bb=s*b; rot=np.angle(np.mean(np.exp(1j*(a-bb))))
        d=np.angle(np.exp(1j*(a-(bb+rot)))); R=np.abs(np.mean(np.exp(1j*d)))
        if best is None or R>best[0]: best=(R,d)
    return best[1]
def phase_of(E,genes=CLK,seed=0,nsplit=14):
    g=[x for x in genes if x in E.index]
    Z=E.loc[g]; Z=((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T.values
    a=unsup(Z); rng=np.random.default_rng(seed); ds=[]
    for _ in range(nsplit):
        pr=rng.permutation(len(g))
        ds.append(align(unsup(Z[pr[:len(g)//2]]),unsup(Z[pr[len(g)//2:]])))
    kh,_=kap(np.concatenate(ds))
    return pd.Series(a,index=E.columns), np.sqrt(1/kh)/2, g

def circ_mean(a): return np.angle(np.mean(np.exp(1j*a)))
def rho_js(a,b):
    """Jammalamadaka-Sarma. Invariante a rotacion y, en valor absoluto, a reflexion.
    Se prefiere al estadistico de rangos: con ventanas de muestreo estrechas el de rangos
    colapsa (0.008 con ventana de 3 h frente a 0.368 de este)."""
    a=a-circ_mean(a); b=b-circ_mean(b)
    return np.sum(np.sin(a)*np.sin(b))/np.sqrt(np.sum(np.sin(a)**2)*np.sum(np.sin(b)**2)+1e-12)
def coupling(pa,pb,nperm=2000,seed=0):
    """|rho| con su P por permutacion DENTRO de la cohorte. El valor absoluto no es
    comparable entre cohortes con ventanas de muestreo distintas; las comparaciones
    validas son internas."""
    r=abs(rho_js(pa,pb)); rng=np.random.default_rng(seed)
    null=np.array([abs(rho_js(pa,pb[rng.permutation(len(pb))])) for _ in range(nperm)])
    return r,(np.sum(null>=r)+1)/(nperm+1),null.mean()
def is_clock(E,ang,genes=None,tolerance_h=4.0):
    """Verifica que la estructura ciclica recuperada ES el reloj: BMAL1 y DBP deben quedar en
    antifase. Sin esta comprobacion, sangre e hipotalamo pasan por relojes.

    Devuelve (cumple, separacion en horas). La separacion entre dos fases en un ciclo de 24 h va
    de 0 a 12 h, y 12 h es la antifase exacta; el criterio es que diste de ella como mucho
    `tolerance_h`, es decir, separacion >= 8 h con el valor por defecto. BMAL1 se busca tambien
    por su simbolo anterior, ARNTL, que es el que llevan las plataformas antiguas."""
    from gene_sets import resolve
    g=[x for x in (genes or CLK) if x in E.index]
    names={k: resolve(k if k!="BMAL1" else "ARNTL",E.index) or (k if k in E.index else None) for k in ("BMAL1","DBP")}
    if names["BMAL1"] is None or names["DBP"] is None: return False,np.nan
    g=list(dict.fromkeys(g+[names["BMAL1"],names["DBP"]]))
    Z=E.loc[g]; Z=((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T
    D=np.column_stack([np.cos(ang),np.sin(ang)]); B,*_=np.linalg.lstsq(D,Z.values.T,rcond=None)
    ph=np.arctan2(B[1],B[0]); gi={x:i for i,x in enumerate(g)}
    sep=abs(np.angle(np.exp(1j*(ph[gi[names["DBP"]]]-ph[gi[names["BMAL1"]]]))))*24/(2*np.pi)
    return bool(12-sep<=tolerance_h), sep


# ---------------------------------------------------------------------------
# ORDENAMIENTO CON PERIODO ANCLADO
#
# El ordenamiento libre no fija el periodo: encuentra el ciclo dominante, y si en un tejido ese
# ciclo es el armonico de 12 h, la fase queda en otra escala y la relacion entre tejidos sale de
# factor dos. Eso produce una nube en forma de X en el plano fase-fase y anula la correlacion
# circular simple, que solo detecta pendiente +1.
#
# QUE METODO USAR PARA QUE:
#   FASE      -> order_anchored. Fija la forma de onda con una plantilla de fases relativas
#                obtenida de una cohorte con horas conocidas, de modo que el periodo queda en
#                24 h. Estable: 0.32-0.47 entre plantillas de 13 a 200 genes.
#   AMPLITUD  -> unsup + amplitude. El ordenamiento libre ajusta la onda a los datos, asi que la
#                base de proyeccion es la del tejido y la amplitud sale bien. Estable: 0.30-0.37.
#
# Usar un solo metodo para las dos cosas da una de ellas mal: el anclado desestabiliza la
# amplitud (0.07-0.38) y el libre desestabiliza la fase (0.21-0.62).
# ---------------------------------------------------------------------------

def select_template(template, n=60, exclude_unannotated=True):
    """Los n genes mas ritmicos de la plantilla, que es el conjunto que ancla la fase.

    El tamano NO es un detalle de implementacion: la estimacion es estable entre 0.32 y 0.47 para
    plantillas de 13 a 200 genes, con declive suave en las grandes al entrar ritmos debiles, y con
    13 genes queda ambiguedad de periodo residual (estadistico de angulo doble 0.23 frente a 0.06
    con 60). Usar la plantilla completa, de miles de genes, no reproduce las cifras publicadas.
    Por defecto n = 60, el valor con el que se calcularon.

    exclude_unannotated descarta identificadores LOC, LINC y ...orf, que son ritmicos en la cohorte
    de referencia pero no transfieren de forma fiable entre plataformas.
    """
    idx = list(template.index)
    if exclude_unannotated:
        idx = [g for g in idx
               if not str(g).startswith(("LOC", "LINC")) and "orf" not in str(g).lower()]
    t = template.loc[idx].sort_values("r2", ascending=False)
    return list(t.index[:n])


def order_anchored(E, template, genes=None, n_template=60, grid=721):
    """Fase con el periodo anclado a 24 h.

    template:   DataFrame indexado por gen con columnas r2, phase_h y amp, obtenido por cosinor en
                una cohorte con horas de muestreo conocidas (ver 02b_phase_template.py).
    genes:      conjunto explicito a usar. Si es None se toman los n_template mas ritmicos
                (ver select_template); NO la plantilla completa.
    """
    if genes is None:
        genes = select_template(template, n_template)
    g = [x for x in genes if x in E.index and x in template.index]
    if E.shape[1] < 3:
        raise ValueError("order_anchored estandariza cada gen entre las muestras de la cohorte: no se puede "
                         "aplicar a una muestra aislada. Para una muestra nueva use fit_reference() e "
                         "infer_from_reference(), que fijan esos parametros en una cohorte de referencia.")
    Z = E.loc[g]; Z = ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T.values
    ph = template.phase_h[g].values * 2 * np.pi / 24
    amp = template.amp[g].values
    gr = np.linspace(0, 2 * np.pi, grid, endpoint=False)
    T = np.array([amp * np.cos(t - ph) for t in gr])
    T = (T - T.mean(1, keepdims=True)) / (T.std(1, keepdims=True) + 1e-9)
    O = (Z - Z.mean(0)) / (Z.std(0) + 1e-9)
    return pd.Series(gr[np.argmax(T @ O, axis=0)], index=E.columns), g

def amplitude(E, ang, genes=None):
    """Amplitud por individuo: norma de la proyeccion sobre la base ortonormal de la forma de
    onda de la cohorte. Independiente de la fase por construccion. Usar con ordenamiento LIBRE.

    Es RELATIVA A LA COHORTE: la estandarizacion y la base se aprenden de las muestras de E, asi que
    el valor de una muestra depende de las demas. Para una muestra que no esta en la cohorte, usar
    fit_reference() e infer_from_reference(). Es una magnitud transcriptomica: no distingue menos
    amplitud en cada celula de mas dispersion de fase entre celulas."""
    g = [x for x in (genes or CLK) if x in E.index]
    Z = E.loc[g]; Z = ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T.values
    D = np.column_stack([np.cos(ang), np.sin(ang)])
    B, *_ = np.linalg.lstsq(D, Z.T, rcond=None)
    Q, _ = np.linalg.qr(B.T)
    return pd.Series(np.sqrt(((Z.T @ Q) ** 2).sum(1)), index=E.columns)

def double_angle(a, b):
    """|rho| sobre el angulo doble. Detecta relaciones axiales (modulo pi) que rho_js no ve.
    Un valor alto delata que persiste la ambiguedad de periodo; con el periodo anclado debe
    quedar por debajo de ~0.12."""
    return abs(rho_js(2 * np.asarray(a), 2 * np.asarray(b)))


# ---------------------------------------------------------------------------
# UNA MUESTRA NUEVA
#
# Los estimadores de arriba son RELATIVOS A UNA COHORTE: estandarizan cada gen entre las muestras
# que reciben y, en el ordenamiento libre, aprenden de ellas la forma de onda. "Muestra unica"
# significa una muestra por persona, no una muestra sin cohorte. Para inferir la amplitud de una
# muestra que no forma parte de la cohorte hay que fijar antes, en una cohorte de referencia, todo
# lo que el estimador aprende: los genes, la media y la desviacion de cada gen, y los coeficientes
# de la forma de onda. Es lo que hace la validacion dejando uno fuera (pasos 04 y 05), y asi el
# valor de una muestra no cambia cuando se anaden o se quitan otras.
# ---------------------------------------------------------------------------

def fit_reference(E, genes=None, ang=None):
    """Aprende de la cohorte E (genes x muestras, ya centrada por muestra) lo necesario para
    puntuar muestras nuevas: genes, media y desviacion por gen, y coeficientes de la forma de onda.
    Las fases de las muestras de la cohorte salen del ordenamiento libre, salvo que se den en `ang`.

    LOS COEFICIENTES SE AJUSTAN POR MINIMOS CUADRADOS, con termino independiente: z = m + a*cos(theta)
    + b*sin(theta). El promedio simple de z*cos(theta) y z*sin(theta), que es como se calcularon los
    valores publicados, solo coincide con ese ajuste (salvo un factor 2) si las fases de la cohorte
    estan repartidas por igual en el ciclo; con un muestreo desigual, por ejemplo una cohorte tomada
    en media jornada, mezcla los dos coeficientes y sesga la fase. Se guardan los dos: `a`, `b` son
    los de minimos cuadrados, y `a_published`, `b_published` los promedios, para poder regenerar lo
    publicado (infer_from_reference con geometry="as published").
    """
    g = [x for x in (genes or CLK) if x in E.index]
    X = E.loc[g]
    Z = ((X.T - X.mean(1)) / (X.std(1) + 1e-9)).T
    ang = unsup(Z.values) if ang is None else np.asarray(ang, float)
    D = np.column_stack([np.ones(len(ang)), np.cos(ang), np.sin(ang)])
    B, *_ = np.linalg.lstsq(D, Z.values.T, rcond=None)
    a_pub = np.array([np.mean(Z.loc[gg].values * np.cos(ang)) for gg in g])
    b_pub = np.array([np.mean(Z.loc[gg].values * np.sin(ang)) for gg in g])
    return dict(genes=g, mu=X.mean(1).values, sd=X.std(1).values + 1e-9, a=B[1], b=B[2], offset=B[0],
                a_published=a_pub, b_published=b_pub, n=X.shape[1],
                phase_uniformity=float(abs(np.mean(np.exp(1j * ang)))))


def infer_from_reference(x, ref, geometry="least squares"):
    """Fase y amplitud de UNA muestra con una referencia fija. `x` es la expresion de la muestra
    (Serie indexada por gen, centrada por muestra igual que la referencia). Devuelve
    (fase en radianes, amplitud). No usa ninguna otra muestra.

    La forma de onda de la referencia es z = a*cos(theta) + b*sin(theta), con un vector a y uno b
    por gen. Esos dos vectores NO son ortogonales ni tienen la misma norma, y eso importa:

    geometry="least squares" (por defecto). Se resuelven los dos coeficientes de z en la base
        (a, b) por minimos cuadrados; la fase es el angulo de esos coeficientes y la amplitud es la
        norma de la proyeccion de z sobre el plano que generan a y b. Recupera theta exactamente
        cuando la muestra sigue la forma de onda, y la amplitud no depende de la fase. Es la misma
        definicion de amplitud que amplitude(), que es su version relativa a la cohorte.
    geometry="as published". Las proyecciones directas z.a y z.b, que es como se calcularon los
        valores de validacion publicados. Con a y b de norma distinta la fase queda sesgada (hasta
        1.7 h en la referencia de miotubos) y la ganancia de la amplitud cambia con la fase (entre
        0.76 y 1.21 veces su media en esa referencia). Se conserva para poder regenerar esos
        valores, no para usarla.
    """
    z = (np.asarray(x.loc[ref["genes"]], float) - ref["mu"]) / ref["sd"]
    if geometry == "as published":
        ap, bp = ref.get("a_published", ref["a"]), ref.get("b_published", ref["b"])
        th = np.arctan2(np.dot(z, bp), np.dot(z, ap))
        amp = np.dot(z, ap * np.cos(th) + bp * np.sin(th)) / np.sqrt(len(ref["genes"]))
        return th % (2 * np.pi), amp
    if geometry != "least squares":
        raise ValueError("geometry debe ser 'least squares' o 'as published'")
    z = z - ref.get("offset", 0.0)                       # el termino independiente del ajuste de la referencia
    B = np.column_stack([ref["a"], ref["b"]])
    if np.linalg.cond(B) > 1e6:
        raise ValueError("la forma de onda de la referencia es degenerada: a y b son casi paralelos, "
                         "y la fase no esta definida. La cohorte de referencia no tiene estructura ciclica.")
    c, *_ = np.linalg.lstsq(B, z, rcond=None)
    Q, _ = np.linalg.qr(B)
    return float(np.arctan2(c[1], c[0]) % (2 * np.pi)), float(np.sqrt(((Q.T @ z) ** 2).sum()))


def fit_anchored_reference(E, template, genes=None, n_template=60, grid=721):
    """Lo mismo para la FASE anclada: fija en la cohorte E los genes, la media y la desviacion de
    cada gen, y deja preparada la forma de onda de la plantilla, de modo que la fase de una muestra
    nueva se obtenga sin usar ninguna otra (ver order_anchored, que es relativa a la cohorte)."""
    if genes is None:
        genes = select_template(template, n_template)
    g = [x for x in genes if x in E.index and x in template.index]
    X = E.loc[g]
    ph = template.phase_h[g].values * 2 * np.pi / 24; amp = template.amp[g].values
    gr = np.linspace(0, 2 * np.pi, grid, endpoint=False)
    T = np.array([amp * np.cos(t - ph) for t in gr])
    T = (T - T.mean(1, keepdims=True)) / (T.std(1, keepdims=True) + 1e-9)
    return dict(genes=g, mu=X.mean(1).values, sd=X.std(1).values + 1e-9, T=T, grid=gr, n=X.shape[1])


def phase_from_anchored_reference(x, ref):
    """Fase (radianes) de UNA muestra contra una referencia anclada fija."""
    z = (np.asarray(x.loc[ref["genes"]], float) - ref["mu"]) / ref["sd"]
    o = (z - z.mean()) / (z.std() + 1e-9)
    return float(ref["grid"][np.argmax(ref["T"] @ o)])
