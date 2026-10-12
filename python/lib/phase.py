"""Inferencia de fase circadiana a partir de una sola muestra.

El problema: ninguna cohorte multiorgano registra la hora de la biopsia, lo que ha impedido medir
la fase entre organos en humanos. La DIFERENCIA de fase entre dos tejidos de una persona, en cambio,
no necesita la hora absoluta. Este modulo aprende la fase del propio transcriptoma.

El ordenamiento no supervisado define la fase SALVO rotacion y reflexion, de modo que las horas
absolutas no son interpretables: es una pseudocuantificacion. Todas las comparaciones usan
estadisticos invariantes a rotacion, y los valores no se comparan entre cohortes con ventanas de
muestreo distintas.

Tres controles, todos obligatorios:
  is_clock()   el ciclo recuperado debe SER el reloj (antifase BMAL1-DBP, 8-16 h). Sangre,
               hipotalamo y corteza fallan; sin este filtro la sangre saldria como el tejido
               mejor acoplado, que es lo contrario de lo que es.
  align()      antes de comparar dos ordenaciones hay que alinearlas. Compararlas sin alinear
               mide la indeterminacion de rotacion y reflexion, no el error por muestra, e infla
               la estimacion unas cinco veces.
  coupling()   P por permutacion DENTRO de la cohorte.
"""
import numpy as np, pandas as pd
CLK=["BMAL1","NPAS2","CLOCK","PER1","PER2","PER3","NR1D1","NR1D2","DBP","TEF","HLF","CIART",
     "CRY1","CRY2","BHLHE40","BHLHE41","NFIL3","RORA","PDK4","CEBPD","FKBP5","KLF15","TXNIP"]
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
def is_clock(E,ang,genes=None):
    """Verifica que la estructura ciclica recuperada ES el reloj: BMAL1 y DBP deben quedar
    en antifase (8-16 h). Sin esta comprobacion, sangre e hipotalamo pasan por relojes."""
    g=[x for x in (genes or CLK) if x in E.index]
    if "BMAL1" not in g or "DBP" not in g: return False,np.nan
    Z=E.loc[g]; Z=((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T
    D=np.column_stack([np.cos(ang),np.sin(ang)]); B,*_=np.linalg.lstsq(D,Z.values.T,rcond=None)
    ph=np.arctan2(B[1],B[0]); gi={x:i for i,x in enumerate(g)}
    sep=abs(np.angle(np.exp(1j*(ph[gi["DBP"]]-ph[gi["BMAL1"]]))))*24/(2*np.pi)
    return (8<=sep<=16), sep
