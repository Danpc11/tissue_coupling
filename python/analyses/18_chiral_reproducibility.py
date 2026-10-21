import sys, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,"python/lib"); from chiral_backend import infer_phase
CLK=["ARNTL","NPAS2","CLOCK","PER1","PER2","PER3","NR1D1","NR1D2","DBP","TEF","HLF",
     "CIART","CRY1","CRY2","BHLHE40","NFIL3","RORA","PDK4","CEBPD","FKBP5","TXNIP","KLF15"]
def conc(a,b):
    d=np.angle(np.exp(1j*(np.asarray(a)-np.asarray(b)))); return np.abs(np.mean(np.exp(1j*d)))
Em0=pd.read_pickle("results/new/GTEx_muscle.pkl"); Ea0=pd.read_pickle("results/new/GTEx_adipose_sc.pkl")
don=sorted(set(Em0.columns)&set(Ea0.columns)); g=[x for x in CLK if x in Em0.index and x in Ea0.index]
print("¿CHIRAL es determinista? GTEx, 5 corridas SIN fijar semilla, sin centrar y centrado\n")
for cen in (False,True):
    Em=Em0.sub(Em0.mean(0),axis=1) if cen else Em0; Ea=Ea0.sub(Ea0.mean(0),axis=1) if cen else Ea0
    R=[]
    for i in range(5):
        a=infer_phase(Em[don].loc[g],g)["phi"]
        b=infer_phase(Ea[don].loc[g],g)["phi"]; R.append(conc(a,b))
    R=np.array(R)
    print(f"  centrado {'sí' if cen else 'no':3s}: R = {R.round(3)}  media {R.mean():.3f}  sd {R.std():.4f}")
