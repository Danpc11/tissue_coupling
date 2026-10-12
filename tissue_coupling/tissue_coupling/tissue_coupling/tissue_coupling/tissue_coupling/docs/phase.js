/* Inferencia de fase circadiana a partir de una matriz de expresion.
 *
 * Implementacion en JavaScript del mismo metodo que python/lib/phase.py, para que la herramienta
 * corra enteramente en el navegador: los datos de quien la use no salen de su maquina.
 *
 * El ordenamiento define la fase SALVO rotacion y reflexion. Las horas absolutas no son
 * interpretables; lo que se reporta es fase relativa, amplitud por individuo y, con dos tejidos,
 * el acoplamiento entre ellos con un estadistico invariante a rotacion.
 */
const CLOCK_GENES = ["BMAL1","ARNTL","NPAS2","CLOCK","PER1","PER2","PER3","NR1D1","NR1D2","DBP",
  "TEF","HLF","CIART","C1ORF51","CRY1","CRY2","BHLHE40","BHLHE41","NFIL3","RORA"];

function zscoreRows(M){
  return M.map(r=>{
    const n=r.length, mu=r.reduce((a,b)=>a+b,0)/n;
    const sd=Math.sqrt(r.reduce((a,b)=>a+(b-mu)**2,0)/n)||1e-9;
    return r.map(v=>(v-mu)/sd);
  });
}
/* Dos primeras componentes por iteracion de potencia con deflacion. Basta con dos: el angulo en
   ese plano es la fase inicial. */
function topTwoComponents(Z){
  const g=Z.length, n=Z[0].length;
  const mulT=v=>{const o=new Float64Array(g);for(let i=0;i<g;i++){let s=0;for(let j=0;j<n;j++)s+=Z[i][j]*v[j];o[i]=s;}return o;};
  const mul =v=>{const o=new Float64Array(n);for(let j=0;j<n;j++){let s=0;for(let i=0;i<g;i++)s+=Z[i][j]*v[i];o[j]=s;}return o;};
  const norm=v=>{let s=0;for(const x of v)s+=x*x;s=Math.sqrt(s)||1e-12;return v.map(x=>x/s);};
  const comps=[]; let R=Z.map(r=>r.slice());
  for(let c=0;c<2;c++){
    let v=norm(Array.from({length:n},()=>Math.random()-0.5));
    for(let it=0;it<200;it++){
      const u=(()=>{const o=new Float64Array(g);for(let i=0;i<g;i++){let s=0;for(let j=0;j<n;j++)s+=R[i][j]*v[j];o[i]=s;}return norm(Array.from(o));})();
      const w=(()=>{const o=new Float64Array(n);for(let j=0;j<n;j++){let s=0;for(let i=0;i<g;i++)s+=R[i][j]*u[i];o[j]=s;}return norm(Array.from(o));})();
      let d=0;for(let j=0;j<n;j++)d+=Math.abs(w[j]-v[j]);
      v=w; if(d<1e-10)break;
    }
    const u=(()=>{const o=new Float64Array(g);for(let i=0;i<g;i++){let s=0;for(let j=0;j<n;j++)s+=R[i][j]*v[j];o[i]=s;}return Array.from(o);})();
    let sv=0;for(let j=0;j<n;j++){let s=0;for(let i=0;i<g;i++)s+=R[i][j]*u[i];sv+=s*v[j];}
    const un=norm(u);
    for(let i=0;i<g;i++)for(let j=0;j<n;j++)R[i][j]-=un[i]*sv*v[j];
    comps.push(v);
  }
  return comps;
}
function lstsqCosSin(Z,ang){
  const n=ang.length;
  let Scc=0,Sss=0,Scs=0;
  for(let j=0;j<n;j++){const c=Math.cos(ang[j]),s=Math.sin(ang[j]);Scc+=c*c;Sss+=s*s;Scs+=c*s;}
  const det=Scc*Sss-Scs*Scs||1e-12;
  return Z.map(r=>{
    let Sc=0,Ss=0;
    for(let j=0;j<n;j++){Sc+=r[j]*Math.cos(ang[j]);Ss+=r[j]*Math.sin(ang[j]);}
    return [(Sss*Sc-Scs*Ss)/det,(Scc*Ss-Scs*Sc)/det];
  });
}
/* Plantilla de fases relativas de los genes del reloj, obtenida por cosinor en una cohorte con
   horas de muestreo conocidas (biopsias seriadas de vasto lateral cada 4 h). Anclar la forma de
   onda fija el periodo en 24 h; sin eso el ordenamiento libre puede engancharse al armonico de
   12 h en un tejido y no en el otro, lo que produce una relacion de factor dos entre ambos y
   anula la correlacion circular simple. */
const TEMPLATE = {
  BMAL1:[13.6,1.07], ARNTL:[13.6,1.07], NPAS2:[14.3,1.10], CLOCK:[15.2,0.62],
  PER1:[9.2,0.88], PER2:[4.0,0.96], PER3:[12.0,0.92], NR1D1:[20.2,0.95], NR1D2:[12.0,0.72],
  DBP:[23.7,1.09], TEF:[12.0,0.83], HLF:[12.0,0.64], CIART:[12.0,0.92], C1ORF51:[12.0,0.92],
  CRY1:[0.0,0.84], CRY2:[12.0,0.66], BHLHE40:[0.0,1.11], BHLHE41:[1.0,1.11],
  NFIL3:[2.0,0.90], RORA:[6.0,0.55]
};
/* Fase con el periodo ANCLADO. Es el metodo correcto para la fase; el libre lo es para la
   amplitud. Usar uno solo para las dos cosas da una de ellas mal. */
function inferPhaseAnchored(Z,genes,grid=721){
  const idx=[],ph=[],am=[];
  genes.forEach((g,i)=>{ const t=TEMPLATE[g.toUpperCase()]; if(t){idx.push(i);ph.push(t[0]*2*Math.PI/24);am.push(t[1]);} });
  if(idx.length<8) return null;
  const n=Z[0].length, out=new Array(n);
  const T=[];
  for(let k=0;k<grid;k++){
    const t=2*Math.PI*k/grid;
    let col=idx.map((_,j)=>am[j]*Math.cos(t-ph[j]));
    const mu=col.reduce((a,b)=>a+b,0)/col.length;
    const sd=Math.sqrt(col.reduce((a,b)=>a+(b-mu)**2,0)/col.length)||1e-9;
    T.push(col.map(v=>(v-mu)/sd));
  }
  for(let j=0;j<n;j++){
    let o=idx.map(i=>Z[i][j]);
    const mu=o.reduce((a,b)=>a+b,0)/o.length;
    const sd=Math.sqrt(o.reduce((a,b)=>a+(b-mu)**2,0)/o.length)||1e-9;
    o=o.map(v=>(v-mu)/sd);
    let best=-Infinity,ba=0;
    for(let k=0;k<grid;k++){ let s=0; const e=T[k]; for(let i=0;i<o.length;i++)s+=e[i]*o[i];
      if(s>best){best=s;ba=2*Math.PI*k/grid;} }
    out[j]=ba;
  }
  return {ang:out,nGenes:idx.length};
}
/* |rho| sobre el angulo doble: delata ambiguedad de periodo residual. Con el periodo anclado
   debe quedar por debajo de ~0.12. */
function rhoDouble(a,b){ return Math.abs(rhoJS(a.map(x=>2*x),b.map(x=>2*x))); }

function inferPhase(Z,iters=60,grid=721){
  const n=Z[0].length;
  const [p1,p2]=topTwoComponents(Z);
  let ang=Array.from({length:n},(_,j)=>Math.atan2(p2[j],p1[j]));
  for(let it=0;it<iters;it++){
    const B=lstsqCosSin(Z,ang);
    const Ecols=[];
    for(let k=0;k<grid;k++){
      const t=2*Math.PI*k/grid;
      let col=B.map(b=>b[0]*Math.cos(t)+b[1]*Math.sin(t));
      const mu=col.reduce((a,b)=>a+b,0)/col.length;
      const sd=Math.sqrt(col.reduce((a,b)=>a+(b-mu)**2,0)/col.length)||1e-9;
      Ecols.push(col.map(v=>(v-mu)/sd));
    }
    const next=new Array(n);
    for(let j=0;j<n;j++){
      let o=Z.map(r=>r[j]);
      const mu=o.reduce((a,b)=>a+b,0)/o.length;
      const sd=Math.sqrt(o.reduce((a,b)=>a+(b-mu)**2,0)/o.length)||1e-9;
      o=o.map(v=>(v-mu)/sd);
      let best=-Infinity,ba=0;
      for(let k=0;k<grid;k++){
        let s=0;const e=Ecols[k];
        for(let i=0;i<o.length;i++)s+=e[i]*o[i];
        if(s>best){best=s;ba=2*Math.PI*k/grid;}
      }
      next[j]=ba;
    }
    let d=0;for(let j=0;j<n;j++)d+=Math.abs(Math.atan2(Math.sin(next[j]-ang[j]),Math.cos(next[j]-ang[j])));
    ang=next; if(d/n<1e-4)break;
  }
  return ang;
}
function amplitude(Z,ang){
  const B=lstsqCosSin(Z,ang), g=Z.length, n=Z[0].length;
  // base ortonormal del plano generado por B (Gram-Schmidt sobre las dos columnas)
  let a=B.map(b=>b[0]), b2=B.map(b=>b[1]);
  const dot=(x,y)=>x.reduce((s,v,i)=>s+v*y[i],0);
  const nrm=x=>Math.sqrt(dot(x,x))||1e-12;
  const q1=a.map(v=>v/nrm(a));
  const proj=dot(b2,q1); const r=b2.map((v,i)=>v-proj*q1[i]);
  const q2=r.map(v=>v/nrm(r));
  const out=new Array(n);
  for(let j=0;j<n;j++){
    let s1=0,s2=0;
    for(let i=0;i<g;i++){s1+=Z[i][j]*q1[i];s2+=Z[i][j]*q2[i];}
    out[j]=Math.sqrt(s1*s1+s2*s2);
  }
  return out;
}
function circMean(a){let c=0,s=0;for(const x of a){c+=Math.cos(x);s+=Math.sin(x);}return Math.atan2(s,c);}
function rhoJS(a,b){
  const ma=circMean(a),mb=circMean(b);
  let num=0,da=0,db=0;
  for(let i=0;i<a.length;i++){
    const sa=Math.sin(a[i]-ma),sb=Math.sin(b[i]-mb);
    num+=sa*sb; da+=sa*sa; db+=sb*sb;
  }
  return num/Math.sqrt(da*db||1e-12);
}
function permP(a,b,nperm=2000){
  const obs=Math.abs(rhoJS(a,b)); let ge=0;
  for(let k=0;k<nperm;k++){
    const p=b.slice();
    for(let i=p.length-1;i>0;i--){const j=Math.floor(Math.random()*(i+1));[p[i],p[j]]=[p[j],p[i]];}
    if(Math.abs(rhoJS(a,p))>=obs)ge++;
  }
  return {rho:obs,p:(ge+1)/(nperm+1)};
}
/* Control obligatorio: el ciclo recuperado debe SER el reloj. Sin el, tejidos cuyo eje ciclico
   dominante no es circadiano (sangre, hipotalamo, corteza) pasan por relojes. */
function isClock(genes,Z,ang){
  const B=lstsqCosSin(Z,ang);
  const ph={}; genes.forEach((g,i)=>ph[g.toUpperCase()]=Math.atan2(B[i][1],B[i][0]));
  const bm=ph["BMAL1"]??ph["ARNTL"], db=ph["DBP"];
  if(bm===undefined||db===undefined) return {ok:null,sep:null,reason:"BMAL1/ARNTL o DBP ausentes"};
  let d=Math.abs(Math.atan2(Math.sin(db-bm),Math.cos(db-bm)))*24/(2*Math.PI);
  return {ok:(d>=8&&d<=16),sep:d,reason:null};
}
