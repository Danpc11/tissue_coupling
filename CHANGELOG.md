# Changelog

## v1.1.0 — reencuadre mecanistico

El analisis pasa de describir el acoplamiento entre tejidos a identificar la senal que lo gobierna.

### Anadido
- `python/analyses/01b_fig1_signature_test.py` — expresion diferencial convencional con t moderada
  estilo limma implementada en el propio script, y prueba de que no recupera la firma. Cierra el
  hueco de reproducibilidad de la Fig. 1, que antes dependia de un paso manual.
- `python/analyses/05_clamp_within_person.py` — clamp y tiazolidinediona dentro de la persona.
- `python/analyses/06_meta_analysis.py` — metaanalisis de efectos aleatorios (DerSimonian-Laird
  sobre Fisher-z) del eje receptor-reloj en 13 estimaciones cohorte-tejido.
- `python/analyses/07_tf_activity.py` — actividad de factores de transcripcion con regulones
  DoRothEA A-C.
- `python/analyses/08_gtex_adjusted.py` — GTEx ajustado por tiempo isquemico y RIN.
- `python/lib/gpl570.py` — anotacion de Affymetrix por Entrez en lugar del simbolo del array.
- `python/model/` — oscilador de Poincare forzado por el ritmo glucocorticoide, barrido de la
  lengua de entrainment y cohortes sinteticas de biopsia unica.
- `python/figures/paper_figures.py` — las cinco figuras del manuscrito.
- Cohortes nuevas: GSE78721, GSE26637, GSE159984, GSE66446, GSE17992, GSE129843, GSE182117,
  GSE87761, GSE11868, GTEx v11 con atributos de muestra.

### Corregido
- `04_perturbations.py` usaba `argmax(|proj|)` para la fase: con un armonico, fase y antifase
  empatan y el empate lo rompe el ruido. En datos reales el 34-36 % de las muestras cambiaban de
  fase mas de 9 h. Ahora usa el argmax con signo, igual que el estimador validado.
- Los scores se calculan sobre matrices centradas por muestra. Sin centrar, en RNA-seq la
  distribucion nula de un gen al azar frente al score de reloj tiene mediana +0,44 en lugar de
  cero, y cualquier correlacion queda inflada.
- GSE159984 estaba etiquetado como musculo; son islotes pancreaticos.
- La anotacion de GPL570 es de 2014 y el 8,3 % de sus sondas con un solo Entrez llevan un simbolo
  retirado. El emparejamiento pasa a hacerse por Entrez.
- La validacion de amplitud usaba la proyeccion de 800 genes, que no valida. Con ordenamiento
  libre sobre genes del reloj: rho 0,90 y 0,88 contra cosinor en miotubos, 0,69 in vivo.

### Retirado del manuscrito
- Menor amplitud del reloj muscular en resistencia a la insulina: no replica con los estimadores
  validados en ninguna de las dos cohortes.
- Dispersion de fase en adiposo: dependia de la regla de fase corregida.
- "La infusion de insulina concentra la fase": los desplazamientos son uniformes y del tamano del
  paso del tiempo.
- HSD11B1 como eslabon causal: no se separa de la adiposidad.

## v1.0.0
Primera version publica, acompanando el envio inicial.
