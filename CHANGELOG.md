# Changelog

## v1.2.0 — cadena de reproducibilidad cerrada

### Corregido
- **Faltaba el paso que ensambla `signatures_log2_common_genes.tsv`**, del que dependen los paneles
  H e I de la Figura 3 y el párrafo final de la sección 3.4. `run_all.sh` fallaba en el paso de
  figuras sobre un clon limpio. Añadido `python/analyses/23_build_signatures_table.py`, que
  regenera tres de las cuatro columnas desde los datos públicos y **comprueba cada una contra la
  referencia versionada**, abortando si alguna discrepa.
- **Rutas absolutas** de un entorno de trabajo en `fig6_mecanismo.py`, sustituidas por rutas
  relativas a `T2D_OUT`.
- **Los paneles redefinían los conjuntos de genes en el sitio**, cada uno de forma distinta y
  añadiendo BMAL1, que pertenece al oscilador central y no al brazo de salida. Ahora importan de
  `gene_sets.py`, como el propio módulo instruye. Esto cambia tres cifras del manuscrito: la
  correlación con tiazolidinediona pasa de −0,32 a −0,02, y los valores de Cushing e
  hidrocortisona estaban intercambiados entre sí.
- **Cadena de citas corrida** al eliminar una referencia sin llamado.

### Añadido
- `tests/verify_inputs.py`: comprueba que cada fichero que lee una figura lo produce el
  repositorio, está versionado en `data/reference/` o es una descarga pública declarada. Conectado
  a `run_all.sh`. Es el test que habría detectado el fallo anterior.
- `data/reference/`: entradas versionadas, entre ellas el contraste resistente frente a sensible
  de GSE13070, que no se regenera de forma consistente por el desbalance de 35 contra 6 y se
  declara así en Métodos.
- `AGENTS.md` y `python/audit_manuscript.py`: quince agentes de revisión con su criterio de
  aprobado, cinco de ellos ejecutables.
- Análisis nuevos: validación leave-one-out a dos niveles, metaanálisis de tres niveles,
  sensibilidad por componentes principales, ruido del estimador por grupo, concordancia con CHIRAL,
  panel ampliado de tejidos de GTEx.

### Archivado
- `python/figures/make_figures.py` pasa a `python/figures/_legacy/`: generaba las figuras del
  envío anterior y leía ficheros que ya no produce nadie.

## v1.1.0 — reencuadre mecanístico
Ver la entrada anterior.
# Changelog

## v1.1.0 — reencuadre mecanistico

El analisis pasa de describir el acoplamiento entre tejidos a identificar la senal que lo gobierna.

### Anadido
- `python/analyses/01b_fig1_signature_test.py` — expresion diferencial convencional con t moderada
  estilo limma implementada en el propio script, y prueba de que no recupera la firma. Cierra el
  hueco de reproducibilidad de la Fig. 1, que antes dependia de un paso manual.
- `python/analyses/05_clamp_within_person.py` — clamp y tiazolidinediona dentro de la persona.
- `python/analyses/06_meta_analysis.py` — metaanalisis de efectos aleatorios (DerSimonian-Laird
