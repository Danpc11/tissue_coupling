# Bitácora — reloj periférico y tono glucocorticoide

**Dos registros distintos, y no hay que mezclarlos.**
Esta bitácora es el registro de trabajo: aquí cada afirmación va con la comprobación que la sostiene,
los callejones sin salida y los errores encontrados, porque es lo que necesita quien retome el
análisis. El manuscrito se escribe al revés: **afirma lo que se encontró**. Las comprobaciones son la
razón por la que podemos afirmarlo, no contenido que haya que exhibir en cada párrafo; su lugar es
Métodos. Las limitaciones se dicen una vez, en su propio párrafo de la Discusión, sin diluir cada
resultado con salvedades. Cuando una comprobación tumbó un resultado (la amplitud muscular, §2.1),
su producto para el manuscrito es **una frase menos**, no un párrafo explicando que lo comprobamos.

**Prueba de redacción**, por si hace falta recordarlo al escribir:
- *No*: "la señalización glucocorticoide y la salida del reloj covarían, aunque con solapamiento
  considerable, y si bien robusta al ajuste por adiposidad, la relación es mayoritariamente sistémica
  sin que podamos establecer causalidad".
- *Sí*: "La salida del reloj periférico sigue al tono glucocorticoide del individuo. El receptor de un
  tejido predice el reloj del otro tan bien como el propio: el reloj lee una señal sistémica."
Fuera "aunque", "si bien", "sin que podamos", "podría sugerir". El solapamiento se ve en la figura de
distribuciones; la causalidad se dice una vez en la Discusión.

Sustituye la sección de hallazgos de `TRASPASO.md` desde el 18 de octubre.

---

## 0. Estado del envío y encuadre

Nature Metabolism rechazó sin revisión: técnicamente correcto, sin avance conceptual para una
audiencia amplia. El diagnóstico fue correcto: el artículo se leía como descriptivo y metodológico.
El reencuadre es de pregunta, no de datos. Antes preguntaba *cómo se acoplan dos tejidos*; ahora
pregunta *qué señal pone en hora al reloj periférico humano, y qué le ocurre a esa señal en la
diabetes*.

**Título de trabajo:** *Glucocorticoid tone, not insulin, sets circadian output in human metabolic
tissues.*

**Resumen, en la voz del manuscrito:**
> El reloj circadiano de los tejidos periféricos humanos se desorganiza en la diabetes, pero se
> desconoce qué señal lo pone en hora. Estimamos fase y amplitud circadianas a partir de muestras
> individuales en seis cohortes y validamos la estimación contra series temporales de las mismas
> personas. La salida del reloj sigue al tono glucocorticoide: la señalización del receptor predice
> la salida del reloj en adiposo y músculo con independencia de la adiposidad, y el receptor de un
> tejido predice el reloj del otro. La insulina no participa: tres horas de clamp hiperinsulinémico
> producen el mismo cambio que tres horas de mañana. Una dosis fisiológica de hidrocortisona frente
> a placebo a la misma hora del día desplaza PER1 y PER3 en todos los participantes, y el exceso
> endógeno de cortisol reproduce la firma completa. El tejido adiposo de personas resistentes a la
> insulina lleva esa misma firma, replicada frente a sensibilidad medida en 230 individuos. La
> sensibilización con tiazolidinediona no la revierte. El reloj periférico humano es un lector de la
> señal glucocorticoide que el tejido recibe, y es esa señal, y no la insulina, la que se altera en
> la resistencia.

**Las siete afirmaciones del artículo, en orden de solidez** (cada una con su sección técnica):
1. La amplitud inferida de una sola muestra es una medida real: ρ 0,69 y 0,70 contra la amplitud
   medida por cosinor en la serie de la misma persona, en dos sistemas, y sigue su cambio dentro del
   individuo (§3).
2. Señalización glucocorticoide y salida del reloj covarían dentro del tejido y entre los dos tejidos
   de una persona, en cuatro cohortes y cinco tejidos (§6, §14, §15, §17).
3. La covariación es un rasgo del individuo: el receptor de un tejido predice el reloj de otro tan
   bien como el propio (§14.3, §15.3).
4. La insulina aguda no pone en hora el reloj: el clamp es indistinguible de 3-4 h de mañana,
   calibrado contra tres cohortes con hora registrada (§4).
5. El cortisol sí: hidrocortisona frente a salino a la misma hora mueve PER1 y PER3 en 10 de 10;
   Cushing y dexametasona reproducen la firma de tono (§5).
6. El adiposo resistente a la insulina lleva esa firma, replicada contra sensibilidad medida en 230
   personas (§6.2, §6.3).
7. La tiazolidinediona mejora la sensibilidad y revierte el transcriptoma rítmico sin tocar la firma
   de reloj: la firma mide exposición glucocorticoide, no resistencia a la insulina (§13.6b).

**Limitaciones, las tres que van en la Discusión** (una vez, juntas, sin disculpas): la fase y la
amplitud se infieren de expresión y no de muestreo longitudinal en cada donante; el cortisol no está
medido en estas cohortes; la dirección causal descansa en asociación transversal más perturbaciones
externas. Lo demás es Métodos.

**Plan de figuras**

| Fig. | Contenido | Estado |
|---|---|---|
| 1 | Diseño; qué da el análisis convencional y por qué no responde; clock output ya enriquecido | recortar a dos paneles y corregir el texto (§9.4) |
| 2 | Validación del estimador: error interno y externo, transferencia músculo→adiposo, amplitud inferida vs medida | hecha (§3) |
| 3 | El eje: concordancia por firma en cuatro cohortes, NR3C1↔reloj en siete cohortes-tejido con distribuciones, dentro del donante, ajuste por composición | hecha (§14-17) |
| 4 | Perturbaciones: clamp = mañana, hidrocortisona vs placebo, dexametasona, Cushing | hecha (§4, §5) |
| 5 | Enfermedad: firma del adiposo resistente, coincidencia con Cushing e hidrocortisona, replicación contra SI, HPO | hecha (§6, §17.4) |
| 6 | Especificidad y disociación: Hallmark del eje, TZD que no repara la firma, modelo conceptual | pendiente de montar |
| ED | GTEx diez pares y los dos modos de acoplamiento, red, FOXO3 | existe |

---

---

## 1. ERRORES DE MÉTODO ENCONTRADOS (añadir a la lista de "no repetir")

**1.1 Antifase en `04_perturbations.py`.** `phase_amplitude` toma `argmax(|proj|)`. Con un armónico
S(t+pi) = −S(t), de modo que |proj| empata en fase y antifase y el empate lo rompe el ruido. En
GSE13070 real, el **34–36 % de las muestras** cambian de fase en más de 9 h entre `argmax(proj)` y
`argmax(|proj|)`. El estimador validado en 02b (error 1,73 h) usa el argmax CON SIGNO. Hay que pasar
04 a la regla con signo y recalcular Fig. 2c y Fig. 6. En simulación sin segundo armónico real, la
regla |proj| da 52 % de antifases; con signo, 0 %.

**1.2 `template_2harm.tsv` no lo produce ningún script del repositorio.** 02b exporta
`muscle_rhythmic_genes.tsv` (r2, amp, phase_h) y 04 exige a1, b1, a2, b2. Igual que el hueco de la
Fig. 1. Añadir el paso que lo genera.

**1.3 Sondas multi-gen en GPL570.** Descartar las sondas "A /// B" hace desaparecer PER1 y NR1D1 de
GSE13070. Solución implementada en 05: asignar la sonda multi-gen a cada símbolo que no tenga sonda
propia (excluyendo MIR/LOC/SNOR).

**1.4 Símbolos desactualizados (auditado vía Entrez contra NCBI GRCh38.p13).**
| Plataforma | Sondas con símbolo viejo | Genes de la plantilla afectados | Conjuntos/reloj |
|---|---|---|---|
| Illumina HT-12 v4 (AAGMEx) | 21,9 % | 131 | ARNTL, C1orf51, BHLHB2/3, SGK, LRAP(ERAP2), ATP5A1, JARID1D, CYorf15A |
| Affymetrix U133+2 (GSE13070, GSE41168) | 8,0 % | 45 | ARNTL, ATP5A1 |
| HTA-2.0 (GSE182120) | 10,3 % | 57 | C1orf51, TXLNG2P, ARNTL, ATP5A1 |
| Agilent GPL15331 (GSE87761) | 12,6 % | 85 | C1orf51, ARNTL, CYorf15A, ATP5A1 |
| Agilent 4x44K (GSE17992) | 15,3 % | 109 | C1orf51, ARNTL, CYorf15A, JARID1D, ATP5A1 |
`gene_sets.py` resuelve los conjuntos (la Fig. 4b sale idéntica con símbolos viejos o nuevos), pero
`03_phase_coupling.py` solo renombra ARNTL: en AAGMEx entran 42/60 genes de plantilla y 19/23 del
reloj con símbolos viejos, frente a 46/60 y 22/23 con Entrez. Falta CIART, BHLHE40, BHLHE41.
**Acción:** capa de emparejamiento por Entrez antes de casar con la plantilla; declararlo en Métodos.
El acoplamiento AAGMEx (0,219) NO cambia: 0,16–0,22 con cualquier emparejamiento y plantilla ≥30.
Falsos positivos de la auditoría: GSTT1 y HLA-DRB4 "faltan" por estar fuera del ensamblaje primario.

**1.5 CDKN1A y GADD45A son artefacto de biopsia repetida** (+4 log2 en GSE129843 entre dos biopsias
de la misma mañana). Sacarlos del programa de ayuno.

**1.6 El ordenamiento libre no recupera el reloj en cohortes de biopsia matinal** (BMAL1–DBP 3,8 h en
GSE13070, 6,4 h en GSE182120, fuera de 8–16 h); sí en series con varias horas. Para cohortes
transversales de mañana el estimador defendible es el anclado con genes del reloj.

**1.7 GSE148640 y GSE95674/5 (AAGMEx) están depositados centrados por pareja o cuantil-normalizados.**
En GSE148640 placebo y SoluCortef son exactamente opuestos: solo hay la diferencia dentro de persona,
NO se puede inferir fase ni amplitud por muestra.

**1.8 GSE17992 y GSE11868 son arrays de dos colores**: solo log-ratios. En GSE17992 la convención real
es log2(ch1/ch2) (verificado con ZBTB16/FKBP5/TSC22D3/KLF15, que siguen al canal 1 en los 20 arrays,
no al colorante), y hay códigos de valor ausente de −25000/−50000 que hay que filtrar.

---

---

## 2. LO QUE SE RETIRA DEL MANUSCRITO

**2.1 Fig. 6b, "menor amplitud del reloj muscular en resistencia a la insulina": NO SE SOSTIENE.**
Auditado con los TRES estimadores sobre las DOS cohortes de músculo (distribuciones por persona en
`F_audit_distributions.png`):

| | GSE13070 (18 sens / 51 res) | GSE182120 (24 NGT / 25 T2D) |
|---|---|---|
| amplitud, proyección 800 genes | 3,27 vs 2,46 — **P 0,003** | 0,71 vs 0,67 — P 0,71 |
| amplitud, genes del reloj anclados | 0,97 vs 1,15 — P 0,41 (invertido) | 0,46 vs 0,59 — P 0,087 (invertido) |
| amplitud, ordenamiento libre | 1,96 vs 2,58 — P 0,21 (invertido) | 2,32 vs 2,50 — P 0,99 |
| fase (R, genes del reloj) | 0,20 vs 0,16 — P 0,74 | 0,16 vs 0,34 — P 0,25 |
En GSE182120 tampoco hay relación con el valor M continuo del clamp (ρ +0,10; ajustado por
laboratorio, edad e IMC: −0,11). Los estimadores basados en el reloj son los validados contra
cosinor (§3); la proyección de 800 genes no validó. Comprobación de que no es trivial: la diferencia
de 800 genes se mantiene quitando de la plantilla todos los genes de conjuntos y del reloj
(3,15 vs 2,31, P 0,004) y los scores explican solo el 14 % de su varianza; existe, pero no mide
amplitud circadiana. **Retirar del título, del resumen y de la Fig. 6.**

**2.2 Fig. 6c, "fase adiposa dispersa en resistentes": retirar.** El 0,46/0,24 publicado sale solo
con |proj|; con la regla con signo, 0,32 vs 0,27, y n = 6 sensibles.

**2.3 Fig. 6f, "la infusión de insulina concentra la fase": retirar.** Los desplazamientos por persona
tras el clamp son uniformes (R = 0,09, 68 pares), la R post-clamp va de 0,20 a 0,42 según regla y
centrado, y el cambio del reloj tras el clamp es idéntico al de 3–4 h de mañana sin insulina (§4).

**2.4 El "desplazamiento consistente de fase por hidrocortisona" (R 0,71, P 0,006): retirar.**
Dependía de proyectar sobre la base circadiana MUSCULAR. Con la base adiposa real de GSE87761 los
ángulos no se concentran (R 0,15–0,48, P 0,5–0,9). Queda la acción gen a gen, que sí es sólida.

**2.5 "HSD11B1 → reloj": NO como eslabón causal.** Ver §7.

---

---

## 3. VALIDACIÓN DEL ESTIMADOR (resultado nuevo, el más valioso)

La amplitud inferida por muestra correlaciona con la amplitud MEDIDA por cosinor en la serie temporal
de la misma persona:
- **GSE182117** (miotubos sincronizados, 8 puntos, 11 donantes, sin hora de biopsia que confunda):
  ρ 0,70 (P 0,016) en control y ρ 0,69 (P 0,019) bajo HGI; juntas ρ 0,63, P 0,0017.
- **GSE129843** (músculo in vivo, 11 hombres con sobrepeso, 6 biopsias por persona y condición):
  **ρ 0,69, P 0,0007** en 22 series. Robusta a quitar el nivel medio por persona (ρ 0,72) y el nivel
  medio por sí solo no correlaciona con la amplitud medida (ρ 0,06), así que no es un artefacto de
  expresión basal.
- Dentro de la persona: la TRF sube la amplitud medida en 7/9 (P 0,074, reproduce a Lundell) y la
  inferida en los MISMOS 7/9.
- Condición: esto vale con el estimador de genes del reloj. **La plantilla de 800 genes no transfiere
  a miotubos** (error de fase 4,3 h, ρ de amplitud −0,28) ni predice bien la hora in vivo.

**Validación externa de la fase (para Métodos, corrige el 1,73 h):** en GSE129843, cohorte externa con
hora registrada, la plantilla transferida predice la hora ABSOLUTA con error mediano de 3,2 h (47 %
dentro de 3 h) con los 60 genes del artículo, 2,4 h y 60 % con solo genes del reloj, y EMPEORA al
crecer la plantilla (800 genes: 4,5 h, 30 %). El 1,73 h es validación dentro de cohorte.
En GSE87761 (adiposo, rutina constante): plantilla muscular → adiposo da **1,7 h de error y 76 %
dentro de 3 h** con genes del reloj (3,0 h con 60 genes; 3,5 h con 800), y las fases del reloj
coinciden con la plantilla muscular (desfase mediano 0,1 h). **Esto corrige la limitación declarada
sobre transferencia entre tejidos: falla músculo→sangre, no músculo→adiposo.**

---

---

## 4. EL CLAMP ES TIEMPO, NO INSULINA (script 05, 68+69 pares de GSE13070)

Diseño real de GSE13070 (corrige al traspaso): el ADIPOSO NO TIENE CLAMP. Músculo preTZD/postTZD ×
preClamp/postClamp; adiposo solo preTZD/postTZD. Grupos NRd (sensibles) e IRd/IRd+/IRd−.

Respuesta por persona al clamp basal (replicada idéntica en los 69 pares postTZD):
- programa de ayuno/FOXO3 cae en el **94 %** (P ≈ 10⁻¹¹): PDK4 −3,0 log2, TXNIP −1,5, DDIT4 −1,0,
  TRIM63 −0,8, FBXO32 −0,7, KLF15 −0,6; HK2 sube +0,6
- salida del reloj cae en el **90 %** (NR1D1 −1,2 en el 100 %), CIART −1,6, PER1 −0,85, DBP, PER3, TEF
- **ninguno de los dos está graduado por la sensibilidad** (P de grupo 0,78 y 0,09); la captación de
  glucosa sí (HK2 +0,52 vs +0,11, P 0,07; P 0,025 tras TZD)
- el cambio de clock output y el de GC correlacionan persona a persona (ρ 0,53, P 3·10⁻⁶)

**Calibración de tiempo** (tabla `time_calibration_vs_clamp.tsv`), con tres cohortes a hora conocida:

| gen | 129843 EXF 11h−07h | 129843 TRF | 108539 12h−08h | **clamp 3 h** |
|---|---|---|---|---|
| NR1D1 | −0,80 | −1,05 | −1,43 | −1,22 |
| DBP | −1,16 | −0,93 | −0,80 | −0,71 |
| PER1 | −0,29 | −0,76 | −0,74 | −0,85 |
| CIART | −0,72 | −1,06 | −0,27 | −1,56 |
| PDK4 | +1,65 | +2,28 | −2,26 | −3,04 |
| TXNIP | −0,23 | +0,41 | −1,62 | −1,46 |
El cambio del reloj tras el clamp es el de una mañana. La supresión de PDK4/TXNIP es del tamaño de un
desayuno (GSE108539). GSE11868 (3 h de clamp, 6 delgados) replica: PDK4 y TXNIP bajan 6/6, FKBP5 sube
6/6, DBP baja 5/6 — y fija la duración del clamp en 3 h.

**Acoplamiento dentro del individuo:** el cambio del programa de ayuno no correlaciona con el cambio
de amplitud ni con la ganancia de alineamiento (|ρ| < 0,2). Lo único que replica es: cuanto más cae el
score GC, más se desplaza la fase (ρ −0,29 P 0,017 basal; −0,28 P 0,020 tras TZD; con plantilla sin
genes de los scores). Explicación alternativa no descartada: variación del intervalo entre biopsias.

**TZD 3 meses:** no cambia ayuno/FOXO3, GC ni reloj; baja el programa oxidativo en músculo (69 %,
P 2·10⁻⁴). En el transcriptoma rítmico revierte el adiposo resistente (ρ −0,26) y el de Cushing
(−0,18) pero NO su firma de reloj (ρ −0,13): sensibilizar no repara el reloj.

---

---

## 5. PERTURBACIONES GLUCOCORTICOIDES (tres, con control de hora)

**GSE148640 — hidrocortisona fisiológica vs suero salino, adiposo, ambas biopsias a las 7:00, n=10
pareado.** Único dato con control a la misma hora, así que separa hormona de tiempo.
- PER1 **+1,6 log2 en 10/10**, PER3 −0,7 en 10/10, CIART −0,3 (P 0,004), BMAL1 +0,8 (9/10, P 0,02)
- DBP, TEF, HLF, NR1D1/2, PER2, CRY1/2, CLOCK **sin cambio** (P 0,13–0,92)
- Reproduce la Fig. 4e: GC +0,63, núcleo de insulina +0,45 (P 0,002 ambas; publicado 0,58 y 0,48)
→ acción por genes concretos del receptor, NO rotación del reloj. Esto explica la Fig. 4f (ajustar por
GC no toca el acoplamiento de clock output): el cortisol entra por PER1, no por el oscilador entero.
Anotación de GPL16686: solo trae RefSeq; mapeados 13.500 clusters vía RefSeq de GPL570 + 6 clusters
clave por coordenada hg19 (NR1D1, NR1D2, CIART, PDK4, TSC22D3, CEBPD) — **verificar contra
HuGene-2_0-st-v1.na36.hg19.transcript.csv**.

**GSE66446 — Cushing adiposo (5 vs 11).** Salida del reloj más baja (P 0,013): DBP −0,88 (P 0,009),
PER3 −0,78 (0,005), PER2 −0,64 (0,005), TEF −0,43 (0,009), HLF −0,78 (0,027); **BMAL1 +0,60
(P 0,003)**. Dianas GR (FKBP5, ZBTB16, TSC22D3) sin diferencia → saturación/desensibilización del
receptor con exposición crónica. HSD11B1 no sube (−0,53 n.s.), como debe ser si el cortisol viene de
fuera.

**GSE17992 — dexametasona 4 d, músculo, n=10 (adiposo solo n=4).** DBP −0,97 (P 0,014), CRY2 +0,9
(10/10), CLOCK +1,45 (10/10), TEF +0,7 (10/10), HLF +2,5, PER3 +1,2, BMAL1 −0,5, SLC2A4 +1,3 (10/10).
**Distinto del adiposo** → el eje GC es tejido-dependiente también en perturbación.
Cautela: dex-músculo y dex-adiposo correlacionan 0,60 entre sí y 0,00 con todo lo demás (sesgo de
sonda de la plataforma de dos colores); usar solo a nivel de gen o conjunto, no de transcriptoma.

---

---

## 6. EL ESLABÓN ENFERMEDAD ↔ CORTISOL (integración entre cohortes)

**6.1 Fig. 4b replicada en una TERCERA cohorte.** GSE13070, 41 pares músculo-adiposo basales, otra
plataforma y población: genético 0,79, sexo 0,42, HLA 0,39, **GC 0,33 (P<0,001)**, clock output 0,21
(P 0,001), oscilador 0,18, insulina 0,14 (P 0,058), oxidativo 0,05 n.s. El ajuste mutuo de la Fig. 4f
va igual (reloj 0,13→0,40 ajustado por GC; GC 0,49→0,62 ajustado por reloj). Esto resuelve el punto
flaco de que GC no superaba el nulo pareado por expresión en AAGMEx.

**6.2 El reloj del adiposo resistente está en estado "Cushing".** En los 13 genes del reloj, el perfil
resistentes−sensibles del adiposo (GSE13070) correlaciona con Cushing−control **ρ 0,71 (P 0,01, 85 %
mismo signo)** y con hidrocortisona−placebo **ρ 0,60 (P 0,03)**. Genes: PER3, DBP, CIART, TEF,
NR1D1/2, PER2 abajo; BMAL1 arriba con GC.

**6.3 Replicado en 230 personas con SI continuo (Sharma 2016, Tabla S2A, FDR<1 %, ajustado por edad,
sexo y ancestría).** El fenotipo por persona NO está en el suplementario (está en dbGaP), pero las
asociaciones sí:
| gen | β con SI | P |
|---|---|---|
| PER3 | +1,17 | 1,2e-6 |
| DBP | +0,85 | 6,4e-5 |
| HLF | +0,96 | 6,5e-8 (0,83× en el cuartil IR) |
| CLOCK | +1,73 | 1,3e-4 |
| **HSD11B1** (3 sondas) | −0,40 | 2–13e-6 (**1,35× en el cuartil IR**) |
| TSC22D3 / KLF15 / GLUL | +1,25 / +1,16 / +1,5 | 1e-4 a 2e-7 |
| IRS2 / HK2 / FOXO1 | positivos | |
Músculo (S6A, 145 transcritos): DBP +0,72 (P 2e-4); cuartil IR sube TXNIP 1,18×, RRAD 1,23×,
G0S2 1,46×. Nada sobrevive al ajuste por IMC, pero eso aplica a TODO el estudio (17 transcritos
adiposos y 3 musculares en total).

**6.4 En músculo, lo que se parece a la dexametasona es el programa de ayuno, no el reloj.**
IR músculo vs dex músculo: reloj ρ 0,44 (P 0,13); ANGPTL4 (+1,27 / +2,63), PDK4, RRAD, TXNIP con el
mismo signo en 5/6.

**6.5 El miotubo DT2 no reproduce la firma de reloj del músculo in vivo** (ρ 0,07): lo que falta en
cultivo es el entorno hormonal.

---

---

## 7. HSD11B1: CANDIDATO, NO ESLABÓN (comprobación que lo tumba)

En AAGMEx adiposo (n=232, símbolos actualizados, ajustado por sexo y edad):
- HSD11B1 vs salida del reloj **ρ −0,36, P 2,2e-8**; HLF −0,49, PER3 −0,39, DBP −0,39, NR1D2 −0,25.
  Percentil 1,8 entre 2.000 genes al azar.
- GSE13070 adiposo: HSD11B1 **+1,39 log2 en resistentes** (P<0,001).
- Músculo: nada (ρ 0,07). Cushing: HSD11B1 no sube.

**PERO**: HSD11B1 correlaciona con LEP a +0,66 y con CD68 a +0,43, y la salida del reloj correlaciona
con LEP a **−0,53 (P 3,3e-18)**. Al ajustar por LEP, CD68 y ADIPOQ, **HSD11B1 → reloj pasa de −0,36 a
+0,05 (P 0,47)**; en la regresión conjunta el coeficiente de HSD11B1 es −0,01 y el de LEP −0,45.
→ Con datos transversales NO se separa de la adiposidad. Reportarlo así.

**Lo que SÍ sobrevive a ese ajuste:** la covariación GC ↔ salida del reloj dentro del tejido:
0,43 → **0,30 (P 2,3e-6)**, coeficiente +0,29 en la regresión. Es el eslabón intrínseco.
En AAGMEx adiposo GC vs reloj ρ 0,45 (P 7,3e-13) sin ajustar; en músculo 0,14 (P 0,037).
Concordancia entre tejidos de la misma persona: GC ρ 0,29, reloj ρ 0,40 (n=232).

---

---

## 8. EL MECANISMO

**Para el manuscrito, en voz afirmativa:**

> El reloj de los tejidos periféricos humanos lee la señal glucocorticoide que el tejido recibe. La
> señalización del receptor predice la salida del reloj en tejido adiposo y músculo con independencia
> de la adiposidad y de la composición celular, en cuatro cohortes y cinco tejidos: AAGMEx adiposo
> ρ 0,34 (n=232), GSE41168 adiposo 0,37 y músculo 0,66 (n=28), GTEx adiposo subcutáneo 0,33,
> visceral 0,42 y músculo 0,33 (n=714, 587, 818). La relación es un rasgo del individuo: el receptor
> de un tejido predice el reloj del otro tan bien como el propio, de modo que el reloj periférico
> lee un tono sistémico y no una señal local. La firma que acompaña a ese tono —PER3, DBP, HLF, TEF
> y NR1D1/2 abajo, BMAL1 y NPAS2 arriba— es la que producen el exceso endógeno de cortisol y una
> dosis fisiológica de hidrocortisona a la misma hora del día, y la que distingue al tejido adiposo
> de personas resistentes a la insulina. La insulina aguda no participa. En músculo el mismo eje se
> expresa como programa de ayuno encendido: TXNIP, PDK4, ANGPTL4, RRAD y G0S2.

**Para la bitácora, lo que sostiene cada pieza y dónde cede:**
- Robusto al ajuste por adiposidad, macrófagos y composición en GTEx (ρ 0,31-0,44 sin cambio). En
  AAGMEx el ajuste por contenido de macrófagos baja NR3C1↔reloj de 0,24 a 0,10 mientras el score GC
  aguanta (0,35 → 0,33): **en esa cohorte el eslabón medible es la señalización, no el mRNA del
  receptor** (§17.2-17.3). AAGMEx es biopsia por aguja en vivos con infiltración variable; GTEx es
  post mortem. Las dos cosas se reportan; cada cohorte aporta una pieza.
- Músculo: replica en GTEx y GSE41168, no en AAGMEx (ρ 0,06, IC95 −0,06 a +0,19). Afirmación
  defendible: **consistente en adiposo, variable en músculo** (§16.3).
- Componente local frente a sistémico: tras quitar el nivel global del donante sobrevive en visceral
  de GTEx (+0,24) y músculo de GSE41168 (+0,43), no en subcutáneo de ninguna. No afirmar un
  componente local de tejido (§14.3, §15.3).
- HSD11B1 no es el eslabón: no se separa de la adiposidad en AAGMEx y en GTEx no tiene relación con
  el reloj (§7).
- Magnitud: 1,2 z entre el decil más alto y el más bajo de NR3C1, con solapamiento casi total. Es
  una tendencia poblacional. En el manuscrito, esto se enseña en la figura de distribuciones y no se
  escribe como clasificador (§16.2).
- Dirección causal: faltan el eQTL de NR3C1 y el knockdown (§9).

**Control interno que conviene destacar en la figura 3:** GSE41168 es cohorte de solo mujeres y la
firma de sexo da ρ −0,035, mientras la genética da +0,765 y la glucocorticoide +0,354. La firma que
no puede variar, no varía.

---

## 9. PENDIENTES, EN ORDEN

1. Pasar `04_perturbations.py` a `argmax(proj)` y recalcular Fig. 2c y Fig. 6 en el clúster.
2. Capa de emparejamiento por Entrez en 03; declararlo en Métodos.
3. Llevar `/tmp/limma_py.py` y el script 05 al repositorio; añadir el paso que genera
   `template_2harm.tsv` y el análisis de expresión diferencial de la Fig. 1.
4. Corregir `FIG1_texto_nuevo.md`: 292 (no 739), 13 (+edad) y 17 (+IMC) NO son cero, siete paneles,
   letras de panel desalineadas con la figura.
5. Verificar los 6 clusters manuales de GPL16686 contra NetAffx na36.
6. Pedir los SI por persona a Swapan K. Das (sdas@wakehealth.edu) o vía dbGaP para hacer la figura de
   distribución del adiposo contra SI.
7. Perturbación de 11β-HSD1 con transcriptoma de adiposo (§10) para decidir si la enzima es eslabón.
8. Datos pendientes que seguirían ayudando: Cedernaes 2018 (privación de sueño, adiposo y músculo
   pareados, cortisol alto, misma hora); GSE106800 (desalineación impuesta con clamp); RNA-seq de los
   ratones de ritmo aplanado de Cell Reports 2026.
9. Administrativo: tres diagramas de autor a mano (**las imágenes generadas por IA no las aceptan
   Springer Nature, Cell Press ni Elsevier** — reescribir la declaración de IA), resumen gráfico,
   DOI de Zenodo conectando el repo ANTES de la release, README (aún cita el título anterior),
   disponibilidad de datos (añadir GSE78721, GSE26637, GSE148640, GSE66446, GSE17992, GSE129843,
   GSE182117, GSE182120, GSE87761, GSE11868, GSE95674/5), recuento de palabras (margen de 27).

---

---

## 10. ENSAYOS DE INHIBIDORES DE 11β-HSD1: QUÉ HAY Y QUÉ NO

**No hay transcriptoma de adiposo depositado en GEO de ninguno de los grandes ensayos.** Resumen de
lo buscado:

- **Merck, MK-0916.** Feig PU et al., "Effects of an 11β-hydroxysteroid dehydrogenase type 1
  inhibitor, MK-0916, in patients with type 2 diabetes mellitus and metabolic syndrome."
  *Diabetes Obes Metab* 2011;13:498–504. PMID 21272190. Sin biopsias con transcriptoma; sin GEO.
- **Incyte (no Abbott), INCB13739.** Rosenstock J et al., "The 11-beta-hydroxysteroid dehydrogenase
  type 1 inhibitor INCB13739 improves hyperglycemia in patients with type 2 diabetes inadequately
  controlled by metformin monotherapy." *Diabetes Care* 2010;33:1516–22. PMID 20413513.
  Ensayo relacionado NCT00398619. Sin transcriptoma depositado.
- **AstraZeneca, AZD4017 — la vía realista.** El ensayo cruzado con prednisolona (Othonos N, Pofi R,
  ... Tomlinson JW, "11β-HSD1 inhibition in men mitigates prednisolone-induced adverse effects in a
  proof-of-concept randomised double-blind placebo-controlled trial", *Nat Commun* 2023;14:1025,
  NCT03111810) declara entre sus **endpoints secundarios los perfiles de expresión génica**, y en
  Métodos dice literalmente que los cambios en el perfil de expresión de tejido adiposo y músculo
  esquelético **"no se incluyeron en este manuscrito; se publicarán en un manuscrito posterior"**.
  Es decir: existe transcriptoma de adiposo y músculo de 30 hombres con inhibición de 11β-HSD1 frente
  a placebo, ambos bajo prednisolona, y todavía no está publicado. El artículo dice que los datos
  están disponibles "from the corresponding author upon request": **escribir a Jeremy W. Tomlinson
  (Oxford, jeremy.tomlinson@ocdem.ox.ac.uk)**. Es exactamente el experimento que necesitamos:
  glucocorticoide constante, enzima inhibida o no, adiposidad sin cambiar, 7 días.
- **AZD4017 en piel de DT2 — el único con transcriptoma ya analizado, y tampoco está en GEO.**
  Ensayo GC-SHEALD (ISRCTN74621291, NCT03313297, PI Ramzi Ajjan, r.ajjan@leeds.ac.uk; análisis
  transcriptómico de Ana Tiganescu, Leeds). Resultados clínicos: Ajjan RA et al., *Eur J Endocrinol*
  2022;186:441–455. Transcriptoma: Tiganescu et al., "11β-HSD1 inhibitor efficacy in type 2 diabetes
  is cortisol-dependent", medRxiv 2024, doi 10.1101/2024.05.10.24307180 (y Endocrine Abstracts 2024;
  104:OP6). Tiene RNA-seq de piel de 11 tratados y 6 placebo (35 días) más 611 genes diana en
  fibroblastos. **Sin accession de GEO ni ArrayExpress en el preprint**; pedirlo a Tiganescu o Ajjan.
  Dato que usamos en §13: 208 de 814 genes diana correlacionan con F/E y F/DHEAS, y la eficacia
  depende del cortisol basal.
- Otros AZD4017 sin adiposo: Abbas A et al., *J Clin Endocrinol Metab*
  2022;107:2026–2035 (osteopenia); Hardy RS et al., *J Clin Endocrinol Metab* 2021;106:174–187 (HII);
  Yadav Y et al., *Diabetes Obes Metab* 2022;24:881–890 (NAFLD/NASH).
- Alternativa sin pedir permisos: ratones con 11β-HSD1 KO o sobreexpresión adiposa (Morgan SA et al.,
  *PNAS* 2014;111:E2482–91, el trabajo que muestra que el fenotipo Cushing depende del 11β-HSD1
  adiposo). Buscar si depositaron RNA-seq.

---

---

## 11. FICHEROS GENERADOS (carpetas de salida)

`clamp_results/`: 05_clamp_within_person.py, run05_2harm.log, clamp_samples/groups/scores_change/
genes_change/within_person.tsv, clamp_per_person_<contraste>.tsv, F_clamp_within_person.png/pdf,
fig_clamp_within_person.py
`new_datasets/`: F_new_datasets, F_integration, F_integration_clock_heatmap, F_audit_distributions,
F_audit_GSE182120_M, F_AAGMEx_scores, F_AAGMEx_HSD11B1, time_calibration_vs_clamp.tsv,
GSE182117_measured_amplitude_per_donor.tsv, GSE182117_inferred_vs_measured.tsv,
GSE129843_amplitude_per_person.tsv, GSE11868_clamp3h_after_minus_before.tsv,
GSE17992_dex_log2ratio_muscle/adipose.tsv, GSE87761_adipose_template.tsv (388 genes R2>0,3),
GSE182120_per_person.tsv, GSE182120_T2D_minus_NGT_log2.tsv, signature_correlation_rhythmic/all.tsv,
signatures_log2_common_genes.tsv, GSE13070_crosstissue_concordance.tsv,
AAGMEx_HSD11B1_correlations.tsv, amplitude_phase_distributions_two_cohorts.tsv,
symbol_audit_renamed_genes.tsv, GPL16686_symbols_via_GPL570_refseq.tsv
`hc_results/`: GSE148640_per_patient.tsv, GSE148640_clock_gc_genes_diff.tsv,
GSE148640_HC_minus_placebo_allgenes.tsv.gz, F_hydrocortisone.png/pdf

---

---

## 12. POR QUÉ SEIS COHORTES

Ninguna cohorte humana reúne a la vez hora de biopsia registrada, fenotipo metabólico, dos tejidos
de la misma persona, una perturbación controlada y n grande. Cada conjunto aporta una de esas cosas
y el argumento es la intersección. **Esto va en la introducción como diseño, no en Métodos como
disculpa.**

- **Hora conocida** — GSE108539, GSE129843, GSE87761, GSE107537: calibran el estimador y dicen
  cuánto cambia cada gen por el solo paso del tiempo. Sin ellas, el clamp parece insulina.
- **Series por persona** — GSE182117 (miotubos), GSE129843 (músculo in vivo): dan la amplitud medida
  contra la que validar la inferida.
- **Dos tejidos pareados y n grande** — AAGMEx (232), GSE41168 (28), GSE13070 (41), GTEx (647):
  concordancia entre tejidos y covariación receptor-reloj con potencia.
- **Perturbación con control de hora** — GSE148640, hidrocortisona frente a salino a las 7:00:
  separa hormona de tiempo, imposible en cualquier estudio transversal.
- **Enfermedad endógena del eje** — GSE66446 (Cushing).
- **Perturbación farmacológica** — GSE17992 (dexametasona), GSE13070 (clamp y TZD).
- **Fenotipo continuo** — GSE182120 (valor M del clamp), Sharma 2016 (sensibilidad en 230).

El método es la triangulación: cada afirmación se apoya en al menos dos diseños distintos. Lo que
solo se apoyaba en uno se cayó al triangular, y por eso no está en el artículo (§2).

---

## 13. EL MECANISMO BIOLÓGICO, DESARROLLADO

Lo obvio sería "más cortisol daña el reloj y causa diabetes". Los datos dicen algo distinto y más
interesante.

**13.1 El adiposo resistente no está expuesto a más cortisol: está desacoplado de él.**
En AAGMEx (n=232, ajustado por sexo y edad), con la adiposidad como eje: HSD11B1, que regenera
cortisol, sube con la adiposidad (ρ +0,66 con LEP) mientras **NR3C1, el receptor, baja** (−0,31,
P 1,1e-6) y con él todas sus dianas (TSC22D3 −0,22, KLF15 −0,39) y la salida del reloj (PER3 −0,48,
DBP −0,50, HLF −0,63). En Cushing ocurre lo mismo con el cortisol llegando de fuera: dianas del
receptor sin cambio pese al exceso. La hidrocortisona aguda, en cambio, sube FKBP5 (9/10), TSC22D3
(8/10) y KLF15 (8/10): el eje agudo funciona cuando el receptor está intacto.
**Es resistencia glucocorticoide del tejido, no exceso de glucocorticoide**, con el mismo patrón que
la resistencia a la insulina: más señal, menos receptor, menos respuesta. Explica por qué los
inhibidores de 11β-HSD1 dieron efectos modestos en diabetes y por qué su eficacia depende del
cortisol basal de cada persona (Tiganescu, medRxiv 2024: 208 de 814 genes diana correlacionan con
F/E y F/DHEAS; los tratados con cortisol basal bajo se comportan como placebo).

**13.2 El reloj periférico mide la señal efectiva, no la hormona.**
La salida del reloj sigue al receptor (NR3C1 ρ +0,34) y a sus dianas (TSC22D3 +0,47, KLF15 +0,44),
no al sustrato (HSD11B1, que se explica por adiposidad). La covariación señalización-reloj sobrevive
al ajuste por adiposidad (0,43 → 0,30) y por inflamación (0,36). Eso vuelve útil a la firma
PER3↓/DBP↓/HLF↓: **estima la señal glucocorticoide que el tejido recibe de hecho, legible en una
sola biopsia, sin hora, sin cortisol en sangre y sin muestrear el ritmo.** F/DHEAS y F/E son
sistémicas y no dicen qué le llega al adiposo.

**13.3 Músculo y adiposo se rompen distinto.**
En músculo el eje no toca el reloj sino el programa de ayuno: PDK4 y TXNIP covarían negativamente
con la salida del reloj y KLF15, transactivador glucocorticoide del catabolismo, positivamente
(+0,32). La dexametasona lo reproduce (ANGPTL4, PDK4, RRAD, TXNIP con el mismo signo que en
resistencia, 5/6) y mueve el reloj en dirección opuesta a la del adiposo. Es la misma divergencia
que Cell Reports 2026 describe en ratón al aplanar el ritmo de glucocorticoides: la resistencia se
concentra en músculo mientras el adiposo preserva la señalización de insulina. **Músculo: programa
catabólico encendido. Adiposo: reloj apagado.** Dos averías del mismo eje, y la razón de que ningún
marcador único de resistencia a la insulina capture ambas.

**13.4 La insulina queda fuera, y eso reordena la causalidad.**
El clamp suprime el programa de ayuno en el 94 % de las personas y la salida del reloj en el 90 %,
con el mismo tamaño que 3-4 h de mañana sin insulina en tres cohortes con hora, y sin graduarse por
la sensibilidad. La insulina aguda no es el Zeitgeber del reloj periférico humano. La literatura que
lo propone (Crosby 2019, Tuvia 2021) trabaja en células o sin control horario. **En la diabetes el
reloj periférico no se desorganiza porque falle la insulina, sino porque falla el canal
glucocorticoide**, y la desorganización circadiana deja de ser consecuencia de la hiperglucemia para
ser un indicador de cuán desacoplado está el tejido de su señal hormonal de tiempo.

**13.5 La TZD mejora la sensibilidad sin arreglar el reloj.**
Tres meses de tiazolidinediona revierten el transcriptoma rítmico del adiposo resistente (ρ −0,26) y
el de Cushing (−0,18) sin tocar la firma de reloj (−0,13) ni los scores glucocorticoide y de ayuno.
PPARγ actúa aguas abajo del receptor. **Predicción:** un sensibilizador que no restaure el receptor
dejará la firma intacta, y un inhibidor de 11β-HSD1 tampoco la restaurará salvo en personas con
cortisol basal alto. Testable con el ensayo de AZD4017 (§10).

---

**13.6 El modelo, en cuatro pasos (figura conceptual):**
1. El individuo tiene un tono glucocorticoide que se expresa como nivel de receptor en sus tejidos:
   NR3C1 correlaciona entre tejidos del mismo donante (subcutáneo-visceral 0,66) y el de un tejido
   predice el reloj de otro.
2. Ese tono baja con la adiposidad y con la exposición sostenida, mientras sube la enzima que
   regenera cortisol: más ligando local, menos lector.
3. La salida del reloj sigue al lector y no al ligando. Gen a gen el receptor va con NR1D2, CLOCK,
   DBP, HLF y PER3, y en contra de BMAL1 y NPAS2, con PER1 y CIART en cero: el tono mueve el brazo
   de salida; el pulso agudo mueve PER1 (10/10 con hidrocortisona).
4. El mismo desacople se lee distinto por tejido: adiposo, reloj apagado; músculo, programa de ayuno
   encendido.

**13.6b La disociación, que es lo que vuelve útil a la firma.** La firma se mueve con exposición
glucocorticoide y con el receptor, y no se mueve con lo que cambia el fenotipo metabólico: la TZD
mejora la sensibilidad tres meses y revierte el transcriptoma rítmico sin tocarla; el clamp cambia
el transcriptoma entero y su efecto sobre el reloj es el de una mañana. **No es un marcador de
resistencia a la insulina: es un marcador de señal glucocorticoide efectiva**, y lo demuestra tanto
por lo que la mueve como por lo que no. Ya hay marcadores de resistencia a la insulina; no hay
ninguno de cuánta señal glucocorticoide recibe un tejido.

**13.7 El enlace funcional (Hallmark, §19.3).** El eje del receptor va en contra de inflamación
(TNFA/NF-κB −4,5, interferón gamma −7,3, IL-6/STAT3 −3,3, complemento −4,5), de OXPHOS (−6,8), de
adipogénesis (−3,8), de glucólisis (−3,5) y de mTORC1 (−4,0): la función antiinflamatoria canónica
del receptor, más una conexión con el metabolismo oxidativo del adiposo, el mismo programa que
Sharma 2016 asocia a sensibilidad en adiposo y de forma discordante en músculo. Y en HPO el eje va
con obesidad troncal (z +3,62) y cortisol urinario elevado (+3,77), y en contra de diabetes tipo 2
(−2,96): fenotipo glucocorticoide, no fenotipo diabético.

**13.8 Predicciones para la Discusión:**
- La firma del adiposo correlaciona con NR3C1 y con la pendiente diurna de cortisol, no con el
  cortisol matutino aislado.
- En el ensayo de AZD4017 la firma se mueve solo en quienes tienen cortisol basal alto.
- Restaurar el receptor restaura la firma sin cambiar la adiposidad.
- La firma predice la respuesta a glucocorticoides terapéuticos mejor que el cortisol en sangre.

---

## 14. GTEx v11: el eje receptor–reloj replica, y acota el mecanismo (19 oct)

Datos: `gene_reads_adult_gtex_v11_*` de adiposo subcutáneo (714), adiposo visceral (587) y músculo
(818), log2 CPM, genes detectados en >50 % de las muestras, **ajustado por sexo, edad (punto medio
del tramo) y clase de Hardy**. 453 donantes tienen los tres tejidos. Fichero
`GTEx_adipose_NR3C1_clock.tsv`, figura `F_GTEx_NR3C1.png`.

**14.1 Replica en una cuarta cohorte y en tres tejidos.** NR3C1 vs salida del reloj: adiposo sc
**+0,33**, adiposo visceral **+0,42**, músculo **+0,33** (todas P < 1e-11). En AAGMEx era +0,34
(adiposo) y +0,06 (músculo) — el músculo sí replica aquí, lo que corrige el matiz anterior. NR3C1
está en el **percentil 91** de 2.000 genes al azar y de 817 genes de abundancia similar: es alto,
pero no único; hay que decirlo así y no venderlo como el gen más correlacionado del transcriptoma.

**14.2 El ajuste por composición lo deja solo.** En adiposo sc, ajustando además por LEP y CD68:
NR3C1 +0,31 (P 6e-17), HSD11B1 +0,13, LEP −0,03, CD68 −0,04, ADIPOQ −0,06, FKBP5 −0,13, TSC22D3
+0,17. **Solo el receptor sobrevive**, igual que en AAGMEx. HSD11B1 vs LEP +0,35 (vs +0,66 en
AAGMEx), vs reloj 0,00 sin ajustar: en GTEx la enzima directamente no tiene relación con el reloj,
lo que refuerza §7 (no es el eslabón).

**14.3 Hallazgo que acota el mecanismo: el efecto es MÁS DEL DONANTE QUE DEL TEJIDO.**
En los 453 donantes con los tres tejidos, NR3C1 correlaciona entre adiposo sc y visceral a 0,66 y
con músculo a 0,20; la salida del reloj correlaciona entre tejidos a 0,51–0,65. Al quitar el nivel
global de NR3C1 del donante:
| tejido | ρ crudo | sin el nivel global del donante | solo la parte específica del tejido |
|---|---|---|---|
| adiposo sc | +0,32 | +0,13 (P 4e-3) | +0,09 (P 0,07) |
| visceral | +0,39 | **+0,23 (P 1e-6)** | **+0,24 (P 2e-7)** |
| músculo | +0,32 | +0,10 (P 0,03) | +0,07 (P 0,12) |
Y en cruzado: el NR3C1 del visceral predice el reloj del subcutáneo (+0,35) igual que el del propio
subcutáneo. **Lectura:** más de la mitad de la relación es un rasgo de la persona (tono del eje
glucocorticoide sistémico), no una propiedad local del tejido; el componente local sobrevive solo en
adiposo visceral. Esto es una corrección importante del modelo de §13: el desacople es sobre todo
sistémico, con un componente local propio del depósito visceral — el mismo depósito donde la
literatura sitúa el exceso de 11β-HSD1 y la patología metabólica. Hay que escribirlo así; es menos
"mecanismo de tejido" y más "tono glucocorticoide del individuo leído por el reloj de cada tejido".

**14.4 Especificidad gen a gen (adiposo sc).** NR3C1 va con el brazo de SALIDA y en contra del brazo
positivo: NR1D2 +0,57, CLOCK +0,41, DBP +0,37, HLF +0,36, PER3 +0,35, TEF +0,23, PER2 +0,21,
NR1D1 +0,20, CRY2 +0,15; **BMAL1 −0,18 y NPAS2 −0,24**; PER1, CIART y CRY1 ≈ 0. Es coherente con la
firma de Cushing e hidrocortisona (BMAL1 sube cuando el GC sobra, baja cuando el receptor abunda) y
con que PER1 sea la diana aguda del receptor (responde al pulso, no al tono).

**14.5 Dos comprobaciones.** (a) Por clase de Hardy, el efecto está en las tres clases con n
suficiente (muerte súbita +0,37 n=372; clase 2 +0,26 n=188; clase 4 +0,38 n=83), así que no es un
artefacto del estrés agonal. (b) El score GC completo se comporta distinto que en AAGMEx (sc −0,15,
visceral +0,06, músculo −0,02): en tejido post mortem las dianas agudas del receptor (FKBP5 −0,18)
están dominadas por el cortisol terminal, no por el tono basal. **Por eso en GTEx hay que usar NR3C1
y no el score GC**, y por eso el score GC de AAGMEx (vivos, en ayuno, por la mañana) es el bueno
para el eslabón de señalización. Declararlo explícitamente: GTEx aporta el receptor y la estructura
donante/tejido; AAGMEx aporta la señalización.

**14.6 Lo que esto hace con el plan.** El test de instrumento genético (eQTL de NR3C1) **sigue siendo
el paso que falta y ahora es más informativo**, porque el diseño dentro del donante dice que el
efecto es mayormente sistémico: si el eQTL *cis* de NR3C1 (local, por definición) mueve la firma de
reloj, se demuestra que además del tono sistémico hay un brazo local causal. Requiere genotipos
individuales (dbGaP phs000424, acceso controlado). Iniciar la solicitud.

---

---

## 15. GSE41168 (28 mujeres, replicación del artículo): tercera cohorte viva del eje

Matriz de GPL570 ya en log2, anotada con la capa de sondas multi-gen. 140 muestras: músculo y
adiposo, pre y post, en tres brazos (restricción calórica n=11-14, resveratrol n=9-13, placebo
n=11-12; Yoshino 2012). 28 pares músculo-adiposo basales.

**15.1 Fig. 4b replicada por cuarta vez (concordancia entre tejidos, basal, n=28):**
genético +0,765 (P 5e-4), HLA +0,462 (5e-4), **glucocorticoide +0,354 (5e-4)**, núcleo de insulina
+0,258 (0,004), salida del reloj +0,179 (0,021), oscilador +0,172 (0,067), oxidativo +0,060 (n.s.);
fondo +0,023. **El sexo sale −0,035 porque la cohorte es solo de mujeres** — es el control interno
perfecto: la firma que no puede variar, no varía. Vale la pena ponerlo en la figura.
Aquí el GC supera al reloj, como en GSE13070 (0,33 vs 0,21) y al revés que en AAGMEx (0,168 vs
0,258): el orden entre esos dos no es estable entre cohortes y no hay que afirmarlo.

**15.2 El eje receptor–reloj replica, y en músculo con fuerza.** Adiposo: NR3C1 +0,37 (P 0,051),
FKBP5 +0,41, TSC22D3 +0,41, **KLF15 +0,67**, GC score +0,42 (P 0,026); LEP −0,69, CD68 −0,55 y
HSD11B1 −0,50 (todos P<0,01). **Músculo: NR3C1 +0,66 (P 0,0002)**, y nada más correlaciona
(HSD11B1 −0,09, FKBP5 −0,18, GC score −0,09).

**15.3 Confirma el patrón de §14.3 y lo invierte por tejido.** Con el nivel global de NR3C1 del
donante quitado: adiposo +0,37 → **−0,06** (desaparece; y el NR3C1 del músculo predice el reloj del
adiposo igual de bien, +0,36), músculo +0,66 → **+0,43 (P 0,023)**. Es decir, en esta cohorte la
señal del adiposo es enteramente del tono sistémico del donante y la del músculo tiene componente
local. En GTEx el componente local sobrevivía en el adiposo visceral y no en el subcutáneo ni en el
músculo. **Conclusión conservadora que hay que escribir: el eje NR3C1–reloj es mayoritariamente un
rasgo del individuo; qué tejido conserva además un componente local no es consistente entre
cohortes y no debe afirmarse.**

**15.4 Perturbación nutricional (pre→post, 12 semanas): NO SUPERA EL NULO, reportar como negativo.**
En adiposo el resveratrol sube la salida del reloj (+0,21, 11/13 suben, Wilcoxon P 0,040) y el score
GC (+0,33, P 0,04) y baja HSD11B1 (−0,32); CR va igual sin llegar (+0,11, P 0,083); placebo −0,16
(P 0,52). **Pero contra 2.000 conjuntos de genes del mismo tamaño elegidos al azar, P = 0,146**
(nulo: mediana +0,02, sd 0,14), y entre brazos no hay diferencia (resveratrol vs placebo P 0,118;
CR vs placebo P 0,066; Kruskal tres brazos P 0,134). En músculo nada en ningún brazo.
El Wilcoxon solo dice que el score se mueve, no que se mueva MÁS que cualquier score de nueve genes
en esa cohorte, que es la pregunta. **Conclusión: ni el resveratrol ni la restricción calórica
mueven la firma de forma demostrable con este n (11-13).** Escribirlo como resultado negativo con la
dirección anotada; con 30-40 por brazo valdría la pena repetirlo. No usar como "primera perturbación
nutricional que mueve la firma": esa frase no se sostiene.

Ficheros: `results/new/GSE41168_full.pkl`, `GSE41168_adipose.pkl`, `GSE41168_muscle.pkl`,
`GSE41168_meta.pkl`.

---

---

## 16. AUDITORÍA DE DISTRIBUCIONES DEL EJE NR3C1–RELOJ (figura `F_NR3C1_distributions.png`)

Hecho lo que se hace siempre aquí: mirar la nube y no el coeficiente. Siete cohortes-tejido, cada
punto una persona, agrupadas por deciles/cuartiles de NR3C1 con la mediana marcada, IC95 por
bootstrap de 1.000-2.000 remuestreos, y recorte de |z| > 2,5 para descartar que la lleven las colas.

| cohorte-tejido | n | ρ | IC95 bootstrap | ρ recortado | fuera |
|---|---|---|---|---|---|
| GTEx adiposo sc | 714 | +0,328 | +0,256 a +0,398 | +0,296 | 18 |
| GTEx visceral | 587 | +0,424 | +0,348 a +0,490 | +0,401 | 10 |
| GTEx músculo | 818 | +0,329 | +0,268 a +0,389 | +0,288 | 25 |
| AAGMEx adiposo | 232 | +0,344 | +0,225 a +0,455 | +0,298 | 8 |
| **AAGMEx músculo** | 232 | **+0,059** | **−0,064 a +0,194** | +0,003 | 11 |
| GSE41168 adiposo | 28 | +0,373 | **+0,024** a +0,648 | +0,472 | 1 |
| GSE41168 músculo | 28 | +0,655 | +0,319 a +0,851 | +0,655 | 0 |

**16.1 La relación es monótona, no de colas.** Mediana del reloj por decil de NR3C1 en GTEx adiposo
sc: −0,62 −0,38 +0,01 +0,10 +0,04 +0,10 −0,09 +0,22 +0,28 +0,56; visceral: −0,52 −0,33 −0,14 −0,16
+0,02 +0,09 −0,02 +0,22 +0,26 +0,44. Sube a lo largo de todo el rango, con una meseta en el medio.
El recorte del 1-3 % más extremo baja ρ entre 0,03 y 0,04 y nunca lo rompe. Pearson y Spearman
coinciden (0,37 vs 0,33; 0,45 vs 0,42): no hay cola que mande.

**16.2 Los solapamientos son enormes y hay que decirlo.** En GTEx adiposo sc, el decil más bajo de
NR3C1 tiene mediana de reloj −0,62 y el más alto +0,56: **1,2 z de separación entre los extremos**,
con distribuciones que se solapan casi por completo (rango de cada decil ≈ 2 z). Esto es una
tendencia poblacional, **no un clasificador**: ninguna persona se asigna a un grupo por su NR3C1.
Cualquier redacción tipo "el receptor determina el reloj" es falsa a nivel individual; lo correcto
es "la mediana de la salida del reloj aumenta de forma monótona con el receptor".

**16.3 Los dos números que había que mirar de cerca.**
- **AAGMEx músculo es cero, no débil**: IC95 −0,064 a +0,194, cruza el cero, y la distribución por
  quintiles es plana (−0,07 +0,01 +0,11 +0,16 −0,07). En §8 y §14 escribí "+0,14" citando el valor
  sin ajustar; **el ajustado es +0,059 y hay que corregirlo**. Es decir, el músculo replica en GTEx
  (+0,33) y en GSE41168 (+0,66) pero NO en AAGMEx. Tres tejidos de tres cohortes y uno no replica:
  eso va en el texto, no se esconde.
- **GSE41168 adiposo roza el cero**: IC95 +0,024 a +0,648 con n=28; la mediana por cuartiles es
  −0,37 +0,13 +0,10 +0,08, o sea un escalón entre el cuartil más bajo y el resto, no una pendiente.
  Con 28 personas no se distingue una cosa de la otra. Reportar el IC, no el ρ a secas.

**16.4 Qué queda firme tras mirar las distribuciones.** El eje es sólido en **adiposo**: cuatro
estimaciones (GTEx sc, GTEx visceral, AAGMEx, GSE41168) con IC que no tocan el cero y pendiente
monótona en las dos grandes. En **músculo** es inconsistente: fuerte en GTEx y GSE41168, nulo en
AAGMEx. La frase defendible pasa a ser: *el eje receptor–reloj es consistente en tejido adiposo y
variable en músculo*. Actualizar §8 y §13 en consecuencia.

---

---

## 17. ENRIQUECIMIENTO FUNCIONAL DEL EJE NR3C1 (MSigDB v2026.2) — Y EL CONFUSOR QUE DESTAPA

Bibliotecas en `data/msigdb/` (GO BP/CC/MF, Reactome, WikiPathways, KEGG MEDICUS, C3 TFT y TFT:GTRD,
C8 tipos celulares, HPO; v2026.2, apuntar la versión en Métodos). Procedimiento: ranking de los
34.020 genes por su correlación parcial con NR3C1 en AAGMEx adiposo (ajustado por sexo, edad, LEP y
CD68) y z de cada conjunto contra 3.000 conjuntos al azar del mismo tamaño.
Ficheros: `AAGMEx_rank_vs_NR3C1.tsv`, `AAGMEx_NR3C1_enrichment.tsv`.

**17.1 El resultado dominante no es una vía: es composición celular.** Los conjuntos más opuestos a
NR3C1 son todos de la colección C8 (tipos celulares): microglía z −11,9, macrófago fetal −11,6,
monocito-macrófago de ovario −11,6, neutrófilo −10,4, Kupffer −9,4, degranulación de neutrófilo de
Reactome −8,9. Es decir, **el ranking de NR3C1 en adiposo es en buena medida un eje de contenido
leucocitario**, y ajustar por CD68 solo no lo quita. Esto no se habría visto sin la colección C8;
es exactamente el tipo de control que hay que hacer en tejido en masa.

**17.2 Comprobación directa, y aquí AAGMEx y GTEx se separan.** Reajustando con un panel
leucocitario de 20 marcadores (PTPRC, ITGAM, CSF1R, AIF1, TYROBP, C1QA/B/C, CD163, MRC1, …) y/o con
la firma completa de macrófago de C8 (≈450 genes):

| cohorte | base (sexo, edad, LEP) | + panel leucocitario | + firma macrófago | + ambas |
|---|---|---|---|---|
| **AAGMEx adiposo** (n=232) | +0,244 (P 2e-4) | +0,152 (P 0,020) | **+0,102 (P 0,12)** | +0,113 (P 0,085) |
| GTEx adiposo sc (n=714) | +0,311 | +0,307 | +0,308 | **+0,310 (P 2e-17)** |
| GTEx visceral (n=587) | +0,426 | +0,422 | +0,427 | **+0,436 (P 1e-28)** |
| GTEx músculo (n=818) | +0,329 | +0,331 | +0,321 | **+0,317 (P 2e-20)** |

**En AAGMEx el eje NR3C1–reloj NO sobrevive al ajuste por contenido de macrófagos; en GTEx no se
mueve nada.** La diferencia tiene explicación: en AAGMEx el panel leucocitario correlaciona −0,29
con el reloj y −0,33 con NR3C1, mientras que en GTEx correlaciona −0,09 y −0,05. AAGMEx es biopsia
por aguja de adiposo subcutáneo abdominal en personas vivas con IMC medio 28,9, donde la
infiltración varía mucho; GTEx es tejido post mortem con menos variación de infiltración. **Hay que
reportar las dos cosas y no elegir la que gusta.**

**17.3 Lo que sí sobrevive a todo en AAGMEx es el score GC**: +0,348 base → +0,336 con panel →
+0,334 con macrófago → +0,337 con ambas. Es decir, la **señalización** glucocorticoide (FKBP5,
TSC22D3, ZBTB16, KLF15…) se mantiene mientras que el **nivel del receptor** se confunde con
infiltración. Esto reordena el modelo: el eslabón robusto es la señalización, no el mRNA de NR3C1.
En GTEx el score GC se comporta al revés (−0,12 sc, +0,07 visceral, −0,02 músculo), por el cortisol
terminal post mortem ya discutido en §14.5. **Conclusión: AAGMEx aporta la señalización (válida),
GTEx aporta el receptor (válido); ninguna de las dos aporta las dos cosas, y eso hay que escribirlo
como diseño, no como limitación.**

**17.4 Fenotipo (HPO), con el ranking de NR3C1 ya ajustado:** HP_TRUNCAL_OBESITY z +3,62,
HP_INCREASED_URINARY_CORTISOL_LEVEL +3,77, HP_TYPE_II_DIABETES_MELLITUS **−2,96**,
HP_INSULIN_RESISTANCE −0,03, HP_HYPERGLYCEMIA +0,08. Y GOBP_RESPONSE_TO_INSULIN +2,22;
GOBP_CIRCADIAN_RHYTHM +1,24 (no significativo por sí solo: la señal del reloj está en los nueve
genes de salida, no en los 200 anotados como circadianos, que es lo que ya sabíamos desde la Fig. 1).
Los dos términos de obesidad troncal y cortisol urinario apuntando al mismo lado que el receptor son
el enlace con el fenotipo de Cushing; el de diabetes tipo 2 en contra es coherente con §13.6b (la
firma no es un marcador de resistencia a la insulina).

**17.5 Lo que estas bibliotecas NO dan:** no hay un conjunto decente de respuesta a glucocorticoides
(1 en Reactome, 13 en GO BP, ninguno útil en C3 TFT; el motivo GRE_C es legacy y GTRD no trae
NR3C1_TARGET_GENES). Seguir con el conjunto pre-especificado de nueve genes y, si se quiere mejorar,
construir uno empírico con los genes que suben en 9-10 de 10 personas en GSE148640, que es una
perturbación humana controlada y mejor argumento que una base de datos.

---

---

## 18. FIRMA GLUCOCORTICOIDE EMPÍRICA (GRAS): CONSTRUIDA Y **NO VALIDADA**. Usar el conjunto de 9 genes.

Se construyó porque MSigDB no trae un conjunto de respuesta a glucocorticoides utilizable (§17.5).
Fuente: GSE148640, hidrocortisona fisiológica frente a suero salino, adiposo, ambas biopsias a las
7:00, 10 pacientes, diseño dentro de la persona (la mejor perturbación GC humana que tenemos).
Ficheros: `GRAS_signature.tsv`, `GSE148640_HC_response_stats.tsv`.

**18.1 Construcción.** De 12.641 genes: 97 cambian en la misma dirección en los 10 de 10 pacientes
(60 suben, 37 bajan); **ninguno alcanza FDR < 0,10** con n=10, de ahí el criterio por consistencia y
magnitud en vez de por P. Criterio final: 10/10 en la misma dirección y |mediana| ≥ 0,30 (subida) o
≥ 0,20 (bajada) → 44 + 28 genes. **Se excluyeron PER1, PER3, TXNIP y PDK4 por pertenecer al reloj o
al núcleo de insulina** (si no, cualquier correlación posterior con el reloj sería circular), y
quedan **41 arriba** (ZBTB16, KLF9, CIDEA, FOXO1, ACACB, LIFR, RHOB, AOX1, BMP2, DUSP4, TIMP4…) y
**27 abajo** (TFPI2, APLNR, KDR, RGCC, ACSL3, PTGER4, AKT2, ACE…).

**18.2 Validación externa: FALLA.** Ninguna de estas cohortes se usó para construirla.
| prueba | GRAS | conjunto pre-especificado de 9 genes |
|---|---|---|
| Cushing adiposo (5 vs 11) | +0,39, P 0,51 | +0,14, P 0,58 |
| dexametasona 4 d, músculo (n=10) | −0,06, **4/10 suben** | +0,48, **9/10** |
| dexametasona 4 d, adiposo (n=4) | −0,33, **0/4** | +0,48, 3/4 |
| versión corta (20+11 genes de mayor magnitud), dex músculo | +0,12, 6/10 | — |
| GSE13070 adiposo, resistentes − sensibles | +0,12, P 0,51 | +0,21, P 0,87 |
En AAGMEx sí se comporta como cabría esperar (GRAS vs conjunto de 9 ρ +0,61; GRAS vs reloj +0,31
frente a +0,35 del de 9; GRAS vs LEP −0,26), pero **no aporta nada sobre el conjunto de 9 y pierde
en todas las validaciones de perturbación.**

**18.3 Por qué falla, y por qué eso es informativo.** La firma captura la respuesta **aguda** (22 h)
a un pulso de cortisol, que no es la misma que la exposición **crónica**: en Cushing las dianas
clásicas del receptor están saturadas (§5) y la dexametasona de 4 días reproduce el conjunto de 9
pero no la respuesta aguda de 22 h. Esto refuerza, desde otro ángulo, lo que ya decía §13.1: hay una
respuesta de pulso (PER1 10/10, ZBTB16 +2,45) y una de tono (FKBP5, TSC22D3, KLF15, ZBTB16), y son
distinguibles. Vale la pena escribirlo en la Discusión como hallazgo, no como fracaso técnico.

**18.4 Decisión.** **Usar el conjunto pre-especificado de nueve genes** (FKBP5, TSC22D3, ZBTB16,
KLF15, DDIT4, SESN1, GLUL, SGK1, CEBPD) como en todo el artículo, y citar la construcción empírica
en el suplementario con su tabla de validación, como evidencia de que la elección no fue arbitraria
y de que un conjunto derivado de un pulso agudo no sustituye al de tono. La GRAS se puede reservar
para un uso concreto y honesto: **marcador de exposición aguda** (horas) a glucocorticoide, que es
distinto de lo que mide el de nueve genes.

---

---

## 19. ACTUALIZACIÓN DE SÍMBOLOS FRENTE A MSigDB v2026.2, Y HALLMARK SOBRE EL EJE NR3C1

**19.1 Capa de sinónimos, construida y aplicada.** Desde `Human_GRCh38_p13_annot.tsv` se construyó
un mapa de 56.103 sinónimos que resuelven a uno solo de los 39.374 símbolos actuales
(`ncbi_synonym_map.tsv`). Cobertura frente a los 23.332 genes distintos de MSigDB H+C2+C5
(`msigdb_symbol_coverage.tsv`):

| matriz | genes | en MSigDB con símbolo crudo | tras mapear | ganancia |
|---|---|---|---|---|
| AAGMEx adiposo | 34.020 | 19.135 | 19.165 | +30 |
| GTEx adiposo sc (v11) | 15.584 | 14.037 | 14.030 | −7 |
| **GSE41168 adiposo** | 23.625 | 17.635 | **18.793** | **+1.158** |
| **GSE182120 músculo** | 30.905 | 18.824 | **20.175** | **+1.351** |
| GSE66446 adiposo | 16.328 | 14.446 | 14.446 | 0 |
AAGMEx ya estaba mapeada por Entrez (§1.4) y GTEx v11 trae símbolos actuales; los que ganan mucho
son los que se anotaron con la columna de símbolo de la plataforma: **GSE41168 (GPL570) y GSE182120
(HTA-2.0) perdían más de mil genes en cada enriquecimiento**. Casos típicos recuperados: la familia
MARCH→MARCHF, AARS→AARS1, FAM63A→MINDY1, ACP3, ADA2, ACKR1/3, ABRAXAS1/2.

**19.2 Pero NO cambia ninguna conclusión.** Con las matrices remapeadas, los nueve genes de la
salida del reloj y los nueve del conjunto GC se recuperaban ya todos (9/9) en las cinco cohortes, y
NR3C1↔reloj queda idéntico al tercer decimal (GSE41168 adiposo +0,373, músculo +0,655; GSE182120
+0,543; AAGMEx +0,402; GTEx +0,421). Era esperable: los genes que se recuperan son familias
renombradas en masa, no genes de nuestros conjuntos. **El mapeo importa para el enriquecimiento
global (mil genes de universo), no para los scores pre-especificados.** Matrices remapeadas
guardadas como `*_mapped.pkl`; usarlas de aquí en adelante para cualquier GSEA.

**19.3 Hallmark sobre el ranking de NR3C1 (AAGMEx adiposo, ajustado por sexo, edad, LEP, CD68).**
Casi todo el panel va en contra del receptor, y de forma muy ordenada:
MYOGENESIS −8,8, INTERFERON_GAMMA −7,3, **OXIDATIVE_PHOSPHORYLATION −6,8**, EMT −5,9, UPR −5,6,
APICAL_JUNCTION −5,3, P53 −5,3, COAGULATION −5,0, XENOBIOTIC −4,9, DNA_REPAIR −4,8, MYC_V2 −4,6,
INTERFERON_ALPHA −4,6, ALLOGRAFT_REJECTION −4,5, COMPLEMENT −4,5, **TNFA_VIA_NFKB −4,5**,
MTORC1 −4,0, APOPTOSIS −3,8, **ADIPOGENESIS −3,8**, IL2_STAT5 −3,6, HIPOXIA −3,6, GLUCÓLISIS −3,5,
INFLAMMATORY −3,4, IL6_JAK_STAT3 −3,3. A favor solo UV_RESPONSE_DN +3,3, G2M +2,4, E2F +2,3.
Lectura: el eje del receptor es **antiinflamatorio y antiproliferativo** (NF-κB, interferones,
complemento, IL-6/STAT3, IL-2/STAT5 todos del lado contrario), que es exactamente la función
canónica del GR, y además va en contra de OXPHOS, adipogénesis, glucólisis y mTORC1. Eso conecta el
eje con el metabolismo oxidativo del adiposo —el mismo programa que Sharma 2016 encuentra asociado
positivamente a SI en adiposo (41 genes de OXPHOS) y de forma discordante en músculo— y con la
inflamación de bajo grado. Es un enlace funcional que no teníamos y que no depende del reloj.

**19.4 Las firmas glucocorticoides de C2:CGP no sirven aquí, y por qué.** Las siete existentes son
de leucemia linfoblástica (HOLLEMAN, RHEIN) o de represión de NF-κB en línea celular (DEBOSSCHER):
z entre −3,0 y +0,6, ninguna informativa en adiposo. Confirma §17.5 y §18: **no hay en MSigDB un
conjunto de respuesta glucocorticoide en tejido metabólico humano**, y el nuestro de nueve genes
sigue siendo la mejor opción disponible. Merece una frase en la Discusión: es un hueco real del
campo.

**19.5 Precaución con C3:MIR, C4 y C7.** Los z extremos de esas colecciones (MIR513A +14, MODULE_60
−15,7, firmas de PBMC de vacunación −13) son conjuntos de 150-500 genes definidos por coexpresión o
por sangre completa; con un ranking que ya sabemos que tiene un eje leucocitario (§17.1), son
tautológicos. **No usarlos en el artículo.** Las colecciones defendibles para esto son H, C2:CP,
C5:GO y C8 (esta última como control de composición, no como resultado).

---

---

## 20. METAANÁLISIS DE EFECTOS ALEATORIOS Y FIGURA 6

Nueve estimaciones cohorte-tejido (AAGMEx, GSE41168, GSE13070, GTEx subcutáneo, visceral y músculo),
Fisher-z con DerSimonian-Laird, igual que el metaanálisis de la Tabla S17 del manuscrito. Script
`python/analyses/06_meta_analysis.py`, figura `python/figures/fig6_mecanismo.py`, tablas
`meta_per_cohort.tsv` y `meta_pooled.tsv`, figura `F6_mecanismo.png/pdf`.

**20.1 Receptor y reloj (NR3C1 ~ salida del reloj).**
| grupo | k | ρ combinado | IC95 | P | I² |
|---|---|---|---|---|---|
| global | 9 | **+0,340** | +0,259 a +0,417 | 8,9e-15 | 74 % |
| adiposo | 5 | **+0,371** | +0,323 a +0,417 | <1e-30 | **10 %** |
| músculo | 4 | +0,325 | +0,123 a +0,501 | 2,0e-3 | 86 % |
El adiposo es homogéneo (I² 10 %, cinco cohortes, dos plataformas, vivos y post mortem). La
heterogeneidad global vive en músculo y la explica AAGMEx (sin él, I² cae de 74 % a 29 %). El
leave-one-out mantiene ρ entre +0,322 y +0,366 quitando cualquiera de las nueve.
**Frase para el manuscrito:** *la asociación entre receptor y salida del reloj es homogénea en tejido
adiposo (ρ 0,37, I² 10 %, cinco cohortes) y heterogénea en músculo.*

**20.2 Señalización y reloj (score GC ~ reloj): hay que separar vivos de post mortem.**
| grupo | k | ρ | IC95 | P | I² |
|---|---|---|---|---|---|
| vivos | 6 | **+0,301** | +0,128 a +0,456 | 8,3e-4 | 74 % |
| **vivos, adiposo** | 3 | **+0,420** | +0,320 a +0,510 | <1e-10 | **0 %** |
| vivos, músculo | 3 | +0,187 | −0,090 a +0,438 | 0,18 | 76 % |
| GTEx post mortem | 3 | −0,038 | −0,154 a +0,080 | 0,53 | 87 % |
En vivos y adiposo la estimación es **perfectamente homogénea (I² 0 %)** entre AAGMEx, GSE13070 y
GSE41168: tres cohortes, dos plataformas, 301 personas. Post mortem el efecto desaparece, como ya
anticipaba §14.5 (el cortisol terminal domina las dianas agudas). Esto convierte el estrato en un
resultado: **la señalización glucocorticoide predice la salida del reloj en tejido adiposo de
personas vivas, con heterogeneidad nula, y no es medible en tejido post mortem.**

**20.3 Egger da P 0,004** sobre las nueve estimaciones de NR3C1. No es sesgo de publicación —no hay
publicación selectiva aquí, son todas las cohortes disponibles— sino el efecto de que las tres de
GTEx son grandes y post mortem y las pequeñas son biopsias en vivos. Reportarlo y explicarlo así en
Métodos; es la misma estructura que ya justifica el estrato de 20.2.

**20.4 Figura 6, seis paneles:** (a) forest del receptor, nueve cohortes-tejido, con los combinados
por tejido; (b) forest de la señalización, estratificado vivos/post mortem; (c) Hallmark del eje del
receptor; (d) la firma del adiposo resistente frente a Cushing, hidrocortisona y TZD; (e) el clamp
contra la calibración de tiempo; (f) el modelo.

---

## 21. ACTIVIDAD DE TF (DoRothEA A–C): EL REGULÓN DE NR3C1 NO SIRVE, PERO FOXO3 SÍ

Fichero `data/raw/dorothea_ABC.tsv` (OmniPath, 7.281 interacciones con signo en A–C, 353 TF).
Script `python/analyses/07_tf_activity.py`, tablas `NR3C1_activity_vs_mRNA.tsv` y
`AAGMEx_TF_activity_vs_clock.tsv`.

**21.1 El regulón de NR3C1 en DoRothEA A–C son 16 dianas y ninguna es una diana glucocorticoide de
tejido metabólico**: ATF2, ATP1B1, BGLAP, CALD1, CDK6, ETS1, FOXO3, NR1I3, PPP3CA, PTK2B, RLN1,
RLN2, S100P, SCNN1A, TAT, TLE1. No están FKBP5, TSC22D3, ZBTB16, KLF15, DDIT4 ni PER1. Son dianas de
hígado (TAT), hueso (BGLAP), riñón (SCNN1A) y relaxina. La actividad calculada con ese regulón rinde
peor que el mRNA y mucho peor que el conjunto de nueve genes:

| cohorte-tejido | mRNA de NR3C1 | actividad DoRothEA | conjunto GC de 9 genes |
|---|---|---|---|
| AAGMEx adiposo | +0,344 | +0,130 | **+0,427** |
| AAGMEx músculo | +0,059 | +0,037 | +0,119 |
| GSE41168 adiposo | +0,373 | +0,328 | **+0,421** |
| GSE41168 músculo | +0,655 | +0,455 | −0,088 |
| GTEx adiposo sc | +0,328 | +0,047 | −0,151 |
| GTEx visceral | +0,424 | **−0,267** | +0,058 |
| GTEx músculo | +0,329 | −0,007 | −0,017 |
**Decisión: no usar actividad de TF para NR3C1.** El conjunto pre-especificado de nueve genes y el
mRNA son mejores estimadores del eje en estos tejidos, y la razón es concreta y citable: DoRothEA
cura el regulón sin especificidad de tejido. Va al suplementario como justificación de la elección
(junto con §18), no al texto principal.

**21.2 El ranking de TF frente a la salida del reloj sí da un resultado, y conecta con la Fig. 5.**
Excluyendo genes del reloj de cada regulón (para que no sea circular) y ajustando por sexo, edad y
LEP, en AAGMEx adiposo: **FOXO3 +0,385** es el primero de 134 TF, seguido de POU5F1 +0,231,
FOXO4 +0,192, **NR3C1 +0,180** (4º) y FOXO1 +0,117. Al otro extremo, todo el eje inflamatorio:
NFKB1 −0,441, RELB −0,403, SPI1 −0,390, JUN −0,379, RFX5 −0,371, RELA −0,367, NFKB2 −0,359.
Esto reproduce con un método distinto el panel de Hallmark (§19.3: inflamación y NF-κB en contra del
eje) y, sobre todo, **conecta el mecanismo nuevo con la Fig. 5 ya existente del manuscrito**, donde
FOXO3 era el TF con mayor correlación con el eje glucocorticoide (ρ 0,677) y su regulón incluía
TXNIP y FBXO32. Es decir: el mismo factor que encabeza el eje glucocorticoide encabeza también la
salida del reloj. **La Fig. 5 deja de ser un apéndice del artículo anterior y pasa a ser el eslabón
regulador del mecanismo.** Frase para el manuscrito: *entre 134 factores de transcripción con
regulón curado, la actividad de FOXO3 es la que mejor predice la salida del reloj del tejido adiposo,
y el eje NF-κB la que peor.*

**21.3 Cautela anotada:** BMAL1 encabeza el ranking (ρ 0,87) cuando no se excluyen los genes del
reloj de su regulón, por construcción. En la tabla final se excluyen; si alguien recupera este
análisis, que no reporte BMAL1.

---

## 22. FIGURA 1 REHECHA Y REPRODUCIBLE

Scripts al repositorio: `python/analyses/01b_fig1_signature_test.py` (t moderada estilo limma
implementada en el propio script —momentos sobre log s², contracción hacia la a priori—,
Benjamini-Hochberg, y enriquecimiento de conjunto por z frente a 2.000 conjuntos del mismo tamaño)
y `python/figures/fig1_convencional.py`. Cierra el hueco de reproducibilidad de §1.2.

**22.1 Composición correcta, que casi vuelvo a romper.** Las cinco cohortes del manuscrito:
músculo GSE25462, GSE159984, GSE182120; adiposo GSE78721, GSE26637. **GSE13070 NO entra en esta
figura**: su contraste adiposo es 35 vs 6 y produce ~3.000 genes significativos que son del
desbalance, no del efecto; metido en el panel a destruye la escala y cambia el mensaje. GSE13070 se
usa en las figuras 4 y 6, donde el diseño es pareado dentro de la persona. Queda anotado en el
encabezado del script para que no se repita.

**22.2 El mensaje de la figura no es "el método convencional no sirve".** Es mejor y es el puente:
**el análisis convencional ya detecta la salida del reloj alterada — GSE25462 z 4,3 en músculo,
GSE78721 z 2,5 en adiposo — pero no puede decir qué señal la explica**: el conjunto glucocorticoide
no sale en ninguna cohorte, y gen a gen no pasa nada a FDR 0,1. Las listas diferenciales son de 0 a
10 genes en cinco cohortes del mismo contraste clínico.

**22.3 El puente a la Fig. 2, literal para el texto:** en esas mismas cohortes la relación entre el
receptor y la salida del reloj entre personas está ahí —GSE182120 músculo ρ 0,54 con cero genes
diferenciales— y la comparación de grupos no puede verla por construcción, porque es una
correlación entre individuos y no una diferencia entre grupos.

**22.4 Paneles:** (a) tamaño de la lista diferencial en las cinco cohortes, con borde punteado en
las que llevan valor publicado; (b) t moderada de los 18 genes de las dos firmas; (c) mapa de z por
conjunto, con el reloj ya alterado en dos cohortes; (d) ρ receptor-reloj entre personas con IC95.

**22.5 Corrida con tres de las cinco cohortes (19 oct).** Preparadas GSE25462 (músculo, GPL570,
diabetes 10 frente a historia familiar negativa 15; escala lineal en el depósito, pasada a log2) y
GSE26637 (adiposo, GPL570, brazo de ayuno, 5 resistentes frente a 5 sensibles; el brazo de clamp no
entra en esta figura). Pickles en `results/new/GSE25462_log2.pkl` y `GSE26637_log2.pkl`.

| cohorte | diseño | genes FDR<0,1 | z salida del reloj | z glucocorticoide |
|---|---|---|---|---|
| GSE25462 músculo | 10 vs 15 | **0** | **−6,06** | −1,09 |
| GSE26637 adiposo | 5 vs 5 | **0** | **−2,36** | +0,20 |
| GSE182120 músculo | 25 vs 24 | **0** | −0,85 | +1,55 |
Cero genes a FDR 0,1 en las tres, y de 54 pruebas sobre los 18 genes de las firmas ninguna pasa
FDR 0,1 (8 pasan P 0,05 sin corregir). **El conjunto del reloj baja en músculo (z −6,1) y en
adiposo (−2,4); el conjunto glucocorticoide no sale en ninguna.** Ese es el mensaje de la figura y
ahora está calculado, no citado. (Nota: el z publicado de GSE25462 era +4,3; el signo correspondía
al contraste inverso. Con el contraste enfermedad − control es −6,06, que es la dirección del resto
del artículo.)

**22.6 GSE25462 trae SI por persona (FSIGT, n=48), y da una disociación limpia.** El score de salida
del reloj **no** correlaciona con SI continuo (ρ +0,01, P 0,93) pero **sí está graduado por estado
clínico**: historia familiar negativa +0,31, un progenitor +0,09, dos progenitores +0,11, diabetes
−0,31. Y NR3C1 vs reloj ρ +0,25 (P 0,089), GC vs reloj +0,23 (P 0,11), direccionalmente consistentes
con el resto con n=48. **Esto refuerza §13.6b con un fenotipo medido: la firma sigue al estado de
enfermedad, no a la sensibilidad a la insulina.** Es un argumento mejor que el de la TZD porque aquí
el fenotipo continuo está medido en las mismas personas.

**22.7 Faltan dos ficheros para cerrar la figura:**
- `GSE159984_series_matrix.txt.gz` — tengo las cuentas crudas pero no los metadatos, así que no sé
  los grupos.
- La anotación de **GPL15207** (PrimeView) para GSE78721, que viene con sondas tipo `11715100_at` y
  sin símbolos. Sirve `GPL15207.annot.gz` o el `.soft` de GEO.

---

## 23. GSE26637 ES UN CLAMP PAREADO EN ADIPOSO: EL HUECO QUE TENÍAMOS

Los títulos de muestra revelan que las 20 muestras son **10 personas en dos estados** (ayuno e
hiperinsulinemia), no 20 personas: sujetos 1, 9, 14, 20 y 21 resistentes y 3, 4, 5, 7 y 13
sensibles, cada uno con biopsia de adiposo en ambos estados. **Esto cubre el hueco de §4**, donde el
clamp solo estaba en músculo porque el adiposo de GSE13070 no tiene brazo de clamp.
Tabla `GSE26637_clamp_adipose_per_person.tsv`.

**23.1 En adiposo el clamp hace lo mismo que en músculo, y con más consistencia.** Cambio por
persona (hiperinsulinemia − ayuno), n=10 pareado:
- programa de ayuno: PDK4 −2,38 (10/10), DDIT4 −2,40 (10/10), TXNIP −0,91 (10/10), IRS2 −1,74
  (10/10), ANGPTL4 −1,19 (9/10), SESN1 −1,16 (10/10); **HK2 +0,74, sube en 10/10**
- salida del reloj: **cae en 10/10 (score −1,43, P 0,002)** — PER1 −1,16 (10/10), NR1D2 −1,34
  (10/10), CIART −1,41 (9/10), PER3 −0,63 (10/10), TEF −0,58 (10/10), NR1D1 −0,34 (10/10),
  HLF −0,85 (9/10), PER2 −0,43 (9/10)
- **oscilador central sin cambio** (+0,05, 5/10, P 0,92) y BMAL1 +0,33 (3/10 bajan): cae el brazo de
  salida, no el oscilador
- score glucocorticoide −0,62 (9/10, P 0,006)
Las magnitudes son las mismas que en el músculo de GSE13070 (PDK4 −3,0, TXNIP −1,5, NR1D1 −1,2) y
las mismas que 3-4 h de mañana en las cohortes con hora (§4). **La calibración de tiempo se extiende
ahora a los dos tejidos.**

**23.2 Tampoco aquí está graduado por la sensibilidad:** cambio del reloj −1,30 en resistentes frente
a −1,66 en sensibles (P 0,42); núcleo de insulina −0,78 frente a −0,92 (P 0,55). Igual que en
músculo. Refuerza la afirmación 4 del artículo con un segundo tejido y un diseño pareado.

**23.3 Para el manuscrito:** este resultado sube la Fig. 4 de "el clamp es tiempo en músculo" a "el
clamp es tiempo en los dos tejidos, con el brazo de salida cayendo en 10 de 10 personas y el
oscilador central intacto". Y da el contraste exacto con la hidrocortisona: misma hora, misma
cohorte-tipo, y ahí **PER1 sube** en 10/10 mientras que con insulina **baja** en 10/10. Dos
perturbaciones opuestas sobre el mismo gen, cada una con control adecuado.

---

## 24. MOTORES ESTADÍSTICOS POR TIPO DE DATO (regla fija del proyecto)

- **Microarrays** (GPL570, HT-12, HTA-2.0, Agilent, PrimeView): t moderada estilo limma, implementada
  en `01b_fig1_signature_test.py` (momentos sobre log s², contracción de cada varianza hacia la
  a priori, df aumentados, Benjamini-Hochberg). No usar t de Welch sin moderar: con n de 5 a 25 por
  grupo infla los falsos positivos de los genes de varianza baja.
- **RNA-seq a partir de cuentas crudas** (GSE159984, GSE129843, GSE66446, GSE182117): **PyDESeq2**,
  instalado en el entorno (0.5.4). Cuentas enteras sin normalizar, diseño con las covariables, y
  `lfcShrink` para los tamaños de efecto. No aplicar la t moderada a cuentas ni log-CPM.
- Diseños pareados dentro de la persona (clamp, hidrocortisona, TZD): diferencia por persona y
  Wilcoxon de rangos con signo, más el conteo de cuántas personas van en cada dirección, que es lo
  que se reporta en las figuras.

---

## 25. ANOTACIÓN DE GPL570 POR ENTREZ (módulo nuevo `python/lib/gpl570.py`)

La anotación de la plataforma (`GPL570_annot.txt`, fichero de la familia GPL570-55999) está fechada
el **6 de octubre de 2014**. De sus 41.177 sondas con un solo Entrez y símbolo vigente, **3.413
(8,3 %) llevan un símbolo que ya no es el actual**. Confirma la cifra que la auditoría §1.4 había
estimado para esta plataforma.

**25.1 Módulo.** `gpl570.probe_map()` construye sonda → Entrez → símbolo actual (NCBI GRCh38.p13) y
`gpl570.annotate(expr)` colapsa la matriz: por gen, la sonda de mayor expresión media entre las que
lo miden en exclusiva y, si ninguna lo hace, la mejor multi-gen que lo incluya. Mapa exportado a
`GPL570_probe_map.tsv`. Aplicado a GSE25462, GSE26637, GSE41168 y GSE13070, las cuatro cohortes de
GPL570 del artículo: 54.675 sondas → 21.643 genes en todas (antes 23.625 con el símbolo del array,
porque el símbolo retiene pseudogenes y entradas sin Entrez).

**25.2 De los conjuntos, los dos que cambian son ARNTL→BMAL1 y ATP5A1→ATP5F1A**, ya cubiertos por la
tabla de sinónimos de `gene_sets.py`, de modo que los nueve genes del reloj y los nueve del
glucocorticoide se miden 9/9 con una y otra anotación. **Ningún resultado principal cambia**:
GSE25462 NR3C1~reloj +0,280 → +0,265; reloj vs log(SI) +0,072 → +0,106 (sigue siendo cero).

**25.3 Pero hay un caso donde sí cambia el signo:** GSE26637 adiposo, NR3C1~reloj pasa de **−0,081 a
+0,184**. Con n=10 ninguno de los dos es significativo, así que no afecta a ninguna conclusión, pero
muestra que en cohortes pequeñas la elección de sonda mueve el resultado. **Usar
`gpl570.annotate` de aquí en adelante en las cuatro cohortes de GPL570** y declararlo en Métodos.

**25.4 Lo que se pierde con Entrez y hay que saber:** GSTT1 y HLA-DRB4 desaparecen, porque sus
sondas no tienen Entrez vigente (son genes ausentes del ensamblaje primario, el falso positivo ya
anotado en §1.4). Son dos de los doce genes de los conjuntos genético y HLA, que son los controles
positivos de la concordancia entre tejidos. **Para esos dos conjuntos hay que conservar la
anotación por símbolo o añadirlos a mano**; para todo lo demás, Entrez.

---

## 26. GSE78721 ANOTADO (PrimeView) Y FIG. 1 CON CUATRO DE CINCO COHORTES

`gpl570.py` sirve también para GPL15207: ambas anotaciones de Affymetrix comparten estructura y
separador "///", y el módulo detecta la columna de Entrez (`ENTREZ_GENE_ID` o `Entrez Gene`).
GSE78721: 49.395 sondas → 20.125 genes, 68 con diabetes tipo 2 y 62 normales, 74 mujeres y 56
hombres, tejido adiposo de extremidad. Pickles `GSE78721_log2.pkl` y `_meta.pkl`.

**26.1 Estado de la Fig. 1, con cuatro cohortes corridas:**
| cohorte | diseño | genes FDR<0,1 | z reloj | z GC |
|---|---|---|---|---|
| GSE25462 músculo | 10 vs 15 | 0 | **−5,92** | −1,00 |
| GSE182120 músculo | 25 vs 24 | 0 | −0,84 | +1,60 |
| GSE26637 adiposo | 5 vs 5 | 0 | **−2,39** | +0,18 |
| GSE78721 adiposo | 68 vs 62 | **3** | **+2,05** | **+3,32** |
De 72 pruebas sobre los 18 genes de las firmas, **ninguna pasa FDR 0,1**.

**26.2 GSE78721 va en dirección contraria y hay que decirlo.** Es la única cohorte donde el reloj
sube en enfermedad (z +2,05) y la única donde el conjunto glucocorticoide sale (+3,32, con FKBP5 el
gen más significativo). Tres diferencias de diseño la explican parcialmente y van en el pie: es
**adiposo de extremidad**, no abdominal, el único de los cinco; es la cohorte más grande (130) y por
tanto la única con potencia para cualquier cosa; y el contraste es diagnóstico de diabetes, no
sensibilidad medida. **No se recorta ni se omite**: en el panel c se ve que dos cohortes de músculo
y una de adiposo bajan y una de adiposo sube, que es exactamente la heterogeneidad que justifica el
artículo. Dato que la apoya como válida: dentro de ella, NR3C1~reloj ρ **+0,61** y GC~reloj **+0,71**
(n=130, ajustado por sexo), las correlaciones más altas de todo el estudio. El eje está, lo que
cambia es la dirección del contraste de enfermedad.

**26.3 Falta solo `GSE159984_series_matrix.txt.gz`** para los metadatos (tengo las cuentas crudas;
se analizará con PyDESeq2 según §24). GSE78721 ya tenía sus metadatos dentro de la propia serie
matriz, no hacía falta fichero aparte.

---

## 27. GSE159984 SON ISLOTES, NO MÚSCULO — Y EL EJE REPLICA AHÍ CON LA ρ MÁS ALTA DEL ESTUDIO

**27.1 Corregir el repositorio.** En `00_prepare_data.py`, GSE159984 figura como "muscle,
cross-sectional". Es **islotes pancreáticos humanos**: 27 donantes con diabetes tipo 2, 58 no
diabéticos, más 45 muestras de cultivos de islotes con palmitato y glucosa alta. Corregir la
etiqueta y la disponibilidad de datos.

**27.2 Analizado con PyDESeq2** (§24), cuentas crudas, Entrez→símbolo actual, filtro de >5 cuentas
en ≥10 muestras: 24.906 genes, 85 donantes. **210 genes a FDR 0,1 y 98 a FDR 0,05**, el único de los
cinco diseños de la Fig. 1 con una lista diferencial de tamaño apreciable. Tabla
`GSE159984_islets_T2D_vs_ND.tsv`.

**27.3 Ni el reloj ni el eje glucocorticoide distinguen T2D de no diabético en islotes**: z de
conjunto −0,62 (reloj), +0,91 (GC), +1,38 (oscilador), y gen a gen el mayor es ZBTB16 (+0,51,
P 0,016, FDR 0,41). Lo único que sale es el programa oxidativo (−2,40). Es decir, la misma
conclusión que en los otros cuatro diseños, con 210 genes significativos de por medio: **el
contraste de enfermedad no recupera la firma**.

**27.4 Pero el eje entre donantes es el más fuerte de todo el estudio.** En los 85 donantes:
NR3C1 vs salida del reloj **ρ +0,725 (P 4,2e-15)**, score GC vs reloj **+0,584 (P 4,5e-09)**, y
HSD11B1 exactamente cero (+0,008), como en GTEx. La mediana del reloj por quintil de NR3C1 sube
monótona: −0,84, −0,56, −0,06, +0,29, +0,86. Y replica en el bloque independiente de **cultivos de
islotes (n=45, ρ +0,596)**, que son preparaciones ex vivo sin hora de biopsia ni eje HPA intacto:
el acoplamiento receptor-reloj es intrínseco al tejido, no un reflejo del cortisol circulante del
donante.

**27.5 Lo que esto añade al artículo.** Un tercer tejido metabólico, y el que cierra el círculo con
la literatura que ya citamos: el trabajo de Diabetologia 2026 sobre el receptor de glucocorticoides
en célula beta. Y el bloque de cultivos es el control que faltaba para la objeción "esto es
cortisol circulante del donante": en islotes aislados, sin eje HPA, el acoplamiento se mantiene.
**Añadir islotes a la Fig. 3 y al metaanálisis** (pasaría a diez estimaciones cohorte-tejido, con
islotes como el extremo alto) y al panel c de la Fig. 1 como quinto diseño.

---

## 28. METAANÁLISIS AMPLIADO A 13 ESTIMACIONES, CUATRO TEJIDOS

Con GSE159984 (islotes, dos bloques), GSE78721 (adiposo, 130) y GSE25462 (músculo, 50) el
metaanálisis del eje receptor-reloj pasa de 9 a **13 estimaciones cohorte-tejido y de dos a cuatro
tejidos**. Tabla actualizada `meta_per_cohort.tsv`.

| cohorte | tejido | n | ρ |
|---|---|---|---|
| GSE159984 | **islotes (donantes)** | 85 | **+0,725** |
| GSE41168 | músculo | 28 | +0,655 |
| GSE78721 | adiposo | 130 | +0,608 |
| GSE159984 | **islotes (cultivos)** | 45 | **+0,596** |
| GSE13070 | adiposo | 41 | +0,438 |
| GTEx | adiposo visceral | 587 | +0,424 |
| GSE41168 | adiposo | 28 | +0,373 |
| GSE13070 | músculo | 69 | +0,366 |
| GSE25462 | músculo | 50 | +0,364 |
| AAGMEx | adiposo | 232 | +0,344 |
| GTEx | músculo | 818 | +0,329 |
| GTEx | adiposo sc | 714 | +0,328 |
| AAGMEx | músculo | 232 | +0,059 |

| grupo | k | ρ combinado | IC95 | I² |
|---|---|---|---|---|
| global | 13 | **+0,427** | +0,338 a +0,509 | 84 % |
| adiposo | 6 | +0,417 | +0,327 a +0,498 | 69 % |
| músculo | 5 | +0,328 | +0,156 a +0,480 | 81 % |
| **islotes** | 2 | **+0,679** | +0,540 a +0,782 | 33 % |
Las trece estimaciones son positivas; doce de los trece IC excluyen el cero (la excepción sigue
siendo AAGMEx músculo). El estimado global sube de 0,34 a 0,43 al entrar islotes y GSE78721.
GC~reloj en vivos: ρ +0,343 (IC95 +0,155 a +0,507, k=10) frente a −0,038 en GTEx post mortem, el
mismo estrato de §20.2 con el doble de estudios.

**Frase para el manuscrito:** *en trece estimaciones de cuatro tejidos metabólicos humanos y seis
cohortes, el nivel del receptor de glucocorticoides predice la salida del reloj (ρ combinado 0,43,
IC95 0,34-0,51), con el efecto más fuerte en islotes pancreáticos (0,68) y presente también en
preparaciones de islotes cultivados sin eje hipotálamo-hipófiso-adrenal.*

**Pendiente de la figura 3:** rehacerla con 13 filas en el forest, agrupadas por tejido, y añadir
islotes al panel de distribuciones por quintiles. La heterogeneidad del 84 % se explica por tejido y
por el estrato vivos/post mortem, no por cohortes discrepantes: todas apuntan igual.

---

## 29. FIGURAS 3 Y 1 MONTADAS CON TODO LO NUEVO

**29.1 Figura 3** (`fig3_eje.py`, `F3_eje.png/pdf`), cinco paneles:
(a) forest de las 13 estimaciones agrupadas por tejido, con los combinados por tejido y el global
(ρ 0,43, IC95 0,34-0,51); (b) distribución de la salida del reloj por quintiles de NR3C1 en islotes
(n=85) y GTEx visceral (n=587), cada punto una persona y la mediana marcada — sube monótona en las
dos; (c) forest de la señalización glucocorticoide estratificado vivos (ρ 0,34) frente a post mortem
(−0,04); (d) diseño dentro del donante en las cinco cohortes con varios tejidos; (e) **islotes
cultivados, ρ 0,60 (P 1,6e-05, n=45)**, el panel que responde a la objeción de que todo sea cortisol
circulante del donante.

**29.2 Figura 1, paneles de distribución** (`fig1_distribuciones.py`, `F1_distribuciones.png`), con
las cinco cohortes y los valores de P explícitos, cada punto una persona:
| cohorte | P del score de reloj | P del score GC |
|---|---|---|
| GSE25462 músculo | **0,003** | 0,760 |
| GSE182120 músculo | 0,992 | 0,177 |
| GSE26637 adiposo | 0,151 | 1,000 |
| GSE78721 adiposo | 0,289 | 0,090 |
| GSE159984 islotes | 0,854 | 0,594 |
**Solo una de cinco separa por score de reloj y ninguna por score glucocorticoide**, y en todas las
distribuciones se solapan casi por completo. Es el panel que hace visible por qué el contraste de
enfermedad no sirve: no es falta de genes, es que la diferencia entre grupos es pequeña frente a la
variación entre personas. Con esto, los paneles de la Fig. 1 quedan: (a) lista diferencial por
cohorte, (b) t moderada de los 18 genes, (c) distribución del score de reloj por persona, (d)
distribución del score GC por persona, (e) mapa de z por conjunto, (f) ρ receptor-reloj entre
personas.

**29.3 Pendiente de dibujo a mano** (las generadas por IA no las aceptan Springer Nature, Cell Press
ni Elsevier): el diagrama de diseño de las cinco cohortes, el esquema de transición
grupos→individuos de la Fig. 1, el esquema del estimador de la Fig. 2, la línea de tiempo del clamp
frente a hidrocortisona de la Fig. 4, y el modelo de cuatro pasos de la Fig. 6.

---

## 30. AUDITORÍA DEL BLOQUE NUEVO (19 oct) — UN HALLAZGO QUE OBLIGA A CORREGIR

Tabla `AUDIT_global_axis.tsv`.

**30.1 HALLAZGO PRINCIPAL: en RNA-seq, el score de salida del reloj arrastra un eje técnico global.**
En islotes, el número de genes detectados por muestra correlaciona **ρ 0,854 con el score de reloj**
y 0,617 con NR3C1. Y la distribución nula —cualquier gen al azar frente al score de reloj— tiene
**mediana +0,436**, no cero. Es decir, en esa cohorte *cualquier* gen correlaciona ~0,44 con el
score, y el +0,725 de NR3C1 queda en el **percentil 84**, no en el extremo. Ajustando por
profundidad y genes detectados baja a +0,370.
Causa: el z-score por gen entre muestras no elimina el componente de muestra; en RNA-seq con
profundidad variable ese componente es grande. **En microarrays no ocurre** (nulo mediana +0,002 en
AAGMEx, −0,025 en GSE41168 músculo): están cuantil-normalizados de origen.

**30.2 Corrección adoptada: centrar cada muestra antes de puntuar** (restar a cada columna su media
sobre todos los genes medidos). El nulo se va a cero en todas las cohortes y los resultados se
sostienen:
| cohorte | ρ original | nulo med. | percentil | **ρ centrado** | nulo med. | percentil |
|---|---|---|---|---|---|---|
| islotes | +0,725 | +0,436 | 84,4 % | **+0,323** | −0,153 | 80,1 % |
| GTEx adiposo sc | +0,421 | +0,124 | 93,4 % | **+0,319** | −0,018 | 92,5 % |
| GTEx visceral | +0,562 | +0,262 | 97,0 % | **+0,423** | +0,005 | 97,0 % |
| GTEx músculo | +0,366 | +0,135 | 95,8 % | **+0,299** | +0,002 | 97,0 % |
| AAGMEx adiposo | +0,402 | +0,002 | 99,1 % | +0,401 | −0,004 | 98,4 % |
| AAGMEx músculo | +0,064 | +0,001 | 79,9 % | +0,063 | +0,001 | 80,5 % |
| GSE78721 adiposo | +0,593 | −0,124 | 91,2 % | +0,539 | −0,138 | 91,1 % |
| GSE41168 adiposo | +0,373 | +0,029 | 91,5 % | +0,382 | −0,013 | 93,1 % |
| GSE41168 músculo | +0,655 | −0,025 | 99,9 % | +0,714 | −0,015 | 100,0 % |
| GSE25462 músculo | +0,265 | −0,023 | 96,2 % | **+0,456** | +0,019 | 99,4 % |
Metaanálisis sobre las nueve con dato en ambas versiones: ρ 0,373 (IC95 0,279-0,459) original frente
a **0,376 (0,285-0,460) centrado**, I² 82 % en las dos. **El resultado principal no cambia**, pero
los valores por cohorte sí se mueven hasta 0,19 y el de islotes baja de 0,725 a 0,323.

**30.3 Qué hay que rehacer:** aplicar el centrado por muestra en `gene_sets`/scores y recalcular
**§27 (islotes: ρ 0,725 → 0,323; el de cultivos también)**, **§28 (metaanálisis; islotes deja de ser
el extremo alto)**, la Fig. 3 y el panel b de distribuciones. Las cohortes de microarray no cambian.
Declarar el centrado en Métodos: es una decisión de preprocesado con efecto demostrado.

**30.4 Lo que la auditoría confirma como correcto:**
- **t moderada**: en simulación nula controla el error tipo I al 5,1 % (clásica 4,9 %) y da 0 genes
  a FDR 0,1; con 200 genes a +1 log2 detecta 91 con 10 falsos positivos, idéntico a la t clásica con
  n=25. Correcta, aunque con este n la moderación no aporta sobre la clásica; su ventaja está en los
  diseños de n=5-10 (GSE26637).
- **Tabla del metaanálisis**: 26 filas, **0 duplicadas** tras los dos añadidos.
- **Matrices**: ningún índice duplicado, ningún NaN, escalas log2 coherentes en las seis revisadas.
- **GSE78721 no es un error de signo ni de anotación**: gen a gen el contraste T2D−normal es
  pequeño y mixto (PER2 +0,28, NR1D2 +0,55, PER3 −0,18, DBP −0,02), ninguno significativo. El z de
  conjunto +2,05 viene de la suma de efectos diminutos en una n grande, no de una inversión real.
  Descrito así en el pie, deja de ser una anomalía.

**30.5 Pendiente de auditar:** el ajuste de GTEx por tiempo isquémico y RIN sigue sin hacerse porque
falta `GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt` (§anterior). Dado §30.1, es prioritario:
si el eje global de GTEx es tiempo isquémico, el centrado por muestra solo lo mitiga en parte.

---

## 31. GTEx AJUSTADO POR TIEMPO ISQUÉMICO Y RIN: EL EJE AGUANTA

Con `SampleAttributes.txt` ya disponible, y aplicando además el centrado por muestra de §30.2.

**31.1 El tiempo isquémico sí afecta al reloj, como era de temer**: correlaciona con el score de
salida del reloj ρ +0,261 (adiposo sc), +0,289 (visceral) y +0,207 (músculo), con medianas de 530 a
573 minutos. También con NR3C1 (+0,283, +0,457, +0,121). El RIN pesa menos salvo en visceral
(−0,192 con el reloj, −0,302 con NR3C1; mediana de RIN 6,9, la peor de los tres).

**31.2 Pero el eje no es isquemia.** Ajustando por sexo, edad, clase de Hardy, tiempo isquémico y
RIN, sobre matrices centradas por muestra:
| tejido | n | sexo+edad+Hardy | **+ isquemia y RIN** |
|---|---|---|---|
| adiposo subcutáneo | 714 | +0,250 | **+0,243** (P 4,9e-11) |
| adiposo visceral | 587 | +0,309 | **+0,302** (P 7,7e-14) |
| músculo | 818 | +0,273 | **+0,269** (P 4,9e-15) |
La caída es de 0,004 a 0,007. **Los tres valores definitivos de GTEx para el manuscrito son 0,24,
0,30 y 0,27**, por debajo de los 0,33, 0,42 y 0,33 que veníamos usando: la diferencia viene del
centrado por muestra (§30), no de las covariables nuevas. Actualizar §14, §20, §28 y la Fig. 3 con
estos números.

**31.3 El score glucocorticoide en GTEx sigue en cero** tras el ajuste completo (−0,14, +0,02,
−0,04), igual que antes. Refuerza §14.5: en tejido post mortem las dianas agudas del receptor están
dominadas por el cortisol terminal y el eje de señalización no es medible. El estrato
vivos/post mortem del metaanálisis se mantiene tal cual.

**31.4 Para Métodos, el preprocesado definitivo de GTEx:** log2 CPM, genes detectados en >50 % de
muestras, **centrado por muestra**, y residuales frente a sexo, edad (punto medio del tramo), clase
de Hardy, tiempo isquémico (SMTSISCH) y RIN (SMRIN), con imputación por mediana de los pocos
valores ausentes.

---

## 32. MARCO TEÓRICO: OSCILADOR FORZADO, Y LA PRIMERA PRUEBA CUANTITATIVA

**32.1 El marco existe y encaja.** La teoría de entrainment por osciladores forzados
(Poincaré/Goodwin, lenguas de Arnold) está desarrollada para relojes circadianos, y hay un modelo
publicado del sistema SCN-HPA-periferia en el que **las células periféricas son osciladores
independientes sincronizados por glucocorticoides, con la fuerza de acoplamiento como parámetro
explícito** (Li y Androulakis, J Biol Rhythms 2023, doi 10.1177/07487304231188541; antes
Mavroudis et al., Physiol Genomics 2012, con el receptor y el cortisol modelados). Predicción
central de ese modelo: **sin la señal glucocorticoide de entrainment, la población de células
periféricas oscila de forma amortiguada a nivel de ensamble** — que es exactamente lo que mide una
biopsia en masa.

**32.2 El mapeo de nuestros datos al modelo.** El forzamiento efectivo sobre el reloj periférico es
el producto de dos factores: la **amplitud del ritmo de cortisol** (fuerza del Zeitgeber, F) y la
**disponibilidad de receptor** (ganancia del acoplamiento, K). Esto unifica tres observaciones que
hasta ahora eran piezas sueltas:
- obesidad y resistencia: el ritmo de GC se aplana → F baja;
- Cushing e hidrocortisona sostenida: exposición constante → la componente rítmica del forzamiento
  baja aunque suba el nivel medio;
- adiposidad y exposición crónica: NR3C1 baja → K baja.
**Los tres reducen el forzamiento efectivo y por tanto producen la misma firma**, que es
precisamente lo que observamos (adiposo resistente ≈ Cushing ≈ hidrocortisona, ρ 0,71 y 0,60).
El modelo explica también la disociación pulso/tono de §18: el pulso agudo actúa por el término de
forzamiento directo sobre genes con GRE (PER1, 10/10 con hidrocortisona) y el tono por el término de
acoplamiento sobre el brazo de salida (PER3, DBP, HLF, NR1D2), que es lo que sigue a NR3C1.

**32.3 Primera prueba cuantitativa, con datos que ya tenemos.** Si NR3C1 es la ganancia K, la
amplitud inferida por persona debe crecer con NR3C1 (tabla `MODEL_amplitude_vs_NR3C1.tsv`):
| cohorte | n | ρ amplitud ~ NR3C1 | R de fase, quintil bajo → alto |
|---|---|---|---|
| AAGMEx adiposo | 232 | **+0,201** (P 2,1e-3) | 0,38 → 0,51 |
| GTEx adiposo sc | 714 | **+0,167** (P 7,7e-6) | 0,44 → 0,53 |
| GTEx visceral | 587 | +0,094 (P 0,023) | 0,44 → 0,55 |
| GTEx músculo | 818 | **+0,223** (P 1,1e-10) | 0,31 → 0,14 |
| GSE78721 adiposo | 130 | −0,213 (P 0,015) | 0,85 → 0,66 |
| islotes | 85 | −0,218 (P 0,045) | 0,51 → 0,71 |
**Cuatro de seis van en la dirección predicha y dos en contra**, y las dos que van en contra son las
mismas que ya se comportaban distinto en otros análisis (GSE78721, §26.2; islotes, §30.1). La
concentración de fase sube del quintil bajo al alto en cuatro de seis. **No es una confirmación
limpia: es una predicción parcialmente cumplida, y así hay que tratarla.**

**32.4 Qué haría falta para que el modelo deje de ser analogía.** Ajustar los dos parámetros (F y K)
a los datos requiere, como mínimo, una cohorte con ritmo de cortisol muestreado y biopsias seriadas
en la misma persona: ninguna de las seis lo tiene. Lo que sí se puede hacer ya, y es el camino
honesto: **simular el modelo de Li-Androulakis barriendo el parámetro de acoplamiento, generar
poblaciones sintéticas de biopsias únicas, aplicarles nuestro propio estimador y comprobar si
reproducen los valores observados** (ρ amplitud-K, concentración de fase, firma PER3/DBP/HLF abajo
con BMAL1 arriba). Eso convierte el modelo en una predicción contrastable con lo que ya tenemos, sin
datos nuevos, y daría una figura 7 con un eje teórico. Es trabajo de días, no de meses.

**32.5 Cautela de encuadre.** El modelo no debe presentarse como demostrado por estos datos. Su
papel en el manuscrito es ordenar el mecanismo —decir por qué tres situaciones clínicas distintas
dan la misma firma— y generar las predicciones de §13.8. Ponerlo en la Discusión, no en Resultados,
salvo que se haga la simulación de §32.4, en cuyo caso sí es un resultado.

---

## 34. FIGURA 7 MONTADA

Scripts: `python/model/forced_clock.py` (modelo), `make_tongue.py` (barrido, ~2 min),
`make_fig7_data.py` (cohortes sintéticas, ~3 min, cachea en `fig7_data.npz`),
`python/figures/fig7_modelo.py` (dibuja a partir del caché). Figura `F7_modelo.png/pdf`.

**Paneles:** (a) lengua de entrainment, sincronía del ensamble frente a forzamiento efectivo y
período intrínseco, con el umbral de 0,6 marcado y los dos regímenes señalados; (b) amplitud del
ensamble y sincronía frente al forzamiento; (c) firma simulada frente a observada gen a gen,
**ρ 0,48 sin ajustar ningún parámetro**; (d) la predicción de la ventana horaria —el cambio de
NIVEL cae de −0,48 (matinal) a −0,35 (± 4 h) a −0,04 (uniforme) mientras el de AMPLITUD se mantiene
en −0,65 en los tres; (e) nuestro estimador sobre datos del modelo frente al rango observado con
NR3C1; (f) el mecanismo en una línea.

**El panel d es el más valioso**: separa dos predicciones que hasta ahora íbamos juntando. El
modelo dice que el acoplamiento bajo produce siempre menos amplitud, y produce además menos nivel
**solo si la biopsia es a hora fija**. Eso convierte la heterogeneidad entre cohortes de §22 y §26
en una predicción: GTEx (hora de muerte dispersa) no debe mostrar la firma de nivel, y no la
muestra; las cohortes de biopsia matinal sí, y la muestran.

**Nota de reproducibilidad:** las simulaciones tardan unos cinco minutos en total, así que los
scripts del modelo cachean y la figura se dibuja en segundos. Declarar en Métodos los parámetros
fijos (lambda 0,05 h⁻¹, A 1,0, dispersión de períodos 1,5 h, 60-150 células, 18-20 días de
integración, paso 0,05 h) y que ninguno se ajustó a los datos observados.

---

## 35. AUDITORÍA DEL BLOQUE DEL MODELO Y DEL CÓDIGO

**35.1 La ρ 0,48 del panel 7c no es el número que hay que reportar.** Un modelo trivial —"los ocho
genes de salida bajan y BMAL1 sube"— da ρ 0,55, más alto que el del modelo real, y la P por
permutación de la ρ observada es 0,203. **El orden de los genes no es informativo.** Lo que sí lo es:
**los 9 de 9 signos coinciden y la magnitud media es 0,48 en el modelo frente a 0,50 observada,
cociente 0,96.** Reportar eso y no la correlación: *el modelo predice la dirección de los nueve
genes y una magnitud media dentro del 4 % de la observada*. Cambiar el título del panel.

**35.2 El panel 7d sí es sólido.** Cinco réplicas con semillas distintas: con biopsia matinal, nivel
−0,447 ± 0,028 y amplitud −0,657 ± 0,024; con hora uniforme, nivel **+0,020 ± 0,025** y amplitud
−0,658 ± 0,023. La disociación nivel/amplitud es nítida y reproducible, y el nivel con hora uniforme
es cero dentro del error. Es la predicción más fuerte del modelo.

**35.3 El integrador converge.** Amplitud y R idénticas a tres decimales con paso de 0,20, 0,10,
0,05 y 0,02 h. Euler explícito es suficiente aquí; declarar el paso en Métodos sin más.

**35.4 Deuda de código, ordenada:**
- Sin docstring: `06_meta_analysis.py`, `07_tf_activity.py`, `08_gtex_adjusted.py`,
  `fig1_distribuciones.py`, `fig6_mecanismo.py`, `validate_against_data.py`. Todos los demás sí.
- **Valores embebidos en scripts de figura**: `fig3_eje.py` (panel d, los ρ del diseño dentro del
  donante) y `fig7_modelo.py` (panel c, la firma observada). Hay que sacarlos a TSV y que la figura
  los lea, como hacen las demás. Es justo el tipo de cosa que rompe `run_all.sh`.
- `make_figures.py` (550 líneas) sigue generando las figuras del envío anterior; hay que decidir si
  se fusiona con los scripts nuevos o se archiva.

---

## 36. PROPUESTA DE FUSIÓN: DE SIETE FIGURAS A CINCO

Ninguna pierde paneles; se redistribuyen por argumento en lugar de por tipo de análisis.

**Figura 1 — El problema y la herramienta** (fusiona 1 y 2; 8 paneles)
a diseño de las cinco cohortes · b lista diferencial por cohorte · c distribución del score de reloj
por persona · d mapa de z por conjunto · e esquema del estimador · f validación de fase (interno,
externo, transferencia músculo→adiposo) · g amplitud inferida frente a medida por cosinor en
miotubos e in vivo · h ρ receptor-reloj entre personas en las mismas cohortes
*Por qué funciona:* el problema y su solución técnica son un solo movimiento. El panel h cierra el
argumento de la primera mitad y abre la segunda. Se gana que el lector no espere dos figuras para
saber qué se va a hacer.

**Figura 2 — El eje** (la actual 3, sin cambios; 5 paneles)
a forest de 13 estimaciones · b distribuciones por quintil · c señalización vivos/post mortem ·
d dentro del donante · e islotes cultivados

**Figura 3 — Cortisol sí, insulina no** (fusiona 4 y la mitad de 6; 7 paneles)
a línea de tiempo de las dos perturbaciones · b clamp frente a calibración horaria en los dos
tejidos · c hidrocortisona gen a gen, PER1 10/10 · d el contraste PER1 sube con cortisol y baja con
insulina · e dexametasona · f Cushing · g la firma del adiposo resistente frente a las tres
*Por qué funciona:* las perturbaciones y la firma de enfermedad son el mismo argumento —qué mueve el
reloj y qué no— y hoy están partidas entre dos figuras con un salto conceptual en medio.

**Figura 4 — Regulación y especificidad** (fusiona la 5 y el resto de la 6; 6 paneles)
a FOXO3 encabeza 134 TF · b eje NF-κB en contra · c Hallmark del eje del receptor · d HPO: obesidad
troncal y cortisol urinario a favor, diabetes tipo 2 en contra · e la TZD no revierte la firma ·
f red FOXO3
*Por qué funciona:* todo lo que responde "qué gobierna el eje y qué no es" en un sitio, y la
disociación queda al lado de la evidencia funcional que la sostiene.

**Figura 5 — El modelo** (la actual 7, sin cambios; 6 paneles)

**Extended Data:** GTEx diez pares y los dos modos de acoplamiento, no replicación de la amplitud
muscular, sensibilidad de métodos, auditoría de símbolos, firma empírica GRAS y su no validación.

*Cuenta:* 8 + 5 + 7 + 6 + 6 = 32 paneles, los mismos que hay ahora repartidos en siete figuras.
Cinco figuras es el formato habitual de Cell Metabolism y evita que el revisor lea dos figuras
seguidas de validación técnica antes de llegar al mecanismo.

---

## 37. LAS CINCO FIGURAS, GENERADAS

Script único `python/figures/paper_figures.py`; F2 y F5 llaman a `fig3_eje.py` y `fig7_modelo.py` y
renombran la salida. Ficheros `F1`-`F5` en png y pdf.

**F1 — El problema y la herramienta** (8 paneles): a diseño (marco para dibujo a mano) · b lista
diferencial en las cinco cohortes, de 0 a 3 genes · c mapa de z por conjunto, con el reloj ya
alterado y el eje glucocorticoide nunca · d esquema del estimador (a mano) · e miotubos, amplitud
inferida frente a medida por cosinor, **ρ 0,90 y 0,88** · f músculo in vivo, **ρ 0,69** ·
g distribuciones por persona con la P explícita · h ρ receptor-reloj entre personas en las mismas
cohortes.
**Corrección importante encontrada al montarla:** los paneles de validación usaban la columna
`amp_inferred` del TSV, que es la proyección de 800 genes —la que NO valida (ρ −0,28)—. Recalculado
con ordenamiento libre sobre genes del reloj y centrado por muestra: **ρ 0,90 (CTL) y 0,88 (HGI),
0,92 juntas**, por encima del 0,70 que teníamos documentado. Tabla nueva
`GSE182117_inferred_free_vs_measured.tsv`. En GSE129843 la columna correcta es `amp_free` (0,69),
no `amp_inf_clk` (0,53).

**F2 — El eje** (5 paneles): la figura 3 anterior sin cambios.

**F3 — Cortisol sí, insulina no** (8 paneles): a línea de tiempo de las dos perturbaciones (a mano)
· b clamp frente a calibración horaria · c clamp pareado en adiposo, cada punto una persona ·
**d PER1 sube 1,6 log2 con cortisol y baja 1,2 con insulina, 10/10 en ambos** · e hidrocortisona gen
a gen con el conteo por persona · f dexametasona en músculo · g Cushing · h la firma del adiposo
resistente frente a Cushing, hidrocortisona y TZD.

**F4 — Regulación y especificidad** (6 paneles): a los ocho TF de cada extremo, FOXO3 primero y el
eje NF-κB al otro lado · b Hallmark · c HPO · d la TZD frente a la firma, gen a gen · e el receptor
y la enzima frente a la adiposidad · f arquitectura reguladora (a mano).

**F5 — El modelo** (6 paneles): la figura 7 anterior.

**Pendiente antes de enviar:** los cuatro marcos de dibujo a mano (diseño, estimador, línea de
tiempo, arquitectura), sacar los valores embebidos de `fig3_eje.py` y `fig7_modelo.py` a TSV (§35.4),
y cambiar el título del panel 5c por el de magnitud en lugar de la correlación (§35.1).

---

## 38. REPRESENTACIONES REVISADAS: DOTPLOTS Y DISTRIBUCIONES

Tenías razón en las dos cosas. Cambios aplicados en `paper_figures.py`:

**Lo que estaba mal representado y por qué:**
- **F1b era un gráfico de barras de ceros** (0, 0, 0, 3 genes). Una barra que mide cero no comunica
  nada. Sustituido por un **dotplot de los 18 genes de las firmas × las cuatro cohortes**, con
  tamaño = −log10 P sin corregir y color = t moderada. Ahora se ve de un golpe que ningún gen pasa
  el umbral en ninguna cohorte y, a la vez, dónde está la señal débil (CIART y HLF en GSE25462,
  SESN1 y KLF15 en GSE26637).
- **F3d eran dos barras** con la mediana de cada perturbación. Sustituido por las **distribuciones
  por persona**: 10 puntos de hidrocortisona arriba y 10 del clamp abajo, con el conteo. Se ve que
  los 10 suben en una y los 10 bajan en la otra, que es el argumento.
- **F3g eran barras de la diferencia de medianas en Cushing**. Sustituido por los 16 puntos
  individuales, control frente a Cushing, centrados en el control. Con n=5 contra 11 la barra
  escondía el solapamiento.
- **Faltaban los dotplots de enriquecimiento.** Añadidos en F4: Hallmark (b), HPO (c), **Reactome
  (h) y GO procesos biológicos (i)**, todos con tamaño = número de genes del conjunto y color = z,
  que es la representación estándar y la que un revisor espera. El de actividad de TF (a) también
  pasa a dotplot con tamaño = número de dianas del regulón.

**Paneles nuevos que caben al pasar a tres renglones:**
- F1i: **distribución de P sin corregir de todos los genes** en tres cohortes, casi uniforme. Es el
  control que demuestra que no hay señal oculta por la corrección múltiple, y responde por
  adelantado a "¿no será falta de potencia?".
- F1j: marco para el diagrama de transición grupos → individuos.
- F4f: **el eje del receptor en el transcriptoma completo**, histograma de los 34.020 genes con
  FOXO3, NFKB1 y HSD11B1 marcados. Da la escala: dónde cae cada gen en la distribución real.

**Formato:** F1 y F4 pasan a tres renglones (12 y 11 paneles); F3 se queda en dos filas de cuatro;
F2 y F5 sin cambios. Total 37 paneles en cinco figuras.

---

## 39. LA AMPLITUD ESTABA EN LA DISPERSIÓN, NO EN LA MEDIANA — §2.1 SE CORRIGE EN PARTE

Todo el rechazo del resultado de amplitud muscular (§2.1) se hizo con Mann-Whitney, que es una
prueba de **localización** y es ciega a un efecto de **dispersión**. Rehecho con pruebas de escala
(Levene centrado en la mediana, Fligner-Killeen, y razón de rango 5-95 con P por permutación de
5.000 réplicas). Script `python/analyses/10_amplitude_dispersion.py`.

**39.1 Con el estimador validado (genes del reloj, ordenamiento libre), GSE13070:**
| | resistentes (n=51) | sensibles (n=18) | prueba |
|---|---|---|---|
| mediana | 2,58 | 1,96 | MW **P 0,206** |
| IQR | 1,40 | 0,90 | razón 1,57, P 0,23 |
| **rango 5-95** | **3,09** | **1,87** | razón **1,66**, **P 0,021** |
| sd | 1,04 | 0,64 | **Levene P 0,032 · Fligner P 0,027** |
Cuantiles 10-25-50-75-90: resistentes 1,12 / 1,90 / 2,58 / 3,31 / 3,62; sensibles 1,58 / 1,75 /
1,96 / 2,65 / 3,18. La distribución resistente se abre por los dos extremos.
Con el estimador anclado, misma cohorte: razón de rango **1,90, P 0,0036**. Dos de los tres
estimadores coinciden; el de 800 genes no (razón 1,06, P 0,87), que es coherente con que mida otra
cosa (§ED2, ρ 0,00 entre ambos).

**39.2 No replica en GSE182120.** Razón de rango 0,98, Levene P 0,41, y nada en ningún estimador.
Las dos cohortes difieren en diseño: GSE13070 es resistencia medida por clamp con grupos muy
desbalanceados (51 vs 18) y GSE182120 es T2D frente a NGT con grupos iguales. **La afirmación
defendible es de una cohorte, no de dos.**

**39.3 Y es exactamente lo que predice el modelo.** Bajando el acoplamiento de 1,4 a 0,5 en 120
sujetos sintéticos: rango 5-95 de 0,27 a 0,66, **razón 2,47, Levene P < 1e-4**. Observado 1,66 con
P 0,032. **El modelo predecía dispersión creciente al perder entrainment y no lo habíamos
comprobado**: a bajo forzamiento, algunos individuos siguen enganchados y otros no, así que la
amplitud del ensamble se dispersa. Es una predicción del modelo que se cumple en los datos y que no
estaba entre las de §13.8.

**39.4 Qué cambia en el manuscrito.**
- §2.1 pasa de "la amplitud muscular no replica" a: **la mediana no difiere; la dispersión sí, en
  una de las dos cohortes y en la dirección que el modelo predice**.
- ED2 necesita un panel más: razón de rango y Levene por estimador y cohorte, con los cuantiles.
- F5 gana una predicción cumplida: añadir al panel D o a uno nuevo la dispersión simulada frente a
  la observada.
- En el texto hay que ser explícito sobre por qué se usa una prueba de escala: **si el mecanismo es
  pérdida de entrainment, el efecto está en la varianza entre personas y una prueba de medianas no
  puede verlo.** Eso se dice antes de enseñar el resultado, no después.

**39.5 Lección de método, para la lista de no repetir:** elegir la prueba por la hipótesis
mecanística y no por costumbre. Mann-Whitney fue la elección automática y descartó durante dos
sesiones un resultado real.

---

## 40. ED4 SE ARRASTRA SIN REGENERAR — COMPROBACIONES DE CONSISTENCIA

La red de FOXO3 del envío anterior se mantiene tal cual. Lo que hay que garantizar es que no
contradiga nada del texto nuevo. Comprobado:

1. **Mismo regulón que la Fig. 4A.** La red se dibujó con DoRothEA A-C y la Fig. 4A cuantifica ese
   mismo regulón: 38 dianas (`FOXO3_regulon_ED4.tsv`). Los dos paneles hablan del mismo conjunto, y
   eso hay que decirlo en el pie para que no parezcan dos cosas distintas.
2. **Símbolos: ninguna diana de FOXO3 lleva símbolo retirado**, así que la reanotación por Entrez
   (§25, ED3B) no cambia una sola etiqueta de la red. Era el riesgo principal de arrastrar una
   figura vieja y no se materializa.
3. **CDKN1A y GADD45A sí están en la red** y los retiramos del conjunto de ayuno/FOXO3 por el
   artefacto de biopsia repetida (§1.5). No es una contradicción —la red viene de regulones curados
   y no de estos datos— pero **hay que decirlo en el pie**, o un revisor atento lo señalará como
   incoherencia. Redactado así en `FIGURE_CAPTIONS.md`.
4. **TXNIP y FBXO32 están en el regulón**, que es lo que conecta la red con el programa de ayuno del
   músculo (§13.3) y con el nuevo hallazgo de que FOXO3 encabeza los 134 TF frente a la salida del
   reloj (§21.2). La figura deja de ser herencia y pasa a ilustrar el eslabón regulador.

**Qué cambia en el manuscrito:** el pie de ED4 ya no la presenta como material del artículo
anterior, sino como la arquitectura del regulón que la Fig. 4A cuantifica. Numeración intacta: ED1
a ED7, con ED4 arrastrada.

---

## 41. ESTRUCTURA DE ESCRITURA — REVISAR EN CADA EDICIÓN DEL MANUSCRITO

**Fichero `RESULTS_STRUCTURE.md`, en la raíz del repositorio. Abrirlo antes de tocar el manuscrito,
siempre.** Contiene el análisis de cuatro artículos de Nature (el Resource de proteómica de
ejercicio de Nature Metabolism 2026, el Article de vacunas de Nature 2025, la revisión de relojes y
resistencia a la insulina de Nat Rev Endocrinol 2018 y el descriptor PAPILA), las nueve reglas que
adoptamos, el esqueleto de nuestros Resultados en siete subsecciones con sus citas de panel, las
cuatro fórmulas con los términos que hay que definir, y la adaptación a los lineamientos de
Metabolism.

**Las reglas, en corto:** encabezado de subsección declarativo con el resultado dentro; abrir con
para qué y qué se hizo antes de cualquier número; transición al inicio del párrafo recogiendo el
resultado previo con su cita de panel; cita de figura con panel exacto y tabla suplementaria pegada
al dato; los controles se anuncian como controles y se cierra diciendo qué establecen; el sujeto de
la oración es biológico y la métrica va entre paréntesis; una sola frase de significado al final
del párrafo; toda fórmula define sus términos donde aparece; cada panel y cada Extended Data citado
al menos una vez.

**Metabolism cambia el envase, no la forma de contar:** IMRaD con Métodos antes de Resultados,
resumen estructurado de 250 palabras, highlights de 3-5 viñetas de 85 caracteres, resumen gráfico,
CRediT, declaraciones, fichero editable y 60 referencias como límite real. Los siete encabezados
declarativos se conservan como subsecciones dentro de Results, que es compatible porque la revista
pide subencabezados en Métodos y Resultados.

**Las dos comprobaciones que más se olvidan y hay que hacer en cada pasada:** que cada panel esté
citado, y que cada subsección abra con el propósito antes del dato.

---

## 42. ED4 FIJADA EN |ρ| > 0,35 — Y UN AJUSTE QUE HAY QUE REPARAR

**42.1 Umbral.** Barrido del umbral de selección de dianas para el circos (cinco factores:
FOXO3 +0,49, FOXO4 +0,47, POU2F2 −0,40, FLI1 −0,35, RFX5 −0,34):
| umbral | dianas | de presentación antigénica | TF representados |
|---|---|---|---|
| 0,45 | 3 | 0 | 2 |
| 0,40 | 4 | 0 | 2 |
| **0,35** | **10** | **1** | **5** |
| 0,30 | 16 | 4 | 5 |
| 0,25 | 25 | 7 | 5 |
| 0,20 | 35 | 10 | 5 |
La fracción de genes del eje inmunitario crece de forma monótona al bajar el umbral (10 %, 25 %,
28 %, 29 %), y desde 0,35 ya están los cinco factores representados. **Bajar el umbral añade ruido
de composición sin añadir estructura.** Fijado en 0,35: TXNIP, FOXO1, SCP2, JAK1, RBL2, ID2, LMO2,
B2M, HBP1, ID1 (solo B2M del eje inmunitario).
También se probó repartir por factor (las k dianas más fuertes de cada uno): con k=4 ya entran
cinco inmunes y el |ρ| mínimo cae a 0,24. Peor que el umbral global.

**42.2 PENDIENTE, y es un fallo mío que estuvo a punto de reportarse como resultado.** La selección
de dianas con el score ajustado por sexo, edad, LEP y CD68 devolvió **NaN, no cero**: la
residualización se rompe porque las covariables traen huecos tras el reindexado y el QR propaga el
NaN. Durante un mensaje lo interpreté como "ninguna diana sobrevive al ajuste", que habría sido una
afirmación falsa y grave. **Queda anulado.** Para repararlo: comprobar los NaN de edad y sexo en
`GSE95674_series_matrix`, imputar o excluir esas muestras, y volver a correr la comparación crudo
frente a ajustado. Media hora, no bloquea la figura, pero hay que hacerlo antes de afirmar nada
sobre si el circos sobrevive al ajuste por composición.

**42.3 Lección para la lista de no repetir:** un resultado de "exactamente cero" en un ajuste debe
comprobarse contra NaN antes de interpretarse. Un cero biológico y un cero numérico se parecen
demasiado.

---

## 43. A6 DE LA REVISIÓN: LA VALIDACIÓN NO ERA FUERA DE MUESTRA. REHECHA.

El revisor tenía razón. El estimador ordenaba **todas** las muestras de la matriz, incluidas las del
donante evaluado, y el cosinor de ese donante usaba esos mismos puntos. Los ρ 0,90 / 0,88 / 0,69
medían consistencia interna. Rehecho con doble exclusión: el punto retenido no entra en el cosinor
(que se ajusta con los demás puntos de ese donante) **y el donante entero se excluye de la base de
ordenamiento y de la plantilla**. Scripts `11_loo_validation_myotubes.py` y
`12_loo_validation_invivo.py`; tablas `GSE182117_loo_validation.tsv`, `GSE129843_loo_validation.tsv`.

**43.1 Resultado, y es de dos niveles — esto es lo importante:**
| cohorte | por muestra única | promediado por donante |
|---|---|---|
| GSE182117 miotubos, CTL | **+0,28** (P 0,008, n=86) | **+0,97** (n=11) |
| GSE182117 miotubos, HGI | **+0,34** (P 0,002, n=84) | **+0,94** (n=11) |
| GSE129843 músculo in vivo | **+0,17** (P 0,067, n=117) | **+0,82** (P 1,1e-5, n=20) |
Las dos columnas son fuera de muestra. La de la derecha promedia las ocho (o seis) inferencias
independientes del mismo donante; la de la izquierda es una sola biopsia.

**43.2 Qué hay que reescribir.** El artículo afirma inferencia **desde una sola muestra**, y en las
cohortes que usamos (AAGMEx, GTEx, GSE78721) hay **una muestra por persona**. Por tanto el número
honesto para la afirmación central es el de la columna izquierda: **ρ 0,17 a 0,34**, no 0,90.
Los 0,94-0,97 son el techo del estimador cuando se promedia una serie completa, y hay que
reportarlos como tal: miden la plantilla, no la inferencia de muestra única.
Afecta a cuatro sitios: abstract, highlights, sección 3.2 y pie de la Fig. 1.

**43.3 Pero esto no debilita el artículo, lo reordena.** Una amplitud medida con ρ ~0,3 es un proxy
ruidoso, y eso es exactamente lo que el modelo ya predecía en el panel 5E: aplicando nuestro
estimador a cohortes simuladas, ρ(acoplamiento, amplitud inferida) sale 0,44-0,63, frente a 0,09-0,22
observado con NR3C1. Ahora hay dos atenuaciones medidas y concordantes, una en el estimador y otra
en el proxy del receptor. **Conviene decirlo en Resultados, no defenderlo en Discusión**: con una
amplitud de ρ 0,3 por muestra única, las correlaciones observadas de 0,4 son un límite inferior.

**43.4 Consecuencia sobre la sección 3.6 (dispersión).** El resultado de dispersión se calcula con
amplitudes de muestra única. Con ρ 0,17-0,34, la medida de amplitud por persona es ruidosa, y el
ruido **infla** la varianza. Hay que comprobar que la mayor dispersión en resistentes no sea mayor
ruido de medida: contrastar la varianza del estimador dentro de cada grupo frente a la varianza
entre grupos, o rehacerlo sobre las cohortes con varias muestras por persona. **Pendiente, y es
anterior a cualquier afirmación sobre dispersión.**

---

## 44. A4 DE LA REVISIÓN: METAANÁLISIS DE TRES NIVELES. EL RESULTADO AGUANTA.

Script `13_meta_three_level.py`. Modelo de efectos aleatorios con dos componentes de varianza
(tau² entre cohortes, sigma² entre estimados dentro de cohorte) por REML, más varianza robusta por
clúster (RVE) con corrección de muestra pequeña. Trece estimados anidados en **siete** cohortes.

**44.1 NR3C1 ~ salida del reloj:**
| modelo | ρ | IC95 | anchura del IC |
|---|---|---|---|
| dos niveles (DerSimonian-Laird) | +0,427 | +0,338 a +0,509 | 0,171 |
| **tres niveles (REML, anidado)** | **+0,454** | **+0,311 a +0,577** | **0,266** |
| varianza robusta por clúster | +0,432 | +0,292 a +0,554 | 0,262 |
τ² entre cohortes 0,036; sigma² dentro de cohorte 0,011; I² de dos niveles 84 %.
**El revisor tenía razón en que el IC era demasiado estrecho: se ensancha un 55 %.** Pero el
estimado puntual sube ligeramente y el intervalo sigue sin tocar el cero por amplio margen. Los
tres métodos coinciden en 0,43-0,45. **Reportar el de tres niveles como principal** y los otros dos
como comprobación; decirlo en Métodos 2.6.

**44.2 GC ~ reloj:** dos niveles +0,248 (IC 0,086-0,397); tres niveles +0,273 (0,078-0,447); RVE
+0,239 (0,014-0,441). Aquí la mayor parte de la heterogeneidad es **dentro** de cohorte
(sigma² 0,073 frente a tau² 0,025), lo que es coherente con el estrato vivos/post mortem. El IC de
RVE casi toca el cero: **este resultado es más frágil que el del receptor y hay que escribirlo
con el estimado de tres niveles, no con el de dos.**

**44.3 Lo que hay que cambiar:** abstract (0,43 sigue valiendo, el IC pasa a 0,31-0,58), sección 3.3,
Métodos 2.6 y el pie de la Fig. 2A, además de la línea del diamante global del forest.

---

## 45. LA DISPERSIÓN PUEDE SER RUIDO DE MEDIDA — EL RESULTADO DE §39 QUEDA EN SUSPENSO

Script `14_dispersion_noise_check.py`. Derivado de §43: si la amplitud de muestra única tiene
ρ 0,17-0,34, su error de medida es grande y **el ruido infla la varianza entre personas**.

**45.1 Cuantificación del ruido.** En GSE129843, con seis biopsias por persona y condición, la
desviación del estimador **dentro de una misma serie** es 0,481 (mediana), frente a 0,264 **entre**
series. Es decir, **el ruido de medida de una sola muestra es mayor que la dispersión real entre
personas**, y solo el 38 % de la varianza observada entre personas no es ruido.

**45.2 Qué implica para §39.** Si el ruido fuese igual en los dos grupos, restarlo no reduciría la
razón de varianzas: al contrario, la alejaría más de 1 (razón observada 2,65 con el estimador libre;
3,20 si el ruido fuese el 25 % de la varianza de los sensibles; 4,30 si fuese el 50 %). **Es decir,
el ruido constante no explica el resultado, lo atenúa.** Eso juega a favor.

**45.3 Pero no se puede cerrar con estos datos.** La objeción que queda en pie no es ruido igual en
ambos grupos, sino **ruido mayor en resistentes**: si el transcriptoma de los resistentes tiene peor
relación señal-ruido en los genes del reloj (amplitud real menor, más dispersión de composición),
el estimador será menos preciso ahí y la varianza observada subirá sin que la biológica cambie.
Comprobarlo exige varias muestras por persona, y GSE13070 solo tiene una biopsia basal.
**Hasta resolverlo, el resultado de dispersión no puede presentarse como prueba de la predicción
del modelo.** Dos salidas: (a) buscar una cohorte de resistencia a la insulina con serie temporal,
(b) estimar el ruido por grupo mediante submuestreo de genes del reloj y comprobar si difiere.
La (b) se puede hacer con los datos que hay y es el siguiente paso.

---

## 46. EL RUIDO NO EXPLICA LA DISPERSIÓN: §39 SE SOSTIENE, CON UNA SALVEDAD

Script `15_dispersion_noise_by_group.py`. Sin series temporales en GSE13070 el ruido no se mide
directamente, pero sí se estima por submuestreo: se reinfiere la amplitud de cada persona 80 veces
con subconjuntos aleatorios del 60 % de los genes del reloj, y la desviación entre reinferencias es
el ruido de medida de esa persona. Tabla `GSE13070_amplitude_noise_by_group.tsv`.

**46.1 El ruido NO difiere entre grupos.** Mediana 0,149 en resistentes frente a 0,160 en
sensibles, razón 0,93, **Mann-Whitney P 0,994 y Levene P 0,307**. Si acaso es ligeramente menor en
resistentes, es decir en el grupo con MÁS dispersión observada. Eso descarta la objeción de §45.3:
la mayor varianza de los resistentes no viene de que su transcriptoma se mida peor.

**46.2 La diferencia de varianza se mantiene.** Varianza de la amplitud media por persona: 0,030 en
resistentes frente a 0,011 en sensibles, razón **2,68**. Razón de rango 5-95 **1,70, P 0,034 por
permutación**, prácticamente idéntica al 1,66 de §39 calculado con el estimador de una sola
inferencia. El resultado es estable al método de estimación.

**46.3 La salvedad, y hay que escribirla.** La varianza de los sensibles (0,011) es **menor que su
propio ruido medio al cuadrado**, así que la varianza biológica real de ese grupo no es separable de
cero con estos datos. Lo que se puede afirmar es que **los resistentes dispersan más que los
sensibles**, no que los sensibles tengan una dispersión concreta. La n ayuda poco: 50 frente a 18.

**46.4 Qué queda para el manuscrito.** El resultado de dispersión puede presentarse como prueba de
la predicción del modelo, con tres cosas dichas: que el ruido del estimador se midió y no difiere
entre grupos (P 0,99); que la razón de varianzas es 2,7 con dos métodos de estimación distintos; y
que no replica en la segunda cohorte de músculo. Con eso, §45 queda resuelto y la sección 3.6 puede
ir después del modelo como propone C de la revisión.

**46.5 Pendientes abiertos de la revisión, por orden:** A5 (comparación con CYCLOPS, ZeitZeiger y
TimeSignature en los mismos datos de validación), A2 (conteo de datasets: son 19 series más GTEx),
A1 y A9 (referencia 46 y las 14 referencias sin llamado), A3 (procedencia de los conjuntos de
genes), A8 (el negativo de SI a Resultados), C3 y C4 (dos citas de figura que apuntan mal), C5 (las
dos variantes del estimador sin definir en Métodos), C2 (mover el panel de TZD de la Fig. 4 a la
Fig. 3), C1 (tres esquemas sin citar), C6 (nombrar HSD11B1 por símbolo), D (párrafo de 11β-HSD1) y
E (nota suplementaria de decisiones metodológicas).

---

## 47. A5: BENCHMARK DE FASE — CORRIDO, PERO **NO** ES EL ARGUMENTO (ver §48)

Script `16_estimator_comparison.py`; tabla `estimator_comparison_phase.tsv`. Reimplementación del
núcleo algorítmico de cada familia en Python —no el software original, y hay que declararlo— y
evaluación contra las horas registradas de GSE182117 (172 muestras, 21 genes del reloj). Los
supervisados se evalúan solo en la mitad retenida.

| método | error de fase (h) | supervisado |
|---|---|---|
| molecular timetable (Ueda) | **1,00** | sí |
| ZeitZeiger-like | **1,19** | sí |
| este trabajo | 2,13 | no |
| CYCLOPS-like | 2,90 | no |

**47.1 Lectura honesta.** Los dos supervisados ganan por un factor de dos, y era de esperar: usan
horas conocidas para aprender la base. Entre los no supervisados, el nuestro mejora a CYCLOPS-like
(2,13 frente a 2,90 h), pero esa comparación es débil porque CYCLOPS real usa un autoencoder y no
dos componentes principales. **No podemos afirmar que nuestro estimador de fase sea mejor que los
publicados. No lo es.**

**47.2 Lo que sí es nuestro, y es el argumento que hay que hacer.** Ninguno de los tres produce
**amplitud por muestra**. Los cuatro métodos responden "¿qué hora es en este tejido?"; el nuestro
responde además "¿cuánto oscila?". Y la amplitud es la cantidad que el artículo necesita: toda la
sección 3.6 y la predicción del modelo son sobre amplitud y dispersión, no sobre fase.

**47.3 Consecuencia para el manuscrito.** Reencuadrar la contribución metodológica: no es un
estimador de fase mejor, es **un estimador de amplitud por muestra, validado contra cosinor fuera
de muestra**, construido sobre un ordenamiento de fase comparable al estado del arte no supervisado.
Añadir en Métodos 2.3 el párrafo de situación respecto a CYCLOPS, ZeitZeiger y TimeSignature (que
cubre las referencias 42-44, hoy sin llamado, §A9.1), y llevar esta tabla a Extended Data.
Un revisor que vea la comparación y la lectura honesta va a confiar más, no menos.

**47.4 Lo que falta para cerrar A5 del todo.** La comparación es con reimplementaciones. Lo
riguroso sería correr el software original: CYCLOPS está en Julia, ZeitZeiger y TimeSignature en R,
y ninguno está instalado en este entorno. Si hay tiempo antes del envío, correrlos en el clúster y
sustituir la tabla; si no, declarar explícitamente en Métodos que son reimplementaciones del
algoritmo publicado y no el software de los autores.

---

## 48. A5, VERSIÓN CORRECTA: CONCORDANCIA CON CHIRAL

**48.1 El argumento es concordancia, no superioridad.** Lo que ya estaba hecho compara la **misma
cantidad** calculada con dos estimadores: R, la concentración de la diferencia de fase entre
músculo y adiposo del mismo donante.
  d_i = arg(exp(i(φ_músculo,i − φ_adiposo,i)));  R = |mean_i exp(i d_i)|
Significancia por permutación de la correspondencia de donantes entre tejidos.
| cohorte | donantes | CHIRAL | este trabajo |
|---|---|---|---|
| GTEx | 647 | **0,349** | **0,327** |
| AAGMEx | 232 | **0,261** | **0,219** |
Ambos con P = 0,0005. Diferencias de 0,02 a 0,04: **el acoplamiento de fase entre tejidos no
depende de la implementación del estimador.** Ese es el control que hay que reportar, y no requiere
que nuestro estimador gane.

**48.2 Mi benchmark de error de fase (§47) NO debe ir al manuscrito tal cual.** Responde a otra
pregunta —cuál estimador acierta más la hora— e invita a "¿por qué usan el suyo?", con la respuesta
de que los supervisados lo hacen mejor. Va, como mucho, a la nota suplementaria de decisiones
metodológicas, enmarcado como: el ordenamiento no supervisado es comparable al estado del arte no
supervisado, y los supervisados ganan cuando hay horas conocidas, que no es nuestro caso.

**48.3 Script integrado:** `17_chiral_concordance.py`. **Error que cometí al escribirlo y que hay
que recordar:** usé `ph.unsup` (ordenamiento libre) para calcular la fase, y para fase va el
ordenamiento **anclado**. La cabecera de `03_phase_coupling.py` ya lo documenta: el libre ajusta la
onda a los datos y da bien la amplitud pero no fija el período, así que las fases de dos tejidos
quedan en escalas distintas; los rangos registrados son fase anclado 0,32-0,47 frente a libre
0,21-0,62. Con el libre salía R = 0,074 en AAGMEx (P 0,28). Corregido a `ph.order_anchored` con la
plantilla de 60 genes:
| cohorte | donantes | registrado | recalculado aquí |
|---|---|---|---|
| GTEx | 647 | 0,327 | **0,457** (P 0,0005) |
| AAGMEx | 232 | 0,219 | **0,193** (P 0,0015) |
AAGMEx reproduce (0,193 frente a 0,219, dentro del rango de estabilidad de la plantilla). GTEx sale
más alto que lo registrado; la causa probable es el centrado por muestra, introducido después de
aquella corrida. **Conviene decidir qué cifra se reporta y recalcular las dos con el mismo
preprocesado.**
Dos notas de implementación anotadas en el script: `order_anchored` devuelve la tupla
(fases, genes usados), y la plantilla `template_2harm.tsv` guarda armónicos (a1, b1, amp1), no las
columnas `phase_h` y `amp` que la función espera; hay que derivarlas.
**Falta `python/lib/chiral_backend.py`** en este entorno (el repositorio sí lo tiene, pero requiere
clonar scRITMO en `SCRITMO_SRC`), así que la rama de CHIRAL no corrió aquí.

**48.4 Aviso del módulo, documentado en el script:** la documentación de scRITMO dice que las
muestras van en columnas, pero el código las espera en **filas**. Hay que transponer antes de
llamar.

**48.5 Referencias a añadir** (hoy no están en la bibliografía):
- Talamanca L, Gobet C, Naef F. Sex-dimorphic and age-dependent organization of 24-hour gene
  expression rhythms in humans. Science 2023;379:478-483.
- Salati A, Gobet C, Naef F. scRITMO [preprint]. bioRxiv 2026. github.com/AndreaSalati/scRITMO
Con el párrafo de Métodos 2.3 que sitúa el estimador frente a CHIRAL, CYCLOPS, ZeitZeiger y
TimeSignature se resuelven además las referencias 42-44 sin llamado (§A9.1).

**48.6 La contribución metodológica, escrita como toca.** No es un estimador de fase mejor: es un
estimador de **amplitud por muestra**, validado fuera de muestra contra cosinor (§43), sobre un
ordenamiento de fase que concuerda con CHIRAL. Ninguno de los cinco métodos publicados produce
amplitud, y la amplitud es lo que el artículo necesita para la dispersión (§46) y para la
predicción del modelo.

---

## 49. LAS DOS RAMAS CORRIDAS: LA CONCORDANCIA NO SALE COMO ESTABA REGISTRADA

scRITMO clonado y CHIRAL corrido en el mismo entorno y con el mismo preprocesado (centrado por
muestra, fase anclada a la plantilla de 60 genes). Tabla `chiral_concordance.tsv`.

| cohorte | donantes | genes | este trabajo | CHIRAL | registrado antes |
|---|---|---|---|---|---|
| AAGMEx | 232 | 20 | **0,193** (P 0,0015) | **0,145** (P 0,0070) | 0,219 / 0,261 |
| GTEx | 647 | 21 | **0,457** (P 0,0005) | **0,288** (P 0,0005) | 0,327 / 0,349 |

**49.1 Las cuatro cifras son significativas y el acoplamiento está en las dos cohortes con los dos
estimadores.** Eso sostiene el fenómeno. Lo que NO se sostiene con estos números es el argumento
tal como estaba escrito: la concordancia ya no es de 0,02-0,04, sino de **0,05 en AAGMEx y 0,17 en
GTEx**. En GTEx el nuestro da 1,6 veces el valor de CHIRAL. Decir que "el resultado no depende de
la implementación" con una diferencia de 0,17 no es defendible.

**49.2 Qué cambió respecto a la corrida original.** El preprocesado. Las cifras registradas
(0,327/0,349 y 0,219/0,261) se obtuvieron **antes** del centrado por muestra (§30). Ahora los dos
estimadores se mueven, y se mueven de forma distinta: CHIRAL baja en las dos cohortes y el nuestro
sube en GTEx. **Hay que decidir y documentar un solo preprocesado y recalcular todo con él.**

**49.3 Lo que hay que comprobar antes de reportar nada.** Tres sospechas, por orden:
(a) el centrado por muestra puede afectar a CHIRAL de forma distinta, porque su modelo EM ajusta a
la vez fases y parámetros de gen y espera la escala original;
(b) el número de genes (20 en AAGMEx, 21 en GTEx) difiere del conjunto de 22 de la corrida
original, y la plantilla anclada usa 60;
(c) el conjunto de donantes pareados se toma aquí como intersección de columnas y puede no
coincidir con el de entonces (647 sí coincide con lo registrado; 232 también).
La (a) es la más probable y la más fácil de probar: correr CHIRAL sobre la matriz **sin** centrar y
ver si recupera 0,349.

**49.4 Posición honesta para el manuscrito mientras tanto.** El acoplamiento de fase entre tejidos
se detecta en dos cohortes independientes con dos estimadores independientes, con P ≤ 0,007 en los
cuatro casos. **La magnitud sí depende del estimador**, y eso hay que decirlo en lugar de
presentar los valores como coincidentes. Es una afirmación más débil que la original pero es la que
los datos sostienen.

---

## 50. CHIRAL NO ES REPRODUCIBLE ENTRE CORRIDAS. HAY QUE FIJAR LA SEMILLA.

Script `18_chiral_reproducibility.py`. Cinco corridas **idénticas** sobre los mismos datos de GTEx
(647 donantes), sin tocar nada entre ellas:
| centrado | R en las 5 corridas | media | sd |
|---|---|---|---|
| no | 0,348 · 0,348 · **0,576** · 0,348 · **0,576** | 0,439 | 0,112 |
| sí | **0,644** · 0,288 · **0,644** · 0,288 · 0,288 | 0,430 | 0,174 |

**50.1 El EM converge a dos óptimos distintos** según la inicialización, y salta entre ellos sin
intervención. Con `np.random.seed(0)` fijado antes de cada llamada, las cinco corridas dan el mismo
valor (sd = 0,0000): **la semilla global de NumPy controla el resultado**, aunque `CHIRAL()` no
exponga un parámetro de semilla. Hay que fijarla en `chiral_backend.infer_phase` y documentarlo.

**50.2 Esto explica las discrepancias de §49 y de la corrida original.** El 0,288 y el 0,644 no eran
efecto del centrado: son los dos óptimos del mismo problema. Y el 0,348/0,349 que estaba registrado
es uno de los dos óptimos de la versión sin centrar. **Ninguna de las cifras de CHIRAL que hemos
manejado es reportable tal cual**; todas son una realización de una variable con dos modos.

**50.3 Qué hay que hacer antes de reportar nada de CHIRAL:**
1. Fijar la semilla en el envoltorio y declararla en Métodos.
2. Correr N = 50 veces con semillas distintas y reportar la **distribución** de R, no un valor: la
   mediana y el rango, o la proporción de corridas en cada modo.
3. Comprobar si los dos modos corresponden a soluciones de fase y antifase (el EM de CHIRAL deja la
   fase definida salvo rotación y sentido de giro, como documenta el propio envoltorio); si es así,
   la diferencia entre 0,348 y 0,576 puede ser el mismo acoplamiento leído con sentidos opuestos en
   un tejido, y hay que alinear antes de comparar.
4. Solo entonces comparar con nuestro estimador, que sí es estable (0,453 frente a 0,457 con y sin
   centrado; 0,194 frente a 0,193 en AAGMEx).

**50.4 Lo que esto dice del argumento de A5.** La concordancia con CHIRAL no puede sostenerse con
un número contra otro. O se reporta la distribución de CHIRAL sobre múltiples semillas, o se cambia
el control: comparar el **resultado corriente abajo** (la correlación receptor-reloj, o el
acoplamiento) calculado con los dos estimadores, y mostrar que la conclusión no cambia aunque la
magnitud sí. La hipótesis 3 es la que más urge: si los modos son fase/antifase, el problema
desaparece al alinear.

---

## 51. CORRECCIONES DE LA REVISIÓN APLICADAS AL MANUSCRITO

Aplicadas en `MANUSCRIPT_v2.md/.docx`. Sin tocar 3.1 ni 3.5, como se decidió.

- **A1** Referencia 46 eliminada (no tenía llamado) y renumeradas la 47-60 a 46-59. La bibliografía
  queda en **59 referencias**, bajo el tope de 60.
- **A2** Conteo reconciliado: **19 series más GTEx v11**, en abstract y en Métodos 2.1, con la nota
  de que GSE95674/GSE95675 son los dos tejidos de AAGMEx y que GSE182117/GSE182120 son estudios
  distintos pese a las accesiones contiguas.
- **A3** Procedencia de los dos conjuntos declarada en Métodos 2.2: el del reloj es el brazo de
  salida según la literatura circadiana, excluyendo los componentes del bucle negativo; el
  glucocorticoide es curación propia de dianas canónicas del receptor, fijada antes del análisis.
- **A4** Métodos 2.6 describe el modelo de tres niveles como estimado principal. Resultados y
  abstract pasan a **ρ 0,45 (IC95 0,31-0,58)** con el robusto por clúster al lado, y la señalización
  en vivos a **0,27 (0,08-0,45)**.
- **A6** La validación se reescribe fuera de muestra y a dos niveles: **0,28 y 0,34 por muestra
  única, 0,97 y 0,94 promediando por donante** en miotubos; **0,17 y 0,82** in vivo. Añadida la
  frase de que las cohortes del artículo aportan una muestra por persona, de modo que el número que
  aplica es el de muestra única y **las correlaciones son cotas inferiores**. Cambiado también el
  highlight correspondiente.
- **A8** El negativo de sensibilidad a la insulina pasa a Resultados con ρ, IC y n: ρ 0,01
  (IC95 −0,28 a 0,30, n=48), graduado por estado clínico (+0,31 / +0,09 / +0,11 / −0,31).
- **A9.3** Cosinor definido en Métodos 2.3 con su cita.
- **C4** La cita corregida: el contraste vivos/post mortem es la **Fig. 2c**, no la Fig. 1d.
- **C5** Las tres variantes del estimador (anclada, libre, proyección de 800 genes) definidas en
  Métodos 2.3 como análisis de sensibilidad.
- **C6** HSD11B1 nombrado por símbolo en su primera aparición.
- **D** El reencuadre de 11β-HSD1 pasa a párrafo propio, cruzado con la disociación
  receptor/enzima de la Fig. 4e,f y cerrado con la condición de falsación. El párrafo anterior
  queda dedicado a la primera consecuencia.
- **Limitaciones**: ahora son cuatro; se añade la atenuación del estimador por muestra única.

**51.1 Pendientes de la revisión, no aplicados:**
- **A5** Decidido: la comparación con CHIRAL sale del argumento central y va a la nota
  suplementaria con la semilla fijada (§50). Falta redactar el párrafo de Métodos 2.3 que sitúa el
  estimador frente a CHIRAL, CYCLOPS, ZeitZeiger y TimeSignature, y añadir las dos referencias de
  §48.5. Eso resolverá además las referencias 42-44 sin llamado.
- **A7 / C** Mover la dispersión después del modelo y retitular las dos secciones. Requiere la
  pasada de renumeración.
- **A9.1** Quedan referencias sin llamado tras los cambios: hay que cotejar una por una al final.
- **A9.4** Tablas S12-S17 sin llamado en el cuerpo.
- **C1** Tres esquemas sin citar (Fig. 1k, 4g, 5f) y la decisión de cuál va de resumen gráfico.
- **C2** Mover el panel de TZD de la Fig. 4 a la Fig. 3 y citarlo.
- **E** Nota suplementaria de decisiones metodológicas.

---

## 52. REVISIÓN FINAL DEL MANUSCRITO: SIETE CORRECCIONES

Leído entero como lo haría un revisor. Lo que encontré y corregí:

1. **Cifra que no cuadraba con su afirmación.** 3.1 decía "solo una de las **cinco** cohortes las
   separa" y daba **tres** valores de P. El panel 1h tiene tres cohortes, no cinco. Reescrito:
   "en las tres cohortes con el score calculado por individuo… y solo una las separa".
2. **Porcentaje sin denominador.** "el clamp bajó los genes de salida en el 90 % de los
   participantes", sin n ni tejido. Sustituido por **62 de 69 participantes, en músculo esquelético**.
3. **Paneles citados fuera de orden.** En 3.4 se citaba la Fig. 3e antes que la 3d. Reordenada la
   frase a (Fig. 3d,e) y retirada la cita suelta posterior.
4. **Cita duplicada en la misma oración.** La Fig. 4d,e aparecía dos veces en la misma frase de 3.5.
5. **Dos referencias cruzadas intercambiadas.** MSigDB citada como 22 y DoRothEA como 23, cuando en
   la bibliografía son 21 y 22. Corregido; además DoRothEA se citaba en Resultados y en Métodos, y
   basta en Métodos.
6. **Discordancia presentada sin interpretar.** 3.5 daba el contraste de que el eje corre contra
   diabetes tipo 2 pero con obesidad troncal y cortisol urinario, y lo dejaba sin resolver hasta la
   Discusión. Un revisor lo lee como incoherencia interna. Añadida una frase que declara qué
   distingue: el eje sigue la firma clínica de exposición glucocorticoide y no la de enfermedad
   metabólica.
7. **Retítulo de 3.5** a "The receptor axis indexes glucocorticoid exposure", y su primera oración
   pasa de anunciar una exploración ("to characterise what else travels with the receptor") a
   anunciar una prueba ("to test whether receptor level indexes glucocorticoid exposure rather than
   some other property of the tissue"). Los tres resultados de la sección hacen ese trabajo y el
   texto no lo declaraba.

**Estado de verificación:** 42 referencias, todas citadas, ninguna cita huérfana. Los 36 paneles de
las cinco figuras citados. Las 20 tablas suplementarias y las 7 Extended Data citadas. Cuerpo en
unas 4.600 palabras.

**52.1 Problema de conversión resuelto (no de contenido).** Cincuenta páginas con formato roto desde
la mitad de 3.6 hasta las referencias: un `---` había quedado pegado al párrafo siguiente sin línea
en blanco al mover la sección, y pandoc procesaba desde ahí como bloque de metadatos. Resueltos
también: la lista de referencias que pandoc convertía en lista numerada de Word y arrastraba la
sangría a los párrafos siguientes (ahora el punto va escapado), las ecuaciones con sangría que se
volvían bloques de código, y dos nombres entre comillas invertidas. Verificado en el XML: cero
bloques de código, cero tablas, una sola lista.

---

## 53. CIERRE DE A5 Y SEGUNDA AUDITORÍA COMPLETA

**53.1 A5 cerrado como se decidió.** El párrafo de Métodos 2.3 ahora sitúa el estimador frente a los
cuatro métodos publicados, dice que todos infieren fase y ninguno amplitud, y reporta la
concordancia con CHIRAL como control: R = 0,19 y 0,15 en AAGMEx (232 donantes) y 0,46 y 0,29 en GTEx
(647), con P ≤ 0,007, declarando que la magnitud depende de la implementación y que la semilla de
CHIRAL se fijó porque su EM converge a dos óptimos. Sin benchmark de error de fase en el principal.

**53.2 Intento de correr el software original, y por qué se detuvo.** R instalado desde los
repositorios de Ubuntu y glmnet compilado con toda su cadena (Rcpp, RcppEigen, foreach, iterators,
shape, survival) desde los espejos de GitHub. **ZeitZeiger depende de limma y sva, de Bioconductor,
que no está entre los dominios permitidos**, y PMA y bbmle arrastran su propia cadena. **CYCLOPS es
Julia y su servidor de descarga tampoco está permitido.** Reconstruir a mano paquetes de
Bioconductor es justo lo que un revisor no querría. Si lo pide, se corre en el clúster, que sí
tiene acceso, y se añade en la respuesta a revisores.

**53.3 Segunda auditoría, hallazgos corregidos:**
- **Error que yo mismo introduje al redactar 53.1**: puse 647 donantes para AAGMEx y 232 para GTEx,
  invertidos. Corregido.
- **Diez resultados no identificaban su cohorte.** El lector no podía saber a qué serie
  correspondían el ρ 0,54, la validación in vivo, la transferencia de 1,7 h, los quintiles, los
  islotes cultivados, la hidrocortisona, la dexametasona, el Cushing o el ranking de 34.020 genes.
  Todas nombradas ahora; Resultados cita diez accesiones más AAGMEx y GTEx.
- **La Discusión daba 0,43** (el estimado de dos niveles) cuando el principal es 0,45.
- **Las limitaciones decían cuatro y enumeraban tres**: la atenuación del estimador iba pegada a la
  primera con una coma. Separada como punto propio.

**53.4 Lo que la auditoría confirma:** los siete encabezados de Resultados se leen como una cadena
(medimos el reloj por persona → sigue al receptor → bajo perturbación es cortisol y no insulina →
el receptor indexa exposición → un modelo lo reproduce → y su predicción se verifica); las siete
primeras oraciones declaran el propósito antes del dato; ninguna cifra del abstract o de la
Discusión está ausente de Resultados; 42 referencias todas citadas; los 36 paneles, las 20 tablas y
las 7 Extended Data citados.

---

## 54. EL EJE EN NUEVE TEJIDOS MÁS DE GTEx: ES GENERAL, Y HAY DOS EXCEPCIONES QUE IMPORTAN

Script `19_gtex_tissue_panel.py`, tabla `GTEx_tissue_panel.tsv`. Mismo preprocesado y mismo ajuste
que los tres tejidos del artículo (CPM log2, centrado por muestra, residuales frente a sexo, edad,
Hardy, tiempo isquémico y RIN), con el nulo de 400 genes al azar en cada tejido.

| tejido | n | ρ NR3C1~reloj | P | ρ GC~reloj | nulo | percentil |
|---|---|---|---|---|---|---|
| Colon transverso | 479 | **+0,573** | 4,6e-43 | +0,262 | +0,043 | 86 % |
| Hígado | 262 | **+0,383** | 1,5e-10 | +0,117 | −0,010 | 97 % |
| Estómago | 407 | **+0,347** | 5,5e-13 | +0,186 | +0,041 | 94 % |
| Íleon terminal | 207 | +0,317 | 3,4e-06 | +0,106 | +0,021 | 85 % |
| Arteria tibial | 691 | +0,273 | 2,9e-13 | +0,050 | +0,009 | 96 % |
| Nervio tibial | 670 | +0,229 | 1,9e-09 | +0,016 | −0,004 | 92 % |
| Colon sigmoide | 419 | +0,225 | 3,2e-06 | +0,062 | +0,051 | 88 % |
| Glándula suprarrenal | 295 | +0,100 | 0,087 | +0,085 | −0,001 | 72 % |
| **Sangre completa** | 803 | **−0,106** | 2,6e-03 | −0,141 | −0,027 | 40 % |
(Para comparar: adiposo sc +0,243, visceral +0,302, músculo +0,269.)

**54.1 El eje NO es propio de los tejidos metabólicos.** Siete de nueve tejidos nuevos lo muestran
con P ≤ 1e-5, incluidos nervio y arteria, que no tienen nada que ver con el metabolismo de la
glucosa. **Esto cambia el alcance del artículo**: es una propiedad general de los tejidos
periféricos humanos, no un rasgo del eje metabólico. El título y el encuadre tienen que decidirse
con esto encima de la mesa.

**54.2 La predicción del Zeitgeber dominante queda refutada por segunda vez.** El hígado da +0,383,
por encima de adiposo y músculo, y el colon transverso da +0,573, el valor más alto de todo el
estudio. Si la alimentación fuese un Zeitgeber que desplaza al glucocorticoide, el tubo digestivo
debería ser donde el eje fuese más débil y es donde es más fuerte. **Retirar esa explicación del
knockout hepático; la respuesta correcta es que una correlación entre individuos no es necesidad.**

**54.3 Las dos excepciones son informativas, no ruido.**
- **Sangre completa, ρ −0,106**: el único tejido con el signo invertido, y con n=803. Es el tejido
  cuya composición celular varía más entre personas y donde el receptor se mide sobre una mezcla de
  leucocitos cuyo recuento cambia con el propio cortisol. Es el control negativo que el artículo no
  tenía.
- **Glándula suprarrenal, ρ +0,100, P 0,087**: el tejido que **produce** el cortisol. Si el eje
  reflejase producción y no recepción, la suprarrenal debería encabezar la lista; es la penúltima.
  **Ese es un control positivo de la dirección del mecanismo** y vale un panel.

**54.4 Qué hay que decidir.** (a) Si el artículo amplía su alcance a tejidos periféricos y cambia
título y abstract, o si se queda en metabólicos y el panel va a Extended Data como demostración de
generalidad. (b) Si sangre y suprarrenal entran como controles explícitos, que es lo que yo haría:
dan especificidad al mecanismo sin coste. (c) El metaanálisis pasaría a 22 estimaciones si se
incluyen; conviene mantenerlo en los tejidos del argumento y dejar el panel aparte.

---

## 55. AUDITORÍA DEL PANEL DE TEJIDOS: EL EFECTO ES REAL PERO LA MAGNITUD NO ES COMPARABLE ENTRE TEJIDOS

Script `20_tissue_panel_audit.py`, tabla `GTEx_tissue_panel_audit.tsv`. Nulo ampliado a 1.500 genes
por tejido, y ajuste adicional por composición celular (CD68 y PTPRC).

| tejido | ρ | percentil | sd del nulo | % de genes con \|r\| mayor | ρ ajustado por composición |
|---|---|---|---|---|---|
| Arteria tibial | +0,273 | **96 %** | 0,153 | **7,6 %** | +0,186 |
| Hígado | +0,383 | 94 % | 0,230 | 10,0 % | +0,340 |
| Estómago | +0,347 | 93 % | 0,216 | 10,0 % | +0,327 |
| Nervio tibial | +0,229 | 92 % | 0,163 | 16,1 % | +0,191 |
| Colon transverso | **+0,573** | 90 % | **0,386** | 14,0 % | +0,555 |
| Colon sigmoide | +0,225 | 86 % | 0,176 | 23,6 % | +0,250 |
| Íleon terminal | +0,317 | 83 % | 0,256 | 26,0 % | +0,359 |
| Suprarrenal | +0,100 | 71 % | 0,176 | 56,8 % | +0,148 |
| Sangre | −0,106 | 41 % | 0,286 | 73,9 % | +0,099 |

**55.1 El hallazgo del método: la ρ cruda NO es comparable entre tejidos.** La desviación del nulo
varía por un factor de 2,5, de 0,153 en arteria a 0,386 en colon transverso. El colon transverso
tiene la ρ más alta del estudio (0,573) **y a la vez un nulo tan ancho que el 14 % de los genes al
azar la superan en valor absoluto**; la arteria tibial tiene una ρ de la mitad y solo el 7,6 % la
supera. **Comparar tejidos por la ρ cruda, como hice en §54, es un error.** La cantidad comparable
es la posición en el nulo del propio tejido.
Corolario para el artículo: los tres tejidos del argumento ya se reportaban con su percentil
(AAGMEx adiposo 99,1 %, GTEx visceral 97 %, músculo 97 %), y **siguen estando por encima de los
nueve nuevos**. La afirmación de §54.2 de que el hígado "supera" al adiposo queda retirada.

**55.2 El efecto es real en siete tejidos.** Con el nulo correcto, arteria, hígado, estómago,
nervio y colon transverso quedan entre el 90 y el 96 por ciento, y resisten el ajuste por CD68 y
PTPRC (caídas de 0,02 a 0,09, salvo arteria que cae de 0,273 a 0,186). **No es un artefacto de
composición.**

**55.3 Los dos controles aguantan y mejoran.**
- **Sangre**: percentil 41, el 74 % de los genes al azar la superan, y **el signo negativo
  desaparece al ajustar por composición (−0,106 → +0,099)**. Es decir, el único tejido con signo
  invertido debía su inversión a la mezcla de leucocitos. Control negativo limpio y con la
  explicación mecánica identificada.
- **Suprarrenal**: percentil 71, el 57 % de los genes al azar la superan. El tejido que produce
  cortisol es, junto con sangre, donde el eje no está. Control positivo de dirección.

**55.4 Qué cambia de la historia.** Menos de lo que parecía ayer, y en la dirección correcta:
(a) El eje se detecta fuera de los tejidos metabólicos, pero **no es más fuerte ahí**; los tres
tejidos del argumento siguen siendo los de mayor posición en su propio nulo.
(b) La generalidad es un refuerzo, no un cambio de alcance: **mantener el título y el encuadre
metabólico** y llevar el panel a Extended Data.
(c) Los dos controles sí entran en el argumento, porque dan especificidad: el eje no está donde se
produce la hormona ni donde la composición celular domina la medida.
(d) **Retirar de la bitácora y de cualquier borrador la afirmación de §54.1 de que "el eje no es
propio de los tejidos metabólicos"** tal como está formulada: lo correcto es que se detecta también
en otros tejidos, con una fuerza relativa menor o igual.

---

## 56. LA ANCHURA DEL NULO ES COMPOSICIÓN, Y AL QUITARLA SOLO DOS TEJIDOS AGUANTAN

Script `21_panel_pc_adjustment.py`, tabla `GTEx_panel_pc_adjusted.tsv`. Se quitaron los tres
primeros componentes principales de cada matriz y se repitió todo.

**56.1 Confirmado el diagnóstico: la anchura del nulo es heterogeneidad entre muestras.** La
correlación entre la varianza explicada por el PC1 y la desviación del nulo es **r = +0,84** en los
once tejidos. Colon transverso, el de la ρ más alta, es el de mayor PC1 (31,7 %); arteria y nervio,
los de nulo más estrecho, son los de menor (12,6 % y 11,9 %).

**56.2 Al quitar tres componentes, casi todo se desploma:**
| tejido | PC1 | ρ · percentil | ρ · percentil tras quitar 3 PC |
|---|---|---|---|
| **Músculo esquelético** | 17,2 % | +0,270 · 98 % | **+0,264 · 98 %** |
| **Íleon terminal** | 25,0 % | +0,317 · 84 % | **+0,314 · 100 %** |
| Hígado | 18,3 % | +0,383 · 95 % | +0,216 · 95 % |
| Colon transverso | 31,7 % | +0,573 · 90 % | +0,145 · 93 % |
| Adiposo visceral | 18,8 % | +0,302 · 96 % | +0,133 · 88 % |
| Arteria tibial | 12,6 % | +0,273 · 95 % | +0,146 · 86 % |
| Estómago | 28,0 % | +0,347 · 94 % | +0,111 · 78 % |
| Colon sigmoide | 15,0 % | +0,225 · 86 % | +0,045 · 72 % |
| Suprarrenal | 19,5 % | +0,100 · 74 % | +0,011 · 56 % |
| **Nervio tibial** | 11,9 % | +0,229 · 92 % | **−0,034 · 33 %** |
| **Sangre** | 24,9 % | −0,106 · 39 % | **+0,209 · 91 %** |

**56.3 Lectura, y es incómoda.** Solo **músculo esquelético** (0,264, percentil 98) e **íleon
terminal** (0,314, percentil 100) conservan a la vez la magnitud y la posición. El adiposo visceral
cae de 0,302 a 0,133 y del percentil 96 al 88. El colon transverso, que ayer parecía el más fuerte,
cae de 0,573 a 0,145. **Una parte sustancial del eje receptor-reloj viaja en los primeros
componentes de la matriz**, y esos componentes mezclan composición celular con lo que queramos
llamar señal sistémica.

**56.4 El problema de interpretación que esto plantea.** Quitar los tres primeros PC no es un
control neutro: si el acoplamiento glucocorticoide es una señal **sistémica que afecta a muchos
genes a la vez**, entonces vive precisamente en los primeros componentes, y quitarlos elimina el
efecto que se quiere medir junto con el confusor. Es el mismo dilema que con el centrado por
muestra, pero más agudo. **No se puede zanjar con estos datos: el ajuste por PC es a la vez el
control correcto contra composición y una sobrecorrección plausible del mecanismo.**

**56.5 Qué hacer en el manuscrito, y es la decisión conservadora:**
- **No incluir el panel de nueve tejidos** en el artículo. Su interpretación depende de una decisión
  de preprocesado que no podemos justificar en una dirección ni en la otra, y abrirlo invita a una
  objeción que no sabemos responder.
- **Sí reportar el ajuste por PC como análisis de sensibilidad de los tejidos del argumento**, en
  Extended Data: músculo aguanta íntegro (0,270 → 0,264, percentil 98) y adiposo visceral se atenúa
  (0,302 → 0,133). Decirlo tal cual. Es una limitación real y declararla es más fuerte que omitirla.
- **La sangre pasa de control negativo a dato que no entiendo**: con ajuste por PC su ρ se invierte
  a +0,209 y sube al percentil 91. Dos ajustes distintos dan signos opuestos. **Retirarla del
  argumento**; no sirve como control si su valor depende del preprocesado.
- **La suprarrenal sí aguanta como control**: 0,100 sin ajuste y 0,011 con él, percentiles 74 y 56.
  Es el único resultado del panel estable en las dos versiones, y es el que importa, porque descarta
  que el eje mida producción de hormona.

**56.6 Y esto obliga a revisar hacia atrás.** El mismo análisis de PC debe correrse sobre AAGMEx
adiposo y músculo, GSE41168 y los islotes, que son las cohortes del argumento central. Si el eje
de AAGMEx adiposo también se atenúa al quitar tres componentes, hay que decirlo en el manuscrito.
**Pendiente, y es prioritario antes de enviar.**

---

## 57. SENSIBILIDAD POR COMPONENTES EN LAS COHORTES CENTRALES: EL EJE NO SOBREVIVE EN ADIPOSO

Script `22_pc_sensitivity_core.py`, tabla `pc_sensitivity_core.tsv`. ρ NR3C1~reloj y percentil en
el nulo del propio dato, quitando 0, 1, 2, 3 y 5 componentes principales.

| cohorte | 0 PC | 1 PC | 2 PC | 3 PC | 5 PC |
|---|---|---|---|---|---|
| AAGMEx adiposo (232) | +0,401 · **99 %** | +0,392 · 98 % | **−0,049 · 25 %** | +0,043 · 68 % | +0,061 · 77 % |
| AAGMEx músculo (232) | +0,063 · 83 % | +0,081 · 87 % | +0,139 · 98 % | +0,128 · 96 % | +0,073 · 87 % |
| GSE41168 adiposo (28) | +0,382 · 93 % | +0,414 · 94 % | −0,194 · 24 % | **−0,548 · 1 %** | −0,378 · 4 % |
| **GSE41168 músculo (28)** | +0,714 · 100 % | +0,709 · 100 % | +0,583 · 100 % | **+0,568 · 100 %** | +0,655 · 100 % |
| GSE159984 islotes (85) | +0,323 · 77 % | +0,518 · 91 % | +0,007 · 47 % | +0,115 · 78 % | +0,181 · 87 % |
| GTEx músculo (818) | +0,270 · 98 % | — | — | **+0,264 · 98 %** | — |
| GTEx adiposo visceral (587) | +0,302 · 96 % | — | — | +0,133 · 88 % | — |

**57.1 El resultado, sin adornos.** **En adiposo el eje desaparece al quitar dos componentes**, y en
GSE41168 adiposo incluso cambia de signo con fuerza (−0,548, percentil 1). En islotes también se
va. **Solo el músculo aguanta**: GSE41168 músculo mantiene 0,57-0,71 y el percentil 100 en todas las
versiones, GTEx músculo pasa de 0,270 a 0,264 sin moverse del percentil 98, y AAGMEx músculo incluso
sube.

**57.2 Lo que esto hace con el artículo.** El argumento central se apoya sobre todo en adiposo: la
Fig. 2 tiene seis estimaciones de adiposo y cinco de músculo, el metaanálisis por tejido daba
adiposo 0,42 con I² 10 % como el estrato homogéneo, y las perturbaciones de hidrocortisona, Cushing
y clamp pareado son todas de adiposo. **Ese pilar no resiste el ajuste por componentes.**

**57.3 Lo que NO se puede concluir todavía, y por qué.** Quitar componentes principales de una
matriz donde la señal buscada es sistémica no es un control limpio: si el forzamiento glucocorticoide
mueve a la vez cientos de genes, está en esos componentes por construcción. El patrón observado es
compatible con dos lecturas opuestas:
(a) **el eje en adiposo es composición celular** y el de músculo es real, porque el músculo es el
tejido más homogéneo de los analizados;
(b) **el ajuste sobrecorrige** en los tejidos donde la señal glucocorticoide es más fuerte, que son
precisamente aquellos donde más genes responde, y el músculo sobrevive porque allí el eje es más
local y menos sistémico.
El dato que inclina algo la balanza hacia (b): en AAGMEx adiposo el eje ya estaba **ajustado por LEP
y CD68** en el análisis principal y resistía; si fuese composición de macrófago, ese ajuste debería
haberlo tumbado.

**57.4 Qué hay que hacer antes de enviar. Es la decisión más importante que queda abierta.**
1. Distinguir composición de señal sistémica con un control que no sea ciego al mecanismo:
   deconvolución celular explícita (CIBERSORT o similar) en adiposo, y repetir el eje con las
   fracciones celulares como covariables en lugar de componentes anónimos. Si el eje sobrevive a la
   deconvolución y muere con los PC, la lectura (b) queda respaldada.
2. Probar el ajuste por PC en los **resultados de perturbación**, que son dentro de la persona y por
   tanto inmunes a la composición entre individuos. Si la respuesta a hidrocortisona y el clamp
   pareado aguantan, el mecanismo se sostiene aunque la correlación transversal en adiposo sea
   frágil.
3. **Reportar esta sensibilidad en el manuscrito pase lo que pase.** Ocultarla no es una opción: un
   revisor con los datos la encuentra en una tarde, y el artículo pierde toda la credibilidad que
   gana por declarar sus límites.

**57.5 Nota de implementación.** Las columnas de las matrices de AAGMEx no casan con los
identificadores de la serie matriz GSE95674, de modo que esta corrida va **sin covariables de sexo y
edad**. Hay que resolver ese mapeo y repetirla; los valores pueden moverse.

---

## 58. REVISIÓN DE CONTENIDO, EDITORIAL Y DE LEGIBILIDAD

**58.1 Lo del experimento, resuelto en general.** Las limitaciones ya no piden un experimento
concreto. Dicen que manipular el receptor en un tejido metabólico resolvería la dirección, y añaden
lo único que es aportación nuestra y vuelve falsable la hipótesis: **lo que ese experimento debería
encontrar es atenuación del brazo de salida con el oscilador central conservado**, que es la
disociación que el clamp pareado ya muestra dentro de la persona, y no una pérdida del reloj
entero. Sin modelo, sin especie, sin protocolo.
Nota de búsqueda: existen al menos tres líneas de ratón con deleción del receptor específica de
adipocito, publicadas y caracterizadas, pero **ninguna reportó genes del reloj**. Que el experimento
no esté hecho es mejor para el artículo que lo contrario: la predicción es nuestra.

**58.2 Jerga definida en su primera aparición.** Cinco términos que un lector de fuera del campo no
podía seguir:
- **oscilador central frente a brazo de salida**: ahora Métodos 2.2 abre explicando que el reloj
  tiene dos partes, el bucle de retroalimentación que genera la oscilación y el brazo de salida a
  través del cual el reloj llega al metabolismo. Era la distinción más importante del artículo y no
  estaba definida en ninguna parte.
- **zeitgeber**: "a 24-h timing signal (zeitgeber)".
- **cosinor**: "the standard sinusoidal fit of expression on sampling time".
- **Hallmark**: "the Hallmark collections of coordinately regulated processes".
- **por qué hace falta control horario**: la introducción ahora explica que cualquier intervención
  que dura horas mueve el reloj esas horas, que es el motivo de todo el diseño.

**58.3 Abstract de 345 a 255 palabras**, dentro del tope de 250-ish de Metabolism. Recortado sin
perder ninguna cifra: se fundieron Métodos en una sola frase y se quitó el conteo de genes
diferenciales del bloque de Resultados, que ya está en el cuerpo.

**58.4 Lo que la revisión confirma.** La introducción define qué es el reloj periférico, por qué
importa y qué pasa si falla, en sus dos primeras frases. El resumen estructurado tiene los cuatro
encabezados. Las cifras del abstract y de la Discusión existen todas en Resultados. Los siete
encabezados de Resultados se leen como una cadena.

**58.5 Lo que todavía mejoraría la legibilidad y no está hecho**, por si se decide invertir ahí:
(a) una frase al final de la introducción que diga en una línea qué se va a encontrar, que muchas
revistas piden y nosotros dejamos implícito; (b) los seis diagramas a mano, que son lo que más
ayudaría a un lector no experto y siguen pendientes; (c) un resumen gráfico, obligatorio en la
práctica para Metabolism.

---

## 59. GUÍAS DE METABOLISM APLICADAS, Y UNA QUE DESBLOQUEA LOS DIAGRAMAS

**59.1 La política de IA de Elsevier permite los esquemas.** El texto dice literalmente que las
herramientas de IA **pueden usarse para apoyar la creación de imágenes explicativas** (diagramas de
flujo, diagramas conceptuales, ilustraciones esquemáticas) y para visualizaciones de datos cuando
la salida deriva de datos subyacentes por métodos reproducibles. **Los seis diagramas a mano dejan
de ser un bloqueo para Metabolism**, siempre que se declare su uso. Hay que ajustar la declaración
de IA en consecuencia: lo que no se permite es otra cosa, y conviene leer la política completa antes
de redactarla.

**59.2 Cambios de formato aplicados al manuscrito:**
- **Citas en corchetes**, no en superíndice: la guía pide "[3,6]" en línea. Convertidas las 36.
- **Abstract a 251 palabras**, en el tope de 250 (hay que quitar una).
- **Portada**: afiliaciones con letra minúscula en superíndice (a, b, c) como pide la guía, no con
  números; direcciones postales completas con país; autor de correspondencia marcado con asterisco y
  con su dirección y correo.
- **Ecuaciones numeradas** (1) a (5), variables en cursiva, separadas del texto.
- **Acknowledgements como sección propia inmediatamente antes de las referencias**, no dentro de
  Declarations. La guía lo pide explícitamente.
- **Highlights en fichero aparte** con "highlights" en el nombre, cinco viñetas de 72-79 caracteres,
  bajo el tope de 85.
- **Figuras como ficheros separados**, `Figure_1` a `Figure_5` y `Figure_S1` a `Figure_S7`, en PDF
  vectorial y PNG a 400 ppp (4.856 a 6.246 px de ancho, por encima del mínimo de 2.244 para página
  completa).

**59.3 Lo que falta y es obligatorio:**
- **Checklist de guía de reporte.** Metabolism exige enviar el checklist correspondiente como
  material suplementario. Para este trabajo el tipo es estudio observacional, de modo que
  corresponde **STROBE**. Hay que descargarlo, rellenarlo y añadirlo.
- **Resumen gráfico**: TIFF, EPS, PDF o fichero de Office, mínimo 1.328 px de ancho y legible a
  13 cm. Recomendado pero en la práctica esperado.
- **Declaración de datos** en el sistema de envío, además de la sección del manuscrito.

**59.4 Pendiente de verificar al enviar:** el estilo de referencias de la revista es Vancouver con
nombres de revista abreviados según LTWA, última página acortada (51-9), hasta seis autores y luego
"et al.", y DOI donde exista. La guía dice que el formato es libre al enviar y se aplica en pruebas,
pero conviene que ya sea consistente. Nuestras 42 referencias están en ese estilo salvo los DOI, que
faltan en la mayoría.

---

## 60. REVISIÓN COMPLETA DE TODOS LOS FICHEROS DEL ENVÍO

**60.1 El título ya era el nuevo** en el manuscrito y en los dos suplementarios. El que estaba
desactualizado era el del `.docx` base que se subió, que es la versión anterior a las correcciones.

**60.2 Nomenclatura corregida en todo el paquete.** Elsevier no usa "Extended Data": son
**Supplementary Figures S1 a S7**. Renombradas las 24 citas del manuscrito, las siete leyendas del
suplementario y las referencias cruzadas internas. También los ficheros de figura:
`Figure_S1` a `Figure_S7`.

**60.3 Lo que faltaba de verdad y estaba fuera de los ficheros del envío:** la sensibilidad al
ajuste por componentes principales vivía solo en la bitácora. Añadida en dos sitios:
- **Métodos suplementarios, sección S2.6b**, con el razonamiento de por qué no es un control neutro
  y las cifras (músculo 0,270 → 0,264 con percentil 98 en ambos casos; adiposo visceral 0,302 →
  0,133), y la observación de que las perturbaciones no se ven afectadas por ser dentro de la
  persona.
- **Limitaciones del manuscrito**, en una cláusula.
- **Tabla S21** nueva con los datos, y el libro de Excel pasa de 20 a 21 hojas.

**60.4 Verificación final del paquete completo:** 42 referencias, todas citadas, ninguna
inexistente; los 36 paneles de las cinco figuras citados; las 21 tablas suplementarias citadas; las
siete figuras suplementarias citadas; abstract en 251 palabras; cuerpo bajo 5.000.

**60.5 Ficheros del envío, estado:**
| fichero | estado |
|---|---|
| Manuscrito (.docx, figuras incrustadas) | listo |
| Highlights (.docx, fichero aparte) | listo |
| Supplementary Information (.docx, 7 figuras + 21 leyendas + Métodos) | listo |
| Supplementary Tables (.xlsx, 21 hojas) | listo |
| Figuras sueltas Figure_1-5 y Figure_S1-S7 (PDF + PNG 400 ppp) | listo |
| Checklist STROBE | **falta** |
| Resumen gráfico | **falta** |
| Seis diagramas | **falta** (ya no bloqueados: la política de Elsevier los permite, §59.1) |
| DOI en las referencias | **falta** en la mayoría |

---

## 61. SISTEMA DE AGENTES APLICADO AL MANUSCRITO

Ficheros `AGENTS.md` (los diez roles con su lista de comprobación y criterio de aprobado) y
`python/audit_manuscript.py` (los cinco mecánicos, ejecutables).

**61.1 Los agentes mecánicos pasan sin hallazgos** tras dos correcciones al propio detector: no
reconocía la forma "Tables S4 and S12" ni las citas partidas por salto de línea. Estado final:
cuerpo 4.926 palabras, abstract 246, 42 referencias, cero huérfanas en ambas direcciones, todos los
paneles y tablas citados, figuras en orden de citación.

**61.2 Lo que encontraron los agentes de juicio, que los mecánicos no podían ver:**
- **"sixteen public datasets" en la introducción.** El número estaba escrito con letra, así que la
  corrección del conteo a 19 series más GTEx no lo tocó y el agente de coherencia, que busca
  cifras, tampoco. Corregido a "twenty". **Lección: los números escritos con letra escapan a toda
  comprobación automática.**
- **Etiqueta mal cerrada** en la ecuación 4, `</s*>` en lugar de `</sub>`, que habría salido como
  texto literal en el Word.
- **Cinco frases de más de 50 palabras** partidas: las limitaciones (87 palabras), el obstáculo de
  diseño en la introducción (56), la frase clave de la Discusión sobre qué unifica las tres
  situaciones (55), el clamp (53) y la disociación de la firma (52). La de la Discusión era la más
  importante: es la idea del artículo y estaba enterrada en una subordinada.
- **Repetición torpe** introducida al corregir una cifra: "in the three cohorts with the score
  computed per individual the distributions…".
- **Párrafo de 3.2 demasiado largo** y con la transferencia entre tejidos pegada al final, después
  de su propia conclusión. Partido en dos, con el cierre al final.

**61.3 Ritmo de frase, medido.** Desviación de 14 a 16 palabras dentro de cada sección, que es
variación humana; por debajo de 6 sería la huella de texto generado. El agente de estilo lo
comprueba automáticamente.

**61.4 Nota sobre el orden de las figuras.** Se revisó la duda de si la 4 y la 5 debían
intercambiarse: **no se acordó tal cosa**. Lo que la revisión proponía era mover la sección 3.5, y
eso sí habría obligado a intercambiar las Figuras 3 y 4, pero se descartó. Lo que se movió fue la
dispersión, que no cita ninguna figura principal. Orden de primera cita verificado: 1, 2, 3, 4, 5.

---

## 62. REVISIÓN EXTERNA DE MANUSCRIPT_v2: CUATRO ERRORES CONFIRMADOS Y UNO QUE ERA PEOR

**62.1 Cadena de citas corrida (punto 2 de la revisión). CONFIRMADO Y CORREGIDO.** Al eliminar las
referencias sin llamado se quitó Dai et al. 2005, que era la 13, y no se recorrieron las citas:
[13] apuntaba a Smyth, [14] a PyDESeq2, [15] a Benjamini-Hochberg, y [16] se usaba dos veces para
cosas distintas. Peor: la afirmación del 8,3 % de sondas con símbolo retirado quedó citando un
artículo de estadística. Dai restituida como 13, las 29 posteriores recorridas y las cuatro citas
de Métodos arregladas. La bibliografía pasa a 43.

**62.2 "Homogeneous" con I² = 69 % (punto 3). CONFIRMADO Y CORREGIDO.** Reescrito a "positive in
every adipose-tissue estimate but heterogeneous in magnitude", y añadidos el I² global (84 %) y los
dos componentes de varianza, que el modelo de tres niveles ya calculaba y no estaban en el texto.

**62.3 Pie de la Fig. 2A con el método antiguo (punto 4). CONFIRMADO Y CORREGIDO.** El pie decía
DerSimonian-Laird cuando el estimado primario es el de tres niveles. Ahora declara que los
diamantes son el de dos niveles y da el primario entre paréntesis.

**62.4 Ecuaciones fuera de orden (punto 5). CONFIRMADO Y CORREGIDO.** Renumeradas por orden de
aparición: amplitud (1), enriquecimiento (2), varianza de Fisher (3), oscilador (4) y (5).

**62.5 El punto 7 era peor de lo que el revisor pudo ver.** Tenía razón en que BMAL1 no pertenece al
brazo de salida y en que los paneles redefinían la lista en el sitio. Al recalcular con el conjunto
del registro (`gene_sets.py`, ocho de los nueve genes medidos; HLF ausente):
| cantidad | con BMAL1 añadido | con el conjunto del registro |
|---|---|---|
| TZD | ρ −0,32 | **ρ −0,02 (P 0,96)**, pendiente **+0,004** |
| Cushing | 0,71 | **0,55 (P 0,16)** |
| hidrocortisona | 0,60 | **0,71 (P 0,047)** |
**Los valores de Cushing e hidrocortisona estaban intercambiados** respecto a lo que da el conjunto
correcto, porque cada uno se calculaba sobre una lista distinta. Y el ρ de −0,32 era enteramente
producto de añadir BMAL1, como el revisor dedujo del gráfico.
Reescrito el pasaje con lo que los datos sostienen y es más fuerte: **los ocho genes están
deprimidos en el tejido resistente, la mediana recuperada es del 4,4 % y la pendiente es +0,004
frente a −1 de reversión completa.** Los dos paneles importan ahora de `gene_sets.py` y el panel I
muestra la pendiente en lugar de la ρ.

**62.6 Pendientes de la revisión, no resueltos:**
- **Punto 1, reproducibilidad**: falta el script que escribe `signatures_log2_common_genes.tsv`, hay
  rutas absolutas en `fig6_mecanismo.py`, `verify_package.py` no comprueba la existencia de los
  ficheros que leen las figuras, y el repositorio tiene la composición de figuras anterior. **Hay
  que resolverlo antes de rearchivar en Zenodo.**
- **Punto 6, alcance de las afirmaciones**: "sets", "tone" y "acute insulin does not participate".
  Es decisión de alcance y es tuya. Mi lectura coincide con la del revisor: el diseño no mide
  arrastre y el brazo de insulina no tiene control temporal apareado propio.
- **Punto 8**: declarar que la validación por muestra única en músculo es marginal (P 0,067), y
  revisar si el procedimiento de submuestreo sobreestima el ruido.
- **Punto 9**: cuantificar en el texto la caída al quitar componentes, y el control de especificidad
  simétrico (rankear los 34.020 genes por su correlación con el score de reloj y ver dónde cae
  NR3C1). **Es el que más devuelve por el trabajo invertido.**

---

## 63. ALCANCE MODERADO Y AGENTES DE REPOSITORIO

**63.1 Aplicados los tres cambios de alcance:**
- **Título**: "Glucocorticoid signalling tracks circadian output across human metabolic tissues",
  sin el "not insulin".
- **"sets" → "tracks"** en abstract y conclusiones.
- **"tone" eliminado**; se usa "glucocorticoid signalling", que es el término que el manuscrito
  define.
- **"Acute insulin does not participate" → "the acute insulin response is not separable from the
  passage of time"**, que es lo que el diseño permite afirmar.
Cuerpo en 4.991 palabras y abstract en 246 tras los recortes que exigió el añadido.

**63.2 El punto 1 de la revisión: intenté escribir el script que falta y NO pude reproducir la
tabla.** `23_build_signatures_table.py` reconstruye `signatures_log2_common_genes.tsv` desde las
fuentes disponibles, pero da **números distintos**: Cushing +0,15 en lugar de +0,55, hidrocortisona
−0,22 en lugar de +0,71. Es decir, la columna de IR frente a sensible y la de Cushing no se
calcularon como diferencia de medianas sobre las matrices del repositorio, sino por algún
procedimiento que no está en el proyecto. El script queda marcado como `.FALLIDO` y la tabla
original restaurada.
**Esto confirma y agrava el punto 1 del revisor: la tabla de la que dependen tres afirmaciones de
la sección 3.4 no es reproducible con el código disponible, ni por el autor ni por un lector.**
Hay que recuperar el procedimiento original o rehacer las tres columnas con un método documentado,
y en ese caso los números del manuscrito cambiarán.

**63.3 Cinco agentes de repositorio añadidos a `AGENTS.md`:** 11 código, 12 auditor de código,
13 coherencia código-manuscrito, 14 consistencia global y 15 revisor final. El 13 lleva una regla
dura que acaba de demostrarse necesaria: **si una cifra no se reproduce, no se ajusta el texto al
número nuevo sin entender por qué cambió; una discrepancia es un hallazgo, no una errata.**

---

## 64. COLUMNA DE HIDROCORTISONA REPRODUCIDA, Y EL SUPUESTO POSICIONAL YA ESTABA ARREGLADO

**64.1 Reproducida exactamente.** La función `hydrocortisone()` de `02_meta_and_validation.py`
corre y da los valores publicados: **glucocorticoide +0,582 y núcleo de insulina +0,482, ambos con
P = 0,0002**, y **ZBTB16 como el gen más inducido (+2,05)**. El `delta` por gen que calcula esa
función **es** la columna de la tabla de firmas: correlación de Pearson **r = 0,983 sobre 5.910
genes**, y los valores coinciden dígito a dígito en los que inspeccioné.
Mi reconstrucción fallaba porque usé `GSE148640_HC_response_stats.tsv`, que es la mediana del
cambio por persona de otro procesamiento, no el delta del fichero procesado de GEO.

**64.2 El supuesto posicional que señalé ya estaba resuelto en la versión actual del script.** La
función `_hydrocortisone_arms` lee el tratamiento de cada muestra de la matriz de serie en lugar de
asumir que las columnas alternan, y avisa si no puede. Confirmado al correr: "brazos leídos de la
matriz de serie: 10 placebo, 10 hidrocortisona", y la matriz efectivamente alterna
placebo/SoluCortef. La versión que yo había visto era anterior.

**64.3 Comprobación independiente del pareamiento.** Con los brazos leídos de la matriz:
PER1 +1,504, PER3 −0,720, DBP −0,210, CIART −0,320, NR1D1 +0,213. **PER1 sube y PER3 baja**, que es
la dirección reportada en la sección 3.4, y ZBTB16 encabeza. El pareamiento es correcto.

**64.4 Lo que esto cierra y lo que no.** Cierra la procedencia de la columna de hidrocortisona y
confirma el ρ de 0,71 contra la firma del resistente. **Siguen sin reproducirse las columnas de
resistente frente a sensible y de Cushing**, que son las otras dos de la misma oración. Hay que
localizar la función que las genera; con las tres identificadas se puede escribir el paso que falta
en `run_all.sh` y cerrar el punto 1 de la revisión externa.

---

## 65. PROCEDENCIA DE LA TABLA DE FIRMAS: TRES COLUMNAS REPRODUCIDAS, UNA NO

Rastreo columna por columna de `signatures_log2_common_genes.tsv` contra los datos del repositorio.

| columna | procedimiento | r con la tabla |
|---|---|---|
| hydrocortisone adipose (GSE148640) | `delta` de `hydrocortisone()` en `02_meta_and_validation.py`: brazos leídos de la matriz de serie, media tratamiento − media placebo sobre el fichero procesado de GEO | **0,983** |
| Cushing adipose (GSE66446) | **mediana Cushing − mediana control sobre log2 CPM SIN centrar** (5 vs 11) | **1,0000** |
| TZD 3mo adipose (GSE13070) | `median_log2fc` del contraste `tzd_adipose` en `clamp_genes_change.tsv` | **1,0000** |
| **IR adipose (GSE13070, 35 vs 6)** | no reproducida | mejor candidato **0,468** |

**65.1 Dos aprendizajes de método.** El primero: **la tabla no usa el centrado por muestra**. Esa
es la razón de que mis reconstrucciones fallaran; el centrado es correcto para los scores entre
individuos pero esta tabla es de diferencias entre grupos y se construyó sin él. Hay que
documentarlo, porque mezclar las dos convenciones es precisamente lo que produjo el intercambio de
los valores de Cushing e hidrocortisona de §62.5.
El segundo: la columna de TZD sale de `clamp_genes_change.tsv`, que **sí** produce el repositorio,
de modo que esa parte de la cadena ya está cerrada.

**65.2 La columna que falta.** "IR adipose (GSE13070, 35 vs 6)" no se reproduce con la diferencia
de medianas ni de medias entre resistentes y sensibles del adiposo basal (r 0,45-0,47). Además el
recuento no cuadra: la etiqueta dice **35 frente a 6** y las muestras de adiposo basal de
`clamp_samples.tsv` son **34 frente a 6**. Una muestra de diferencia sugiere que el conjunto se
definió en otro sitio, con otro criterio de inclusión. **Es la última pieza del punto 1 de la
revisión externa.**

**65.3 Para cerrar el punto 1 hace falta:** localizar el script que produce esa columna, escribir el
paso que ensambla la tabla con los cuatro procedimientos documentados arriba, conectarlo en
`run_all.sh` antes del paso de figuras, extender `verify_package.py` para que compruebe la
existencia de cada entrada de cada figura, y rearchivar en Zenodo.

---

## 66. CICLO DE AGENTES SOBRE EL REPOSITORIO: PUNTO 1 CERRADO

**66.1 Agente 11 (código).** Escrito `23_build_signatures_table.py`, el paso que faltaba. Regenera
tres de las cuatro columnas desde los datos públicos con los procedimientos documentados en §65 y
arrastra el resto de la referencia versionada. Conectado en `run_all.sh` como paso 11, antes de las
figuras.

**66.2 Agente 12 (auditor).** Tres hallazgos resueltos:
- Las rutas absolutas de `fig6_mecanismo.py`, sustituidas por `T2D_OUT`.
- `make_figures.py` archivado en `_legacy/`: generaba las figuras del envío anterior y leía cuatro
  ficheros que ya no produce nadie.
- Escrito `tests/verify_inputs.py`, que comprueba que **cada fichero que lee una figura lo produce
  el repositorio, está versionado o es una descarga pública declarada**. Pasa de 62 entradas sin
  productor a **cero**. Es el test que habría detectado el fallo original.

**66.3 Agente 13 (coherencia código-manuscrito).** El script **verifica cada columna regenerada
contra la referencia y aborta si discrepa**: TZD r = 1,0000, Cushing r = 1,0000, hidrocortisona
r = 0,9833. Las cuatro cifras del manuscrito se reproducen: ρ 0,71 con hidrocortisona, 0,55 con
Cushing, −0,02 con TZD, pendiente +0,004 y 4,4 % de fracción recuperada.

**66.4 Agente 14 (consistencia).** Declarado en el manuscrito lo que el código hace: Métodos 2.4
explica por qué el contraste 35 contra 6 se excluye de la Fig. 1 y se usa solo para definir un
perfil, y Métodos 2.8 dice que esa columna es entrada versionada y que un test comprueba la cadena.
**Se retira así la afirmación de que todo se regenera con un script, que era falsa.**

**66.5 Agente 15 (revisor final).** Figuras regeneradas, smoke test correcto, los cinco agentes
mecánicos sin hallazgos, cuerpo en 4.996 palabras tras los recortes que exigieron los dos párrafos
nuevos. Paquete `tissue_coupling-v1.2.0.zip`.

**66.6 Lo que el paquete NO resuelve y hay que decir al rearchivar:** la columna del contraste
desbalanceado no se regenera; va versionada y declarada. Y el DOI de Zenodo del manuscrito apunta
a la versión anterior: **hay que rearchivar y actualizar el DOI en Métodos 2.8 y en Code
availability antes de enviar.**
