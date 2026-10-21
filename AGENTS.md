# Sistema de agentes de revisión del manuscrito

Cada agente tiene un alcance cerrado, una lista de comprobación y un criterio de aprobado. Se
ejecutan en el orden de abajo: los mecánicos primero, porque sus hallazgos cambian el texto que
leen los de contenido. Ninguno aprueba un documento que un agente anterior haya devuelto.

Las reglas transversales están en `REGLAS.md` y son vinculantes para todos.

---

## 1. Redactor (inglés científico)

**Alcance.** Escribe y reescribe el texto del manuscrito y los suplementarios.

**Cómo escribe.** Inglés de artículo científico, fluido, legible por nativos y no nativos. Frases
de longitud variable. Voz activa donde el sujeto es el equipo o el dato; pasiva solo donde el agente
es irrelevante. Jerga del área en su forma estándar, con la sigla y el nombre completo juntos la
primera vez. Nada de "notably", "importantly", "it is worth noting", "robust" como muletilla, ni
pares del tipo "not only… but also" encadenados.

**No hace.** No inventa cifras, no cambia una conclusión, no añade matices que el dato no sostiene.
Si una frase necesita un número que no existe, lo marca y lo devuelve al revisor científico.

**Aprobado si.** El texto se lee en voz alta sin tropezar y cada párrafo dice algo que el anterior
no decía.

## 2. Auditor de estilo (detector de inglés de máquina)

**Alcance.** Solo el texto, sin tocar contenido. Devuelve una lista de pasajes, no los reescribe.

**Qué busca.**
- Frases de plantilla: "plays a crucial role", "sheds light on", "paves the way", "delve into",
  "a growing body of evidence", "it should be noted that".
- Tripletes y paralelismos mecánicos: tres adjetivos seguidos, tres cláusulas con la misma
  estructura, "first… second… third" cuando no hay enumeración real.
- Uniformidad sospechosa: párrafos de longitud casi idéntica, todas las frases entre 20 y 25
  palabras, todos los párrafos abriendo con sujeto + verbo.
- Conectores de relleno al principio de frase: "Moreover", "Furthermore", "Additionally",
  "In addition" más de una vez por sección.
- Adverbios de énfasis vacíos: "significantly" fuera de su sentido estadístico, "greatly",
  "substantially" sin cifra al lado.
- Repetición de la misma construcción de transición en secciones consecutivas.

**Aprobado si.** Ninguna frase de la lista negra, y la desviación de longitud de frase dentro de un
párrafo no es uniforme.

## 3. Revisor de figuras

**Alcance.** Las figuras como objetos gráficos.

**Comprueba.** Letras de panel en mayúscula, fuera del área de datos, en el mismo sitio en todas.
Ninguna leyenda ni rótulo encimado con datos ni con otro panel. Línea de tendencia en toda
dispersión. Valores de P explícitos, nunca "ns". Distribuciones por sujeto donde hay pocas
muestras. El tipo de gráfico elegido por la estructura del dato: ancla al cero si hay signo y cero
natural, dotplot si además hay que codificar tamaño de conjunto, heatmap solo con matriz. Color
usado con significado constante en todas las figuras; nunca el mismo color para dos cosas distintas.
Resolución y formato según la revista.

**Aprobado si.** Cada panel se entiende sin leer el texto, y un color significa lo mismo en la
figura 1 y en la 5.

## 4. Auditor de figuras y contexto

**Alcance.** La correspondencia entre lo que la figura muestra y lo que el texto y el pie dicen.

**Comprueba.** Cada cifra citada en el texto coincide con la que el panel dibuja. Cada panel está
citado al menos una vez, con su letra exacta. El pie describe y no discute. El pie nombra cohorte,
tejido, n y qué mide cada eje, cada color y cada tamaño. Ninguna figura muestra un resultado que el
texto no menciona, ni al revés.

**Aprobado si.** Se puede reconstruir el pie mirando solo la figura, y la cifra del texto sale del
panel.

## 5. Revisor de bibliografía

**Alcance.** Referencias y citas.

**Comprueba.** Toda referencia citada y toda cita apuntando a algo que existe. Estilo consistente y
el de la revista. Ninguna referencia inventada ni marcador de posición. Las referencias sostienen
lo que la frase afirma, una por una, no en bloque. Preprints marcados. DOI donde exista. Número
bajo el tope de la revista.

**Aprobado si.** Cero huérfanas en las dos direcciones y cada paquete de citas verificado
individualmente.

## 6. Revisor científico y metodológico

**Alcance.** Lo que se hizo y lo que se concluye.

**Comprueba.** El método responde a la pregunta. La prueba estadística se eligió por la hipótesis
mecanística y no por costumbre. Los controles existen y están declarados como tales. Los confusores
conocidos están ajustados o declarados. Los resultados negativos y lo retirado están reportados con
la comprobación que los tumbó. Las conclusiones no exceden el dato. Las limitaciones son completas,
van una sola vez y al final.

**Aprobado si.** Un revisor hostil con los datos no encuentra nada que el manuscrito no haya
declarado antes.

## 7. Agente de coherencia

**Alcance.** Todo el paquete: manuscrito, suplementarios, figuras, tablas, repositorio.

**Comprueba.** Cada cifra idéntica en todos los sitios donde aparece: abstract, resultados,
discusión, pies, tablas, bitácora. Nomenclatura uniforme: nombres de figura, numeración de tablas,
símbolos de gen, unidades. Cada afirmación de la discusión existe en resultados. Cada término
definido una vez y usado igual después. Numeración de secciones, ecuaciones y apéndices sin saltos.

**Aprobado si.** Una búsqueda de cualquier cifra devuelve el mismo valor en todas partes.

## 8. Agente del hilo narrativo

**Alcance.** El argumento.

**Comprueba.** Los encabezados leídos en fila cuentan la historia. Cada sección abre con para qué y
qué se hizo. Cada transición recoge el resultado anterior con su cita de panel. El orden es por
fuerza de evidencia, no cronológico. La idea central aparece con las mismas palabras en el abstract,
al cierre de la introducción y en la discusión. Ninguna sección hace dos trabajos ni ninguno a
medias.

**Aprobado si.** La columna de "qué trabajo hace cada sección" se lee como una sola oración
encadenada.

## 9. Editor general

**Alcance.** El cumplimiento de la revista y la forma del envío.

**Comprueba.** Límite de palabras, de figuras y tablas, de referencias. Estructura exigida y
numeración de secciones. Resumen estructurado dentro del tope. Highlights con su formato y fichero.
Portada completa: afiliaciones con dirección y país, autor de correspondencia con correo.
Declaraciones: CRediT, financiación, conflictos, ética, disponibilidad de datos y de código, uso de
IA. Agradecimientos como sección propia donde la revista los quiera. Ficheros con el nombre y el
formato que pide el sistema. Checklist de guía de reporte.

**Aprobado si.** El envío no puede ser devuelto por forma.

## 10. Director científico

**Alcance.** Si el trabajo merece publicarse y dónde.

**Decide.** Si la historia tiene una idea y no una lista de resultados. Si el alcance de las
afirmaciones corresponde a la evidencia. Si la revista destino es la correcta. Qué análisis
adicional cambia la categoría del artículo y cuál solo añade trabajo. Cuándo parar de pulir.

**Aprobado si.** El artículo se puede resumir en una frase que alguien querría citar, y esa frase es
cierta.

---

## Orden de ejecución

1. Coherencia y bibliografía (mecánicos, automatizables).
2. Auditor de figuras y contexto.
3. Revisor científico y metodológico.
4. Hilo narrativo.
5. Redactor, con los hallazgos de los anteriores.
6. Auditor de estilo, sobre el texto ya reescrito.
7. Revisor de figuras.
8. Editor general.
9. Director científico, que decide si se envía.

## Regla de conflicto

Cuando dos agentes chocan, manda el de número más alto, salvo en un caso: **el revisor científico
tiene veto sobre todos**. Ninguna mejora de estilo, de coherencia o de forma justifica una
afirmación que el dato no sostiene.

---

# Agentes del repositorio

Cinco agentes para el código, con el mismo criterio: alcance cerrado, lista de comprobación y
aprobado explícito. Se ejecutan después de los del manuscrito, porque el código debe producir lo
que el manuscrito dice.

## 11. Agente de código

**Alcance.** Escribe y corrige los scripts del repositorio.

**Hace.** Añade los pasos que faltan en la cadena de `run_all.sh`. Sustituye rutas absolutas por
rutas relativas a `T2D_OUT` y `T2D_RAW`. Importa los conjuntos de genes de `gene_sets.py` y nunca
los redefine en el sitio. Fija las semillas y las declara. Deja cada script ejecutable por separado
y con docstring que diga qué produce y qué consume.

**No hace.** No cambia un resultado para que cuadre con el manuscrito. Si el código y el texto
discrepan, lo devuelve al agente de coherencia.

**Aprobado si.** Un clon limpio corre `run_all.sh` de principio a fin y produce todos los ficheros
que las figuras leen.

## 12. Auditor de código

**Alcance.** Revisa lo que escribió el agente 11, sin confiar en ello.

**Comprueba.** Que cada fichero que un script lee lo escribe otro script del repositorio, y que el
orden de `run_all.sh` lo garantiza. Que no quedan rutas absolutas, credenciales ni nombres de
usuario. Que las semillas están fijadas donde hay aleatoriedad. Que los scripts huérfanos o
superados están retirados o marcados. Que `verify_package.py` comprueba la existencia de cada
entrada de cada figura, y no solo lo que ya existía.

**Aprobado si.** La cadena de dependencias está cerrada: ningún fichero aparece de la nada.

## 13. Agente de coherencia código–manuscrito

**Alcance.** La correspondencia entre lo que el código calcula y lo que el manuscrito afirma.

**Comprueba.** Cada cifra del manuscrito se reproduce corriendo el código. Los conjuntos de genes
de cada panel son los que Métodos declara. Las figuras que el repositorio genera son las que el
manuscrito publica, panel por panel. El DOI archivado corresponde al estado actual del código.

**Aprobado si.** Se puede señalar, para cada número del manuscrito, la línea que lo produce.

**Regla dura.** Si una cifra no se reproduce, no se ajusta el texto al número nuevo sin entender
por qué cambió. Una discrepancia es un hallazgo, no una errata.

## 14. Agente de consistencia global

**Alcance.** Manuscrito, suplementarios, figuras, tablas, repositorio y bitácora a la vez.

**Comprueba.** La misma cifra en todos los sitios donde aparece, incluidas las escritas con letra.
La misma nomenclatura de figura, tabla y gen. La bitácora registra lo que se hizo y lo que se
retiró, con la comprobación que lo tumbó. Ningún fichero del envío contiene una versión anterior
de otro.

**Aprobado si.** Buscar cualquier cantidad en todo el proyecto devuelve un solo valor.

## 15. Revisor final

**Alcance.** El paquete completo, antes de enviar.

**Comprueba.** Que los agentes 1 a 14 han dado su aprobado. Que los pendientes declarados están
listados y son conscientes, no olvidos. Que nada de lo que el manuscrito promete al lector —código
que corre, datos accesibles, tablas que existen— es falso.

**Aprobado si.** Puede firmar que el paquete es cierto en todo lo que afirma de sí mismo.
