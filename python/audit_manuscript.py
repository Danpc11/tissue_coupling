#!/usr/bin/env python3
"""Agentes mecanicos del sistema de revision (AGENTS.md): coherencia, bibliografia, figuras y
contexto, estilo y editor general. Los que dependen de juicio quedan fuera.

  python python/audit_manuscript.py manuscript/MANUSCRIPT_v2.md
"""
import re, sys, collections

PLANTILLA = ["plays a crucial role", "sheds light on", "paves the way", "delve into",
             "a growing body of evidence", "it should be noted", "it is worth noting",
             "importantly,", "notably,", "crucially,", "in conclusion,",
             "has garnered", "underscores the", "highlights the importance"]
CONECTORES = ["moreover", "furthermore", "additionally", "in addition"]
VACIOS = ["greatly", "substantially", "markedly", "dramatically", "remarkably"]

def bib(s):
    """Agente 5: referencias y citas."""
    body = s[:s.index("## References")]; refs = s[s.index("## References"):]
    cited = set()
    for m in re.finditer(r'\[([0-9,\u2013]+)\]', body):
        for p in m.group(1).split(","):
            if "\u2013" in p:
                a, b = p.split("\u2013"); cited |= set(range(int(a), int(b) + 1))
            elif p.isdigit(): cited.add(int(p))
    total = set(int(m.group(1)) for m in re.finditer(r'^(\d+)\\?\. ', refs, re.M))
    out = []
    if total - cited: out.append(f"referencias sin citar: {sorted(total - cited)}")
    if cited - total: out.append(f"citas sin referencia: {sorted(cited - total)}")
    for bad in ("verify", "placeholder", "TODO", "XXXX"):
        if bad.lower() in refs.lower(): out.append(f"marcador de posicion en la bibliografia: '{bad}'")
    return len(total), out

def figuras(s):
    """Agente 4: paneles, tablas y figuras suplementarias citados."""
    body = re.sub(r'\n+', ' ', s[:s.index("## References")])
    out = []
    paneles = collections.defaultdict(set)
    for m in re.finditer(r'## Figure (\d) \|(.*?)(?=\n## |\Z)', s, re.S):
        paneles[m.group(1)] = set(re.findall(r'\*\*([A-I])\*\*', m.group(2)))
    for f, letras in sorted(paneles.items()):
        citadas = set()
        for m in re.finditer(r'Fig\. ' + f + r'([a-i](?:,[a-i])*)', body):
            citadas |= set(m.group(1).split(","))
        faltan = sorted({l.lower() for l in letras} - citadas)
        if faltan: out.append(f"Fig. {f}: paneles sin citar {faltan}")
    for tag, rng in (("Supplementary Fig", 8), ("Supplementary Table", 22)):
        # reconoce "Table S4", "Tables S4 and S12", "Tables S3, S9 and S18"
        presentes = set()
        for m in re.finditer(re.escape(tag) + r's? ((?:S\d+(?:,| and |\s)*)+)', body):
            presentes |= set(int(x) for x in re.findall(r'S(\d+)', m.group(1)))
        faltan = [i for i in range(1, rng) if i not in presentes]
        if faltan and len(faltan) < rng - 1: out.append(f"{tag}: sin citar {faltan}")
    orden = []
    r = body[body.index("## 3. Results"):]
    for m in re.finditer(r'Fig\. ([1-5])', r):
        if m.group(1) not in orden: orden.append(m.group(1))
    if orden != sorted(orden): out.append(f"figuras citadas fuera de orden: {orden}")
    return out

def estilo(s):
    """Agente 2: ingles de plantilla."""
    body = s[:s.index("## References")].lower()
    out = []
    for p in PLANTILLA:
        if p in body: out.append(f"frase de plantilla: '{p}'")
    for c in CONECTORES:
        n = len(re.findall(r'(?:^|\. )' + c, body))
        if n > 2: out.append(f"conector de relleno '{c}' al inicio de frase {n} veces")
    for v in VACIOS:
        if v in body: out.append(f"adverbio de enfasis vacio: '{v}'")
    # uniformidad de longitud de frase
    fr = [len(x.split()) for x in re.split(r'(?<=[.!?]) ', body) if 5 < len(x.split()) < 80]
    if fr:
        import statistics
        sd = statistics.pstdev(fr)
        if sd < 6: out.append(f"longitud de frase demasiado uniforme (sd={sd:.1f} palabras)")
    return out

def coherencia(s):
    """Agente 7: la misma cifra en todas partes."""
    out = []
    body = s[:s.index("## References")]
    sec = {}
    for name, a, b in (("abstract", "## Abstract", "## 1."), ("resultados", "## 3. Results", "## 4."),
                       ("discusion", "## 4. Discussion", "## Acknowledg")):
        try: sec[name] = body[body.index(a):body.index(b)]
        except ValueError: pass
    nums = lambda t: set(re.findall(r'(?:ρ = |= )([\-\u22120-9]\.\d{2})', t))
    for origen in ("abstract", "discusion"):
        if origen in sec and "resultados" in sec:
            faltan = nums(sec[origen]) - nums(sec["resultados"])
            if faltan: out.append(f"cifras en {origen} que no estan en resultados: {sorted(faltan)}")
    return out

def editor(s):
    """Agente 9: limites y piezas obligatorias."""
    out = []
    body = s[:s.index("## References")]
    cuerpo = sum(len(body[body.index(a):body.index(b)].split())
                 for a, b in (("## 1. Introduction", "## 2. Methods"), ("## 2. Methods", "## 3. Results"),
                              ("## 3. Results", "## 4. Discussion"), ("## 4. Discussion", "## Acknowledg")))
    ab = body[body.index("## Abstract"):body.index("**Keywords")]
    nab = len(re.sub(r'\*\*[^*]+\*\*', '', ab).replace("## Abstract", "").split())
    if cuerpo > 5000: out.append(f"cuerpo {cuerpo} palabras, tope 5000")
    if nab > 250: out.append(f"abstract {nab} palabras, tope 250")
    for pieza in ("CRediT", "Funding", "Conflicts of interest", "Data availability",
                  "Code availability", "artificial intelligence", "Acknowledgements"):
        if pieza.lower() not in s.lower(): out.append(f"falta la seccion/declaracion: {pieza}")
    return cuerpo, nab, out

if __name__ == "__main__":
    s = open(sys.argv[1] if len(sys.argv) > 1 else "manuscript/MANUSCRIPT_v2.md").read()
    nref, b = bib(s); f = figuras(s); e = estilo(s); c = coherencia(s); cuerpo, nab, ed = editor(s)
    print(f"cuerpo {cuerpo} palabras · abstract {nab} · referencias {nref}\n")
    fallos = 0
    for nombre, lista in (("5 bibliografia", b), ("4 figuras y contexto", f), ("2 estilo", e),
                          ("7 coherencia", c), ("9 editor general", ed)):
        print(f"agente {nombre}: {'OK' if not lista else str(len(lista)) + ' hallazgos'}")
        for x in lista: print("   -", x); 
        fallos += len(lista)
    print(f"\n{'SIN HALLAZGOS MECANICOS' if not fallos else str(fallos) + ' hallazgos a resolver'}")
