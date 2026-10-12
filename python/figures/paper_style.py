"""Estilo y utilidades comunes a las figuras del articulo.

Dos decisiones que resuelven la mayor parte de los problemas de maquetacion:

1. `constrained_layout` en lugar de hspace/wspace manuales. Matplotlib calcula el espacio que
   necesita cada panel segun sus etiquetas reales, de modo que los solapes de titulo con etiqueta
   de panel, o de etiquetas de eje con el panel contiguo, desaparecen sin ajuste a mano. Es lo
   recomendado para figuras de revista y evita tener que recolocar cada vez que cambia un dato.

2. Un unico fichero de estilo, journal.mplstyle, para todas las figuras. Tipografia, grosores y
   paleta identicos en todos los paneles.

Paleta: Okabe-Ito para categorias (segura para daltonismo y razonable en escala de grises) y
cividis/RdBu para continuo. Se evita el rojo-verde, que es indistinguible para cerca del 8 % de
los varones, y las escalas tipo jet.
"""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt

MM = 1 / 25.4
WIDTH_SINGLE = 85 * MM      # columna simple
WIDTH_DOUBLE = 180 * MM     # doble columna

# Okabe-Ito con nombres utiles para este articulo
OI = dict(blue="#0072B2", vermillion="#D55E00", green="#009E73", purple="#CC79A7",
          orange="#E69F00", skyblue="#56B4E9", yellow="#F0E442", black="#000000", grey="#999999")
AXIS_COLORS = dict(Genetic="#000000", Sex="#CC79A7", HLA="#E69F00",
                   **{"Clock output": "#0072B2", "Glucocorticoid": "#D55E00",
                      "Insulin core": "#E69F00", "Core oscillator": "#56B4E9", "Oxidative": "#999999"})

def use_style():
    plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), "journal.mplstyle"))

def figure(width=WIDTH_DOUBLE, height=None, **kw):
    """Figura con constrained_layout activo y ancho en unidades de revista."""
    return plt.figure(figsize=(width, height or width * 0.62), layout="constrained", **kw)

def panel_label(ax, letter, dx=-0.02, dy=0.02):
    """Etiqueta de panel segun el estilo de Elsevier: MAYUSCULA, Arial, 10.5 pt, negrita,
    en la esquina superior izquierda. Los titulos de panel NO van dentro de la figura: el pie
    lleva un titulo breve y la descripcion de cada panel."""
    ax.text(dx, 1 + dy, str(letter).upper(), transform=ax.transAxes, fontsize=10.5,
            fontweight="bold", family="Liberation Sans", va="bottom", ha="right")

def stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
