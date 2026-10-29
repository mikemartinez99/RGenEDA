"""Build the RGenEDA GitHub social preview (1280x640) as SVG with outlined K2D text."""
import sys
from textpath import text_elem, text_width
from build_logo import logo_markup, CORAL, SKY, GOLD, BORDER, INK

OUT = sys.argv[1] if len(sys.argv) > 1 else "social-preview.svg"
W, H = 1280, 640
SAGE = "#8FB0A6"
WHITE = "#FFFFFF"

parts = []
A = parts.append
A(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">')
A('''  <defs>
    <linearGradient id="sbg" x1="0" y1="0" x2="0" y2="1">
      <stop offset="0" stop-color="#153F40"/>
      <stop offset="1" stop-color="#071D20"/>
    </linearGradient>
    <radialGradient id="glow" cx="0.3" cy="0.5" r="0.7">
      <stop offset="0" stop-color="#1E5A58" stop-opacity="0.55"/>
      <stop offset="1" stop-color="#071D20" stop-opacity="0"/>
    </radialGradient>
  </defs>''')
A(f'  <rect width="{W}" height="{H}" fill="url(#sbg)"/>')
A(f'  <rect width="{W}" height="{H}" fill="url(#glow)"/>')


# ---- hex logo ------------------------------------------------------------------
A('  <g transform="translate(62 94) scale(0.745)">')
parts.extend(logo_markup())
A('  </g>')

# ---- headline -------------------------------------------------------------------
X0 = 470
A('  ' + text_elem("Streamlined, reproducible", "Bold", 48, X0, 252, WHITE))
A('  ' + text_elem("genomic EDA", "Bold", 48, X0, 314, WHITE))

# ---- pipeline pills ----------------------------------------------------------
steps = [("Normalized counts", BORDER), ("HVGs", GOLD), ("PCA", CORAL), ("Eigencorrelations", SKY), ("DEGs", BORDER)]
x = X0
y = 360
ph = 38
for i, (label, col) in enumerate(steps):
    tw = text_width(label, "Regular", 16)
    pw = tw + 30
    A(f'  <rect x="{x}" y="{y}" width="{pw:.1f}" height="{ph}" rx="{ph/2}" fill="{col}" fill-opacity="0.16" stroke="{col}" stroke-width="2"/>')
    A('  ' + text_elem(label, "Regular", 16, x + pw / 2, y + 25, WHITE, anchor="middle"))
    x += pw
    if i < len(steps) - 1:
        ax = x + 6
        A(f'  <path d="M{ax},{y+ph/2} H{ax+12} M{ax+8},{y+ph/2-4} L{ax+13},{y+ph/2} L{ax+8},{y+ph/2+4}" fill="none" stroke="{SAGE}" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/>')
        x += 26

# ---- url ---------------------------------------------------------------------------
A('  ' + text_elem("github.com/mikemartinez99/RGenEDA", "Regular", 17, X0, 446, SAGE, spacing=0.3))

A('</svg>')
open(OUT, "w").write("\n".join(parts) + "\n")
print("wrote", OUT)
