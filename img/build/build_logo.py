"""Build the RGenEDA hex logo SVG (text outlined to paths, K2D).

Concept: a magnifying glass over a muted expression matrix, revealing PCA structure.
"""
import math, random, sys
from textpath import text_elem

# palette
BG_TOP, BG_BOT = "#124443", "#071D20"
BORDER = "#31C46F"
CORAL, SKY, GOLD = "#F26D5B", "#5DB7E8", "#F2B035"
AXIS = "#9FB8B0"
WHITE, TAG = "#FFFFFF", "#A9C9BE"
INK = "#06181A"

HEX = "259,8 511,154 511,446 259,592 7,446 7,154"


def pca(x, y, w, h, scale=1.0):
    """PCA scatter with two sample groups inside a w x h box at (x, y)."""
    s = scale
    out = []
    out.append(f'<path d="M{x},{y} V{y+h} H{x+w}" fill="none" stroke="{AXIS}" stroke-opacity="0.55" '
               f'stroke-width="{2*s:.1f}" stroke-linecap="round"/>')
    cx1, cy1 = x + w * 0.33, y + h * 0.38
    cx2, cy2 = x + w * 0.68, y + h * 0.66
    rx, ry = w * 0.22, h * 0.2
    out.append(f'<ellipse cx="{cx1:.1f}" cy="{cy1:.1f}" rx="{rx:.1f}" ry="{ry:.1f}" transform="rotate(-18 {cx1:.1f} {cy1:.1f})" '
               f'fill="{CORAL}" fill-opacity="0.15" stroke="{CORAL}" stroke-opacity="0.65" stroke-width="{2*s:.1f}"/>')
    out.append(f'<ellipse cx="{cx2:.1f}" cy="{cy2:.1f}" rx="{rx:.1f}" ry="{ry:.1f}" transform="rotate(-14 {cx2:.1f} {cy2:.1f})" '
               f'fill="{SKY}" fill-opacity="0.15" stroke="{SKY}" stroke-opacity="0.65" stroke-width="{2*s:.1f}"/>')
    r = 5.5 * s
    pts1 = [(-0.12, 0.05), (-0.04, -0.07), (0.05, 0.05), (0.13, -0.04), (0.0, 0.11), (-0.08, -0.02)]
    pts2 = [(-0.12, 0.04), (-0.04, -0.06), (0.06, 0.06), (0.13, -0.04), (0.0, 0.11), (0.1, 0.1)]
    for dx, dy in pts1:
        out.append(f'<circle cx="{cx1+dx*w:.1f}" cy="{cy1+dy*h*1.3:.1f}" r="{r:.1f}" fill="{CORAL}" stroke="{INK}" stroke-width="{1.5*s:.1f}"/>')
    for dx, dy in pts2:
        out.append(f'<circle cx="{cx2+dx*w:.1f}" cy="{cy2+dy*h*1.3:.1f}" r="{r:.1f}" fill="{SKY}" stroke="{INK}" stroke-width="{1.5*s:.1f}"/>')
    return out


def logo_markup():
    """Inner markup of the logo in a 518x600 coordinate space (no <svg> wrapper)."""
    random.seed(9)
    parts = []
    A = parts.append
    A(f'''  <defs>
    <linearGradient id="bg" x1="0" y1="0" x2="0" y2="1">
      <stop offset="0" stop-color="{BG_TOP}"/>
      <stop offset="1" stop-color="{BG_BOT}"/>
    </linearGradient>
    <clipPath id="hexclip"><polygon points="{HEX}"/></clipPath>
    <clipPath id="lens"><circle cx="250" cy="240" r="84"/></clipPath>
  </defs>''')
    A(f'  <polygon points="{HEX}" fill="url(#bg)"/>')
    A('  <g clip-path="url(#hexclip)">')

    # ---- muted expression matrix (hex clip lets it run off the edges) ----
    greys = ["#2E5A58", "#3C6E6B", "#4F8682", "#1F4846", "#3C6E6B"]
    cell, gap = 16, 3
    for i in range(10):
        for j in range(20):
            x = 70 + j * (cell + gap)
            y = 150 + i * (cell + gap)
            roll = random.random()
            if roll < 0.07:
                fill, op = CORAL, 0.35
            elif roll < 0.14:
                fill, op = SKY, 0.35
            else:
                fill, op = random.choice(greys), 0.75
            A(f'    <rect x="{x}" y="{y}" width="{cell}" height="{cell}" rx="2" fill="{fill}" fill-opacity="{op}"/>')

    # ---- magnifying glass ----
    cx, cy, r = 250, 240, 88
    A(f'    <circle cx="{cx}" cy="{cy}" r="{r}" fill="{BG_TOP}"/>')
    A('    <g clip-path="url(#lens)">')
    for line in pca(cx - 62, cy - 58, 128, 112, scale=1.0):
        A('      ' + line)
    A('    </g>')
    A(f'    <circle cx="{cx}" cy="{cy}" r="{r}" fill="none" stroke="{WHITE}" stroke-width="7"/>')
    A(f'    <circle cx="{cx}" cy="{cy}" r="{r}" fill="none" stroke="{BORDER}" stroke-width="3"/>')
    ang = math.radians(40)
    hx1, hy1 = cx + (r + 3) * math.cos(ang), cy + (r + 3) * math.sin(ang)
    hx2, hy2 = cx + (r + 62) * math.cos(ang), cy + (r + 62) * math.sin(ang)
    A(f'    <line x1="{hx1:.1f}" y1="{hy1:.1f}" x2="{hx2:.1f}" y2="{hy2:.1f}" stroke="{WHITE}" stroke-width="18" stroke-linecap="round"/>')
    A(f'    <line x1="{hx1:.1f}" y1="{hy1:.1f}" x2="{hx2:.1f}" y2="{hy2:.1f}" stroke="{BORDER}" stroke-width="10" stroke-linecap="round"/>')
    A(f'    <path d="M{cx-58},{cy-40} A70,70 0 0 1 {cx-25},{cy-72}" fill="none" stroke="{WHITE}" stroke-opacity="0.45" stroke-width="4" stroke-linecap="round"/>')

    A('  </g>')  # end hex clip

    # ---- wordmark + tagline ----
    A('  ' + text_elem("RGenEDA", "Bold", 84, 259, 448, WHITE, anchor="middle", spacing=1))
    A('  ' + text_elem("genomic exploratory data analysis", "Regular", 16, 259, 485, TAG, anchor="middle", spacing=1.5))

    # ---- border ----
    A(f'  <polygon points="{HEX}" fill="none" stroke="{BORDER}" stroke-width="14" stroke-linejoin="round"/>')
    return parts


def logo_svg():
    return ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 518 600" width="518" height="600">\n'
            + "\n".join(logo_markup()) + "\n</svg>\n")


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "logo.svg"
    open(out, "w").write(logo_svg())
    print("wrote", out)
