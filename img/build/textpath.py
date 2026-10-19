"""Outline K2D text to SVG path data with fontTools (no font dependency in the SVG)."""
from fontTools.ttLib import TTCollection
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen

K2D = "/System/Library/AssetsV2/PreinstalledAssetsV2/InstallWithOs/com_apple_MobileAsset_Font7/322738609e9863571011b800fb8f5e191be23ac5.asset/AssetData/K2D.ttc"
_coll = TTCollection(K2D)
_fonts = {f["name"].getDebugName(2): f for f in _coll.fonts}


def text_width(text, style, size, spacing=0.0):
    f = _fonts[style]
    upem = f["head"].unitsPerEm
    cmap = f.getBestCmap()
    hmtx = f["hmtx"]
    w = 0.0
    for i, ch in enumerate(text):
        g = cmap.get(ord(ch))
        if g is None:
            continue
        w += hmtx[g][0] * size / upem
        if i < len(text) - 1:
            w += spacing
    return w


def text_path(text, style, size, x, y, anchor="start", spacing=0.0):
    """Return an SVG path 'd' string for text rendered at (x, y baseline)."""
    f = _fonts[style]
    upem = f["head"].unitsPerEm
    scale = size / upem
    cmap = f.getBestCmap()
    gs = f.getGlyphSet()
    hmtx = f["hmtx"]
    width = text_width(text, style, size, spacing)
    if anchor == "middle":
        x -= width / 2
    elif anchor == "end":
        x -= width
    d = []
    cx = x
    for i, ch in enumerate(text):
        g = cmap.get(ord(ch))
        if g is None:
            continue
        pen = SVGPathPen(gs)
        tpen = TransformPen(pen, (scale, 0, 0, -scale, cx, y))
        gs[g].draw(tpen)
        seg = pen.getCommands()
        if seg:
            d.append(seg)
        cx += hmtx[g][0] * scale + (spacing if i < len(text) - 1 else 0)
    return " ".join(d)


def text_elem(text, style, size, x, y, fill, anchor="start", spacing=0.0, extra=""):
    return f'<path d="{text_path(text, style, size, x, y, anchor, spacing)}" fill="{fill}" {extra}/>'
