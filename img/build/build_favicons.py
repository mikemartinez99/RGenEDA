"""Render favicons from ../logo.svg into ../../pkgdown/favicon and ../../docs."""
import os, subprocess, shutil, json
from PIL import Image

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(here, "..", ".."))
out = os.path.join(root, "pkgdown", "favicon")
os.makedirs(out, exist_ok=True)

# hex centred on a transparent square canvas
svg = open(os.path.join(here, "..", "logo.svg")).read()
svg = svg.replace('viewBox="0 0 518 600" width="518" height="600">', 'viewBox="-41 0 600 600" width="600" height="600">', 1)
fav_svg = os.path.join(out, "favicon.svg")
open(fav_svg, "w").write(svg)

for name, px in [("web-app-manifest-512x512.png", 512), ("web-app-manifest-192x192.png", 192),
                 ("apple-touch-icon.png", 180), ("favicon-96x96.png", 96)]:
    subprocess.run(["rsvg-convert", "-w", str(px), "-h", str(px), fav_svg, "-o", os.path.join(out, name)], check=True)

tmp = os.path.join(out, "_256.png")
subprocess.run(["rsvg-convert", "-w", "256", "-h", "256", fav_svg, "-o", tmp], check=True)
Image.open(tmp).save(os.path.join(out, "favicon.ico"), sizes=[(48, 48), (32, 32), (16, 16)])
os.remove(tmp)

json.dump({
    "name": "RGenEDA", "short_name": "RGenEDA",
    "icons": [
        {"src": "/web-app-manifest-192x192.png", "sizes": "192x192", "type": "image/png", "purpose": "maskable"},
        {"src": "/web-app-manifest-512x512.png", "sizes": "512x512", "type": "image/png", "purpose": "maskable"},
    ],
    "theme_color": "#ffffff", "background_color": "#ffffff", "display": "standalone",
}, open(os.path.join(out, "site.webmanifest"), "w"), indent=2)

# mirror into the built site (pkgdown copies these on build_site too)
docs = os.path.join(root, "docs")
for f in os.listdir(out):
    shutil.copy(os.path.join(out, f), os.path.join(docs, f))
shutil.copy(os.path.join(here, "..", "logo.png"), os.path.join(docs, "logo.png"))
print("favicons written to", out, "and docs/")
