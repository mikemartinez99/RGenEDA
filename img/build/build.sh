#!/usr/bin/env bash
# Rebuild RGenEDA graphics (logo, social preview, overview). Needs python3 + fontTools + Pillow,
# rsvg-convert (brew install librsvg), Google Chrome, and the K2D font (ships with macOS).
set -e
cd "$(dirname "$0")"
python3 build_logo.py ../logo.svg
rsvg-convert -w 1036 -h 1200 ../logo.svg -o ../logo.png
python3 build_social.py social-preview.svg
rsvg-convert -w 1280 -h 640 social-preview.svg -o ../social-preview.png
FIG="file://$(cd ../../docs/articles/Snail1_Vignette_files/figure-html && pwd)/"
sed -e "s#LOGO_PNG#file://$(cd .. && pwd)/logo.png#" -e "s#FIG/#$FIG#g" overview.html > overview_abs.html
"/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" --headless=new --disable-gpu --hide-scrollbars \
  --force-device-scale-factor=1 --window-size=3300,1875 --screenshot="$(pwd)/../overview.png" "file://$(pwd)/overview_abs.html"
rm -f overview_abs.html social-preview.svg
cp ../logo.png ../overview.png ../../man/figures/
python3 build_favicons.py
