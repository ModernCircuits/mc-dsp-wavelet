#!/bin/sh

window=3

skip_start=5
skip_end=5

min=100
max=200

output="output/$min-$max"
mkdir -p $output

flags="--output $output --window $window --skip-start $skip_start --skip-end $skip_end --min $min --max $max"

python scripts/bpm.py $flags                 --filename ~/Downloads/BPM\ Detection/Asfalt\ Tango\ \(\~120\ bpm,\ 2:4\).wav
python scripts/bpm.py $flags                 --filename ~/Downloads/BPM\ Detection/Dead\ Limit\ \(172\ bpm,\ 4:4\).wav
python scripts/bpm.py $flags                 --filename ~/Downloads/BPM\ Detection/Grenade\ \(112\ bpm,\ 4:4\).wav
python scripts/bpm.py $flags                 --filename ~/Downloads/BPM\ Detection/Someone\ Like\ You\ \(135bpm,\ 4:4\).wav
python scripts/bpm.py $flags --beat "5/4"    --filename ~/Downloads/BPM\ Detection/Take\ Five\ \(\~90\ bpm,\ 5:4\).wav
python scripts/bpm.py $flags --beat "6/8"    --filename ~/Downloads/BPM\ Detection/Take\ it\ to\ the\ Limit\ \(\~91,\ 6:8\).wav
python scripts/bpm.py $flags                 --filename ~/Downloads/BPM\ Detection/The\ Rock\ Show\ \(192\ bpm,\ 4:4\).wav
python scripts/bpm.py $flags                 --filename ~/Downloads/BPM\ Detection/Whats\ Poppin\ \(145\ bpm,\ 4:4\).wav

qpdf --empty --pages output/"$min-$max"/*.pdf -- "output/$min-$max.pdf"


min=60
max=200

output="output/$min-$max"
mkdir -p $output

flags="--output $output --window $window --skip-start $skip_start --skip-end $skip_end --min $min --max $max"

python scripts/bpm.py $flags                 --filename ~/Downloads/BPM\ Detection/Asfalt\ Tango\ \(\~120\ bpm,\ 2:4\).wav
python scripts/bpm.py $flags                 --filename ~/Downloads/BPM\ Detection/Dead\ Limit\ \(172\ bpm,\ 4:4\).wav
python scripts/bpm.py $flags                 --filename ~/Downloads/BPM\ Detection/Grenade\ \(112\ bpm,\ 4:4\).wav
python scripts/bpm.py $flags                 --filename ~/Downloads/BPM\ Detection/Someone\ Like\ You\ \(135bpm,\ 4:4\).wav
python scripts/bpm.py $flags --beat "5/4"    --filename ~/Downloads/BPM\ Detection/Take\ Five\ \(\~90\ bpm,\ 5:4\).wav
python scripts/bpm.py $flags --beat "6/8"    --filename ~/Downloads/BPM\ Detection/Take\ it\ to\ the\ Limit\ \(\~91,\ 6:8\).wav
python scripts/bpm.py $flags                 --filename ~/Downloads/BPM\ Detection/The\ Rock\ Show\ \(192\ bpm,\ 4:4\).wav
python scripts/bpm.py $flags                 --filename ~/Downloads/BPM\ Detection/Whats\ Poppin\ \(145\ bpm,\ 4:4\).wav

qpdf --empty --pages output/"$min-$max"/*.pdf -- "output/$min-$max.pdf"


min=60
max=130

output="output/$min-$max"
mkdir -p $output

flags="--output $output --window $window --skip-start $skip_start --skip-end $skip_end --min $min --max $max"

python scripts/bpm.py $flags                 --filename ~/Downloads/BPM\ Detection/Asfalt\ Tango\ \(\~120\ bpm,\ 2:4\).wav
python scripts/bpm.py $flags                 --filename ~/Downloads/BPM\ Detection/Dead\ Limit\ \(172\ bpm,\ 4:4\).wav
python scripts/bpm.py $flags                 --filename ~/Downloads/BPM\ Detection/Grenade\ \(112\ bpm,\ 4:4\).wav
python scripts/bpm.py $flags                 --filename ~/Downloads/BPM\ Detection/Someone\ Like\ You\ \(135bpm,\ 4:4\).wav
python scripts/bpm.py $flags --beat "5/4"    --filename ~/Downloads/BPM\ Detection/Take\ Five\ \(\~90\ bpm,\ 5:4\).wav
python scripts/bpm.py $flags --beat "6/8"    --filename ~/Downloads/BPM\ Detection/Take\ it\ to\ the\ Limit\ \(\~91,\ 6:8\).wav
python scripts/bpm.py $flags                 --filename ~/Downloads/BPM\ Detection/The\ Rock\ Show\ \(192\ bpm,\ 4:4\).wav
python scripts/bpm.py $flags                 --filename ~/Downloads/BPM\ Detection/Whats\ Poppin\ \(145\ bpm,\ 4:4\).wav

qpdf --empty --pages output/"$min-$max"/*.pdf -- "output/$min-$max.pdf"