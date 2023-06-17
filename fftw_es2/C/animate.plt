set term png
set output "images/diffusivity1.png"
minT=0
maxT=0.9
set cbrange [minT:maxT]
plot 'data/diffusivity1.dat' matrix with image

set term png
set output "images/diffusivity2.png"
minT=0
maxT=0.9
set cbrange [minT:maxT]
plot 'data/diffusivity2.dat' matrix with image

set term png
set output "images/diffusivity3.png"
minT=0
maxT=0.9
set cbrange [minT:maxT]
plot 'data/diffusivity3.dat' matrix with image

set term png
set output "images/concentration_init.png"
minT=0
maxT=0.9
set cbrange [minT:maxT]
plot 'data/concentration_init.dat' matrix with image

set term gif animate delay 25
set output "images/animate.gif"
frames = 61
minT=0
maxT=0.1
set cbrange [minT:maxT]
plot 'data/concentration_init.dat' matrix with image
do for [i=1:frames] {
  plot 'data/concentration_'.i.'.dat' matrix  with image
}

