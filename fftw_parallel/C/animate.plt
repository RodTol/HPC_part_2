set term png
set output "diffusivity.png"
minT=0
maxT=0.9
set cbrange [minT:maxT]
plot 'diffusivity.dat' matrix with image

set term png
set output "concentration_init.png"
minT=0
maxT=0.9
set cbrange [minT:maxT]
plot 'concentration_init.dat' matrix with image

set term gif animate delay 50
set output "animate.gif"
frames = 31
minT=0
maxT=0.1
set cbrange [minT:maxT]

plot 'concentration_init.dat' matrix with image

do for [i=1:frames] {
  plot 'concentration_'.i.'.dat' matrix  with image
}

