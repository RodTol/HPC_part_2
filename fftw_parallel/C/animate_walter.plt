#### Initial diffusion ###############################
set term png
set output 'diffusivity.png'
minT = 0
maxT = 0.9
set cbrange [minT:maxT]
plot 'diffusivity.dat' matrix with image

#### Initial concentration ###########################
set term png
set output 'initial_concentration.png'
minT = 0
maxT = 0.9
set cbrange [minT:maxT]
plot 'initial_concentration.dat' matrix with image

#### Concentration evolving in time ###################
set term gif animate
unset key
set title 'Concentration in time'
set output 'animation.gif'
frames = 17
minT = 0
maxT = 0.9
set cbrange [minT: maxT]

do for [i = 1:frames] {
  plot 'concentration_'.i.'.dat' matrix  with image
}