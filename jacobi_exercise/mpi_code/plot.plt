set term png
set output "images/initial.png"
unset colorbox
set palette rgb 33,13,10
set size square
plot 'initial.dat' with image

set term png
set output "images/result.png"
unset colorbox
set palette rgb 33,13,10
set size square
plot 'solution.dat' with image
