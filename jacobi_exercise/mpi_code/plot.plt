set term png
set output "initial.png"
unset colorbox
set palette rgb 33,13,10
set size square
plot 'initial.dat' with image

set term png
set output "result.png"
unset colorbox
set palette rgb 33,13,10
set size square
plot 'solution.dat' with image
