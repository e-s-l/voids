# usage: gnuplot -e "file='file.dat'" plot_data.gnuplot

# set output
set terminal pngcairo enhanced
set output 'test_plot.png'

# set labels
set title 'Test'
set xlabel 'R'
set ylabel 'Data'
set xtic auto
set ytic auto

# set grid
set grid

set style line 1 lt rgb "#9E77A1" lw 2

# plot
plot for [file in files] file with lines linestyle 1 title ''

#finish data write
unset output

