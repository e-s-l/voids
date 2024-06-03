# usage: gnuplot -e "file='file.dat'" plot_data.gnuplot

# set output
set terminal pngcairo enhanced
set output 'test_plot.png'

# set labels
set title 'Test Plot'
set xlabel 'X'
set ylabel 'Y'
set xtic auto
set ytic auto

# set grid
set grid

# plot
plot file

#finish data write
unset output

