#set xrange[-2:102]
#set yrange[-2:102]
plot '/tmp/best_tour.dat' u 2:3 w l, 'cities.dat' u 2:3 pt 7 ps 1, '/tmp/best_tour.dat' every ::0::0 u 2:3 pt 7 ps 1
pause 0.1
reread
