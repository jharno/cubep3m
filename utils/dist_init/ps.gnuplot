#set term postscript enhanced color
#set output "power.ps"
set logscale x
set logscale y
plot "pk0.init" u 1:2 with lines
replot "pk.init" u 1:2 with lines
replot "pk.init" u 1:3 with lines
replot "pk.init" u 1:4 with lines
replot "pk.init" u 1:5 with lines
replot "pk.init" u 1:6 with lines
replot "pk.init" u 1:7 with lines
