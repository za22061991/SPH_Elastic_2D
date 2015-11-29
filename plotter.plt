graphName = sprintf("scr/scr%05i.png",i)
dataName = sprintf("dat/dat%05i.dat",i)
titleName = sprintf("Solid Bar n=%5i",i)
set output graphName
set title titleName
plot dataName w p pt 6 ps 2 notitle
i=i+1
if(i<max) reread