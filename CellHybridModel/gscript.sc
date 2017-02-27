set term png size 1920,1080 enhanced font 'Verdana,24'
set output "short_reports.png"
set grid
set key outside
set termoption lw 5
set logscale y
set xlabel 'Time (min)'
set ylabel 'Count (cells/particles)'
plot "short_reports.txt" using 1:2 title 'CD4+ cells' with lines, \
"short_reports.txt" using 1:3 title 'CD4+ infected cells' with lines, \
"short_reports.txt" using 1:4 title 'CD8+ cells' with lines, \
"short_reports.txt" using 1:5 title 'Growth factor' with lines, \
"short_reports.txt" using 1:6 title 'Inflammatory factor' with lines, \
"short_reports.txt" using 1:7 title 'Free virus population' with lines
