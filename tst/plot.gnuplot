set term png

set output "comparison0.png"
c=2

set xlabel "Snapshot"
set ylabel "<T^2>"

p "proj-comparison-F.asc" u 1:c w l lw 5 title "TOTAL",\
"proj-recomparison3-F.asc" u 1:c w l lw 1 title "K=3",\
"proj-recomparison6-F.asc" u 1:c w l lw 2 title "K=6",\
"proj-recomparison9-F.asc" u 1:c w l lw 3 title "K=9",\
"proj-recomparison12-F.asc" u 1:c w l lw 4 title "K=12"


set output "comparison1.png"
c=3

set xlabel "Snapshot"
set ylabel "<U^2>"

p "proj-comparison-F.asc" u 1:c w l lw 5 title "TOTAL",\
"proj-recomparison3-F.asc" u 1:c w l lw 1 title "K=3",\
"proj-recomparison6-F.asc" u 1:c w l lw 2 title "K=6",\
"proj-recomparison9-F.asc" u 1:c w l lw 3 title "K=9",\
"proj-recomparison12-F.asc" u 1:c w l lw 4 title "K=12"


set output "comparison2.png"
c=4

set xlabel "Snapshot"
set ylabel "<V^2>"

p "proj-comparison-F.asc" u 1:c w l lw 5 title "TOTAL",\
"proj-recomparison3-F.asc" u 1:c w l lw 1 title "K=3",\
"proj-recomparison6-F.asc" u 1:c w l lw 2 title "K=6",\
"proj-recomparison9-F.asc" u 1:c w l lw 3 title "K=9",\
"proj-recomparison12-F.asc" u 1:c w l lw 4 title "K=12"

set output "comparison3.png"
c=8

set xlabel "Snapshot"
set ylabel "<(a*T)^2 + (b*U)^2 + (b*V)^2 >"

p "proj-comparison-F.asc" u 1:c w l lw 5 title "TOTAL",\
"proj-recomparison3-F.asc" u 1:c w l lw 1 title "K=3",\
"proj-recomparison6-F.asc" u 1:c w l lw 2 title "K=6",\
"proj-recomparison9-F.asc" u 1:c w l lw 3 title "K=9",\
"proj-recomparison12-F.asc" u 1:c w l lw 4 title "K=12"

