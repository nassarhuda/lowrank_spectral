using Plots
using StatPlots
numrows = 6
numcols = 6
T = [ x<y ? true : false for x in 1:numrows, y in 1:numcols][:]

# plots:
using Plots
font = Plots.font("sans-serif", 11)
pyplot(size=(400,400),guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)

nplots = div(7*6,2)

series1 = GS3_EigenAlign_lr[:][T]
series1p = GS3_EigenAlign[:][T]
series2 = GS3_Netalignbp[:][T]
series3 = GS3_Netalignmr[:][T]
series4 = GS3_Netalignbp_lr[:][T]
series5 = GS3_Netalignmr_lr[:][T]

y = hcat(series1,series1p,series2,series3,series4,series5)
f = violin([" " "  " "   " "    " "     " "      "],y,leg=false,ylim=(0,1))
ylabel!(f,"NCV-GS3 score")

savefig(f,"nvc218.pdf")

series1 = mharmonic1[:][T]
series1p = mharmonic_EA[:][T]
series2 = mharmonicBP_before[:][T]
series3 = mharmonicMR_before[:][T]
series4 = mharmonic3[:][T]
series5 = mharmonic2[:][T]

y = hcat(series1,series1p,series2,series3,series4,series5)
g = violin(["LR" "EA" "BP" "Klau" "LR+\nBP" "LR+\nKlau"],y,leg=false,ylim=(0,1))
ylabel!(g,"F-NC score")


l = @layout([a;b])
p = plot(
      f,
      g,
      leg=false,
      layout=l,
      border=:black)

savefig("GS3_NCV_Magna4.pdf")
