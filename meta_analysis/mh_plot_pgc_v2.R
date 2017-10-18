args <- commandArgs(trailingOnly = TRUE)
 scriptloc <- args[1]
 results <- args[2]
 outfile <- args[3]
 colorchoice <- args[4]
 highlight_p <- args[5] 
 highlightboundary <- args[5] 
 
 
 =5e-8,=50000

 blue=rgb(153,191,254,max=255)
 red=rgb(230,185,184,max=255)
 green=rgb(147,205, 221,max=255)
 purple=rgb(179,162,199,max=255)

library(data.table)

#read file
 dat1 <- fread(results, header=T,data.table=F)

#Load manhattanplot script
 print("Plotting data")
 source(scriptloc)

 attach(dat1)

    ManhattanPlot_AJS_cut(CHR, BP, P, SNP, genomesig = 5*(10^-8), genomesug = 0,photoname = '', 
    outname = outfile, colors = c(rgb(78,78,77,max=255),eval(parse(text = colorchoice))), sigcolor = 'darkred', sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = 'FALSE', pvalue_miss = 'NA', sigsnpcolor = "red",highlight_p=highlight_p, highlightboundary)
