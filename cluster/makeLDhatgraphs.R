#!/usr/local/bioinfo/R/3.2.2/bin/Rscript --vanilla
# -*- coding: utf-8 -*-
# @author Lea Picard

library("optparse")

## define options
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="vcftools .hap.ld file", metavar="filename"),
  make_option(c("-o", "--out"), type="character", default="out", help="Output filename\n\t\t[default = %default.png]", metavar="filename")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

file = opt$file
out = opt$out

## program...
df = read.table(file, header=TRUE, na.strings="-nan", sep="\t", dec=".")
df1 = df[-1,]

png(paste(out, "_recombinationplots.png", sep=""), width = 1500, height = 1500, res=200)
plot(df1$Loci, df1$Mean_rho, type = "l", main="Mean Rho", xlab="Loci", ylab="Mean rho", ylim=c(0,12))
dev.off()
