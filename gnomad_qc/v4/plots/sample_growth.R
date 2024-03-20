# R script to get the sample growth of gnomAD v2, v3, and v4 exomes and genomes.
# To run this script, simply copy and paste the code into an R console or use `Rscript sample_growth.R` in the terminal.

# change the working directory as needed, note: the folder should already exist.
setwd("~/Documents")

if (!require("tidyr")) install.packages("tidyr")

color_amr = k_amr = '#ED1E24'
color_eur = k_eur = '#6AA5CD'
color_afr = k_afr = '#941494'
color_sas = k_sas = '#FF9912'
color_eas = k_eas = '#108C44'
color_oth = k_oth = '#ABB9B9'
color_mde = k_mde = '#000080'
color_asj = k_asj = 'gold2'

exac.colors <- c(color_eas, color_sas, color_eur, color_mde, color_afr, color_amr, color_asj, color_oth)

ancestries <- c("East Asian", "South Asian", "European", "Middle Eastern", "African", "Admixed American", "Ashkenazi Jewish", "Remaining")
exac.pop <- c(4327, 8256, (3307+33370), 0, 5203, 5789, 0, 454)
# v2 includes exomes and genomes
gnomadv2.pop <- c(9977, 15308, (12562+64603), 0, 12487, 17720, 5185, 3614)
gnomadv3.pop <- c(2604, 2419, 5316+34029, 158, 20744, 7647, 1736, 456+1047)
gnomadv4.pop <- c(19850, 43129, 556006+26710, 2884, 16740, 22362, 13068, 30198)

# plot 1: ExAC, v3, v2 ordered by sample size
pdf('fordist_gnomad_datasets_prev4_sizeorder_2023_09_20.pdf', height = 5, width=8)
barplot(as.matrix(cbind(exac.pop, gnomadv3.pop, gnomadv2.pop)),
        col=exac.colors,
        names.arg=c("ExAC", "gnomADv3", "gnomADv2"),
        border=NA, beside=FALSE,  main="", axes=FALSE, ylim=c(0,145000))
axis(side=2, at=seq(0,140000,10000), labels=formatC(seq(0,140000,10000),big.mark=',',format='fg'), las=2, cex.axis=1, lwd.ticks=1, lwd=0)
par(las=0)
legend('topleft', bty="n", ncol=1, cex=0.9,
        legend=rev(ancestries), fill=rev(exac.colors),
        border=rev(exac.colors), text.col=rev(exac.colors))
dev.off()

# plot 2: ExAC, v3, v2, v4 exomes ordered by sample size
pdf('fordist_gnomad_datasets_v4exomes_sizeorder_2023_09_20.pdf', height = 5, width=8)
barplot(as.matrix(cbind(exac.pop, gnomadv3.pop, gnomadv2.pop, gnomadv4.pop)),
        col=exac.colors,
        names.arg=c("ExAC", "gnomAD v3", "gnomAD v2", "gnomAD v4"),
        border=NA, beside=FALSE,  main="", axes=FALSE, ylim=c(0,750000))
axis(side=2, at=seq(0,750000,100000), labels=formatC(seq(0,750000,100000),big.mark=',',format='fg'), las=2, cex.axis=1, lwd.ticks=1, lwd=0)
par(las=0)
mtext("individuals", side=2, line=4)
legend('topleft', bty="n", ncol=1, cex=0.9,
        legend=rev(ancestries), fill=rev(exac.colors),
        border=rev(exac.colors), text.col=rev(exac.colors))
dev.off()

# plot 3: ExAC, v3, v2, v4 exomes+genomes ordered by sample size
pdf('fordist_gnomad_datasets_v4exomegenome_sizeorder_2023_09_20.pdf', height = 5, width=8)
barplot(as.matrix(cbind(exac.pop, gnomadv3.pop, gnomadv2.pop, gnomadv4.pop+gnomadv3.pop)),
        col=exac.colors,
        names.arg=c("ExAC", "gnomAD v3", "gnomAD v2", "gnomAD v4"),
        border=NA, beside=FALSE,  main="", axes=FALSE, ylim=c(0,820000))
axis(side=2, at=seq(0,800000,100000), labels=formatC(seq(0,800000,100000),big.mark=',',format='fg'), las=2, cex.axis=1, lwd.ticks=1, lwd=0)
par(las=0)
legend('topleft', bty="n", ncol=1, cex=0.9,
        legend=rev(ancestries), fill=rev(exac.colors),
        border=rev(exac.colors), text.col=rev(exac.colors))
dev.off()

# plot 4: ExAC, v2, v4 exomes+genomes ordered by sample size
pdf('fordist_gnomad_datasets_v4exomegenome_nov3_sizeorder_2023_09_20.pdf', height = 5, width=8)
barplot(as.matrix(cbind(exac.pop, gnomadv2.pop, gnomadv4.pop+gnomadv3.pop)),
        col=exac.colors,
        names.arg=c("ExAC", "gnomAD v2", "gnomAD v4"),
        border=NA, beside=FALSE,  main="", axes=FALSE, ylim=c(0,820000))
axis(side=2, at=seq(0,800000,100000), labels=formatC(seq(0,800000,100000),big.mark=',',format='fg'), las=2, cex.axis=1, lwd.ticks=1, lwd=0)
par(las=0)
legend('topleft', bty="n", ncol=1, cex=0.9,
        legend=rev(ancestries), fill=rev(exac.colors),
        border=rev(exac.colors), text.col=rev(exac.colors))
dev.off()
