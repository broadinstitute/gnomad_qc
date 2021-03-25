library(umap)
library(ggplot2)
genomes=read.delim("/Users/mwilson/gnomad/v3.1/data/gnomad_v3.1_metadata_v3.1.tsv.gz", header=T)
genomes=genomes[genomes$release=="true",]

#Grab population PCs
ancestry=rbind(genomes[,(grep("PC", names(genomes)))])

# Project PCs onto UMAP with defaults
ancestry_umap=umap(ancestry[,c(1:6,8:16)], n_neighbors=20, min_dist=0.75)
save(ancestry_umap, file="/Users/mwilson/gnomad/v3.1/data/gnomad.ancestry.20nn.md0.75.UMAP.RData")

# Add coloring labels and plot
ancestry_umap[is.na(ancestry_umap)] <-"Missing"
color_amr = k_amr = '#ED1E24'
color_nfe = k_nfe = '#6AA5CD'
color_afr = k_afr = '#941494'
color_sas = k_sas = '#FF9912'
color_eas = k_eas = '#108C44'
color_oth = k_oth = '#ABB9B9'
color_mid = k_mid = '#33CC33'
color_asj = k_asj = 'coral'
color_fin = k_fin = '#002F6C'
pop_colors = c('afr' = color_afr,
               'amr' = color_amr,
               'eas' = color_eas,
               'fin' = color_fin,
               'nfe' = color_nfe,
               'oth' = color_oth,
               'sas' = color_sas,
               'mid' = color_mid,
               'asj' = color_asj,
               'ami' = 'black',
               'Missing' = 'gray'
)

colScale <- scale_color_manual(name = 'pops', values = pop_colors)

plot_ancestry=data.frame(ancestry_umap$layout)
plot_ancestry$pop=as.factor(as.character(genomes$population_inference.pop))
plot_ancestry$plot_pop=as.character(plot_ancestry$pop)
png("/Users/mwilson/gnomad/v3.1/plots/gnomad.ancestry.20nn.UMAP.RData", res=600, width=8, unit='in', height=8)
ggplot(plot_ancestry, aes(X1, X2)) + geom_point(aes(color=plot_pop), size=1, alpha=.5, shape=16) + colScale + theme(axis.line=element_blank(),
                                                                                                                    axis.text.x=element_blank(),
                                                                                                                    axis.text.y=element_blank(),
                                                                                                                    axis.ticks=element_blank(),
                                                                                                                    axis.title.x=element_blank(),
                                                                                                                    axis.title.y=element_blank(),
                                                                                                                    panel.background=element_blank(),
                                                                                                                    panel.border=element_blank(),
                                                                                                                    panel.grid.major=element_blank(),
                                                                                                                    panel.grid.minor=element_blank(),
                                                                                                                    plot.background=element_blank(),
                                                                                                                    legend.position = "none") # for legend, remove this line and add + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3)))
dev.off()

