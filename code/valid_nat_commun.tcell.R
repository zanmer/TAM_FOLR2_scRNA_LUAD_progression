rm(list=ls())
library(ggplot2)
library(reshape2)
library(foreach)
library(dplyr)
library(patchwork)
library(tidyverse)
library(ggsci)
library(grid)
library(ggrastr)
library(RColorBrewer)
library(ggalluvial)
library(viridis)
library(cowplot)
library(ggpubr)

rast.scale = 0.2

add_cluster_labels <- function(p, meta, label="seurat_clusters", size=2, stat="mean"){
	labels = as.character(unique(meta[,label]))
	for(x in labels){
		indices = which(meta[,label] == x)
		if(stat=="mean"){
			p = p + ggplot2::annotate("text", x=mean(meta$umap.1[indices]), y=mean(meta$umap.2[indices]), label=x, size=size)
		}else{
			p = p + ggplot2::annotate("text", x=median(meta$umap.1[indices]), y=median(meta$umap.2[indices]), label=x, size=size)
		}
	}
	return(p)
}

sample_cols = unique(c(brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3")))
defined_cols = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', sample_cols)


origin_cols = c("nLung"="#CA6627", "tLung"="#75A43A", "mBrain"="#456BAB", "mLN"="#D53E88", "nLN"="#8F6894", "tL/B"="#4C9C7A", "PE"="#DEAD3B")


celltypesub_cols = c("Naive CD4+ T"="#5E4FA2", "CD4+ Th"="#3288BD", "Exhausted Tfh"="#66C2A5", "Treg"="#bcf60c", "CD8+/CD4+ Mixed Th"="#FEE08B", "CD8 low T"="#FDAE61", "Exhausted CD8+ T"="#F46D43", "Cytotoxic CD8+ T"="#D53E4F", "Naive CD8+ T"="#9E0142", "NK"="#B1C6DD", "Undetermined"="#888888")

# "NK cells"="#B1C6DD", "T lymphocytes"="#4B7CB3", 
celltype_new2_cols = c("CD4-FOXP3"="#bcf60c", "CD4-NR4A3"="#4DAF4A", "CD8-GZMK"="#4363d8", "CD8-ZNF683"="#f032e6")

load("../data/Tcell.dotplot.Rdata")
marker.data = data.matrix(marker.data)

meta = read.table("../data/NC.Tcell.meta_data.txt",head=T,sep="\t")

meta$celltype = meta$Cell_subtype
meta$seurat_clusters = meta$RNA_snn_res.0.4

meta$CD3D = marker.data["CD3D",meta$barcode]
meta$CD3E = marker.data["CD3E",meta$barcode]
meta$PTPRC = marker.data["PTPRC",meta$barcode]
meta$CD4 = marker.data["CD4",meta$barcode]
meta$CD8A = marker.data["CD8A",meta$barcode]
meta$CD8B = marker.data["CD8B",meta$barcode]
meta$FOXP3 = marker.data["FOXP3",meta$barcode]
meta$NR4A3 = marker.data["NR4A3",meta$barcode]
meta$CCR7 = marker.data["CCR7",meta$barcode]
meta$CXCR6 = marker.data["CXCR6",meta$barcode]
meta$GZMK = marker.data["GZMK",meta$barcode]
meta$ZNF683 = marker.data["ZNF683",meta$barcode]
meta$CXCL13 = marker.data["CXCL13",meta$barcode]
meta$KLRD1 = marker.data["KLRD1",meta$barcode]
meta$CCL4L2 = marker.data["CCL4L2",meta$barcode]
meta$NKG7 = marker.data["NKG7",meta$barcode]


cell_type_defined = c(
  "0" = "CD4",
  "1" = "CD4",
  "2" = "CD8-GZMK",
  "3" = "NK",
  "4" = "CD4-FOXP3",
  "5" = "CD4",
  "6" = "CD8-ZNF683",
  "7" = "CD4-NR4A3",
  "8" = "CD8-GZMK",
  "9" = "CD4",
  "10" = "CD4",
  "11" = "CD4",
  "12" = "NK",
  "13" = "CD8-GZMK",
  "14" = "CD4",
  "15" = "Undetermined"
)

meta$celltype_this = cell_type_defined[as.character(meta$seurat_clusters)]
meta$celltype_new = meta$celltype
#meta$celltype_new[meta$celltype_this=="CD4-CCR7"] = "CD4-CCR7"
meta$celltype_new[meta$celltype_this=="CD4-FOXP3"] = "CD4-FOXP3"
meta$celltype_new[meta$celltype_this=="CD4-NR4A3"] = "CD4-NR4A3"
meta$celltype_new[meta$celltype_this=="CD8-GZMK"] = "CD8-GZMK"
meta$celltype_new[meta$celltype_this=="CD8-ZNF683"] = "CD8-ZNF683"
meta$celltype_new2 = meta$celltype_new
meta$celltype_new2[!meta$celltype_new2 %in% c("CD8-ZNF683","CD8-GZMK","CD4-NR4A3","CD4-FOXP3")] = NA



theme_marker = theme_classic() + 
	theme(plot.title=element_text(hjust=0.5,size=9),
		legend.position=c(0.03,0.98),
		legend.background=element_blank(),
		legend.direction="horizontal",
		legend.key.height=unit(3,"mm"),
		axis.text=element_text(size=7),
		axis.title=element_text(size=7),
		legend.key.width=unit(4,"mm"),
		legend.justification=c(0,1))	


p1 = meta %>%
	ggplot(aes(x=umap.1, y=umap.2, color=factor(seurat_clusters))) + 
	geom_point_rast(size=0.1, scale=rast.scale, raster.dpi = getOption("ggrastr.default.dpi", 500)) + 
	guides(color=guide_legend(title="",nrow=4,byrow=TRUE,override.aes = list(size=1))) + 
	xlab("UMAP 1") + ylab("UMAP 2") + 
	ggtitle("Clusters") + 
	theme_classic() + 
	theme(plot.title=element_text(hjust=0.5),
		legend.position=c(0.02,1.04),
		legend.background=element_blank(),
		legend.key.height=unit(2,"mm"),
		legend.key.width=unit(1,"mm"),
		legend.justification=c(0,1))
p2 = meta %>%
	ggplot(aes(x=umap.1, y=umap.2, color=celltype)) + 
	geom_point_rast(size=0.1, scale=rast.scale, raster.dpi = getOption("ggrastr.default.dpi", 500)) + 
	guides(color=guide_legend(title="",ncol=2,byrow=TRUE,override.aes = list(size=1))) + 
	scale_color_manual(values= celltypesub_cols) + 
	xlab("UMAP 1") + ylab("UMAP 2") + 
	ggtitle("Cell types") + 
	theme_classic() + 
	theme(plot.title=element_text(hjust=0.5),
		legend.position=c(0.05,0.4),
		legend.background=element_blank(),
		legend.key.height=unit(2,"mm"),
		legend.key.width=unit(1,"mm"),
		legend.justification=c(0,1))

p3 = meta %>%
	ggplot(aes(x=umap.1, y=umap.2, color=celltype_new2)) + 
	geom_point_rast(size=0.1, scale=rast.scale, raster.dpi = getOption("ggrastr.default.dpi", 500)) + 
	guides(color=guide_legend(title="",nrow=5,byrow=TRUE,override.aes = list(size=1))) + 
	scale_color_manual(values=celltype_new2_cols) + 
	xlab("UMAP 1") + ylab("UMAP 2") + 
	ggtitle("Cell subtypes") + 
	theme_classic() + 
	theme(plot.title=element_text(hjust=0.5),
		legend.position=c(0.05,1),
		legend.background=element_blank(),
		legend.key.height=unit(2,"mm"),
		legend.key.width=unit(1,"mm"),
		legend.justification=c(0,1))

meta.ZNF683 = meta[order(meta$ZNF683),]
p4 = meta.ZNF683 %>%
	ggplot(aes(x=umap.1, y=umap.2, color= ZNF683)) + 
	geom_point_rast(size=0.1, scale=rast.scale, raster.dpi = getOption("ggrastr.default.dpi", 500)) + 
	scale_color_viridis(option="A") + 
	guides(color=guide_colorbar(title="")) + 
	xlab("UMAP 1") + ylab("UMAP 2") + 
	ggtitle("ZNF683 (CD8-ZNF683)") + 
	theme_marker
meta.GZMK = meta[order(meta$GZMK),]
p5 = meta.GZMK %>%
	ggplot(aes(x=umap.1, y=umap.2, color= GZMK)) + 
	geom_point_rast(size=0.1, scale=rast.scale, raster.dpi = getOption("ggrastr.default.dpi", 500)) + 
	scale_color_viridis(option="A") + 
	guides(color=guide_colorbar(title="")) + 
	xlab("UMAP 1") + ylab("UMAP 2") + 
	ggtitle("CD8 (CD8-GZMK)") + 
	theme_marker
meta.NR4A3 = meta[order(meta$NR4A3),]
p6 = meta.NR4A3 %>%
	ggplot(aes(x=umap.1, y=umap.2, color=NR4A3)) + 
	geom_point_rast(size=0.1, scale=rast.scale, raster.dpi = getOption("ggrastr.default.dpi", 500)) + 
	scale_color_viridis(option="A") + 
	guides(color=guide_colorbar(title="")) + 
	xlab("UMAP 1") + ylab("UMAP 2") + 
	ggtitle("NR4A3 (CD4-NR4A3)") + 
	theme_marker
meta.FOXP3 = meta[order(meta$FOXP3),]
p7 = meta.FOXP3 %>%
	ggplot(aes(x=umap.1, y=umap.2, color=FOXP3)) + 
	geom_point_rast(size=0.1, scale=rast.scale, raster.dpi = getOption("ggrastr.default.dpi", 500)) + 
	scale_color_viridis(option="A") + 
	guides(color=guide_colorbar(title="")) + 
	xlab("UMAP 1") + ylab("UMAP 2") + 
	ggtitle("FOXP3 (CD4-FOXP3)") + 
	theme_marker



celltype_class = c("CD4"="CD4", "CD4-CCR7"="CD4", "CD4-CXCR6"="CD4", "CD4-FOXP3"="CD4", "CD4-NR4A3"="CD4", "CD8"="CD8", "CD8-CXCL13"="CD8", "CD8-GZMK"="CD8", "CD8-KLRD1"="CD8", "CD8-KLRD1-CCL4L2"="CD8", "CD8-ZNF683"="CD8", "NK"="Other", "T.cell"="Other")
pdata2 = subset(pdata, cell_type %in% names(celltype_class))
pdata2$class = celltype_class[pdata2$cell_type]

pdata.genes = unique(pdata2$features.plot)
pdata.clusters = unique(as.character(pdata2$id))
pmat = foreach(g=pdata.genes, .combine="rbind") %do% {
	foreach(cc=pdata.clusters, .combine="c") %do% {
		subset(pdata, features.plot==g & id==cc)$avg.exp
	}
}
colnames(pmat) = pdata.clusters
rownames(pmat) = pdata.genes


bar.data = data.frame(table(meta[,c("Sample_Origin", "seurat_clusters")])) 
bar.data$size = table(meta$seurat_clusters)[as.character(bar.data$seurat_clusters)]
bar.data$ratio = bar.data$Freq/bar.data$size


group_tmp = unique(meta[,c("Sample","Sample_Origin")])
sample_group = group_tmp[,2]
names(sample_group) = group_tmp[,1]
bar_data = data.frame(table(meta[,c("celltype_new","Sample")]))
bar_data = subset(bar_data, celltype_new!="Undetermined")
bar_data$group = factor(sample_group[as.character(bar_data$Sample)], levels=c("nLung","tLung","tL/B","nLN","mLN","PE","mBrain"))
bar_data$size = table(meta$Sample[meta$celltype_new!="Undetermined"])[bar_data$Sample]
bar_data$ratio = bar_data$Freq/bar_data$size
bar_data = bar_data[order(bar_data$Sample),]
bar_data = bar_data[order(bar_data$group),]
bar_data$Sample = factor(bar_data$Sample, levels=unique(bar_data$Sample))


test_celltype <- function(test.data, treated, control, celltype="TAM-FOLR2"){
	test.data = subset(test.data, celltype_new==celltype)
	pvalues = foreach(tt = treated, .combine="c") %do% {
		tmp = t.test(test.data$ratio[as.character(test.data$group)==tt], test.data$ratio[as.character(test.data$group)==control])
		tmp$p.value
	}
	test.res = data.frame(celltype=treated, pvalue=pvalues, label="")
	test.res$label[test.res$pvalue<0.05] = "*"
	test.res$label[test.res$pvalue<0.01] = "**"
	return(test.res)
}


ZNF683_nLung_test = test_celltype(bar_data, levels(bar_data$group), "nLung", celltype="CD8-ZNF683")
ZNF683_nLN_test = test_celltype(bar_data, levels(bar_data$group), "nLN", celltype="CD8-ZNF683")
GZMK_nLung_test = test_celltype(bar_data, levels(bar_data$group), "nLung", celltype="CD8-GZMK")
GZMK_nLN_test = test_celltype(bar_data, levels(bar_data$group), "nLN", celltype="CD8-GZMK")
NR4A3_nLung_test = test_celltype(bar_data, levels(bar_data$group), "nLung", celltype="CD4-NR4A3")
NR4A3_nLN_test = test_celltype(bar_data, levels(bar_data$group), "nLN", celltype="CD4-NR4A3")
FOXP3_nLung_test = test_celltype(bar_data, levels(bar_data$group), "nLung", celltype="CD4-FOXP3")
FOXP3_nLN_test = test_celltype(bar_data, levels(bar_data$group), "nLN", celltype="CD4-FOXP3")
nLung.y = 1.02
nLN.y = 0.95

celltype_new_cols = c(celltype_new2_cols, celltypesub_cols)
p8 = bar_data %>%
	ggplot(aes(x=Sample, y=ratio, fill=celltype_new)) + 
	geom_bar(position="stack", stat="identity") + 
	scale_fill_manual(values=celltype_new_cols) + 
	facet_grid(~group, scales="free_x", space="free") +
	xlab("") + ylab("Proportion") + 
	guides(fill="none") + 
	theme_classic() + 
	theme(
		axis.text.x=element_blank(),
		axis.title=element_text(size=8),
		axis.ticks.x=element_blank(),
		axis.title.x=element_blank(),
		strip.background=element_blank()
	)
p9 = bar_data %>%
	ggplot(aes(x=Sample, y=Freq, fill=celltype_new)) + 
	geom_bar(position="stack", stat="identity") + 
	scale_fill_manual(values=celltype_new_cols) + 
	facet_grid(~group, scales="free_x", space="free") +
	xlab("") + ylab("Number") + 
	guides(fill="none") + 
	theme_classic() + 
	theme(
		axis.text.x=element_text(angle=90,size=5,hjust=1),
		axis.title=element_text(size=8),
		strip.background=element_blank(),
		strip.text=element_blank()
	)
p10 = bar_data %>%
	filter(celltype_new == "CD8-ZNF683") %>%
	ggplot(aes(x=group, y=ratio)) + 
	geom_jitter(size=0.5) + 
	scale_fill_manual(values= origin_cols) + 
	geom_boxplot(aes(fill=group), outlier.colour="#FFFFFF") + 
	guides(fill="none") + 
	xlab("") + ylab("Proportion") + 
	ggtitle("CD8-ZNF683") + 
	theme_classic() + 
	ylim(-0.01,1.05) +
	theme(plot.title=element_text(size=9),
		axis.title=element_text(size=7),
		axis.text.x=element_text(size=6,angle=60,hjust=1),
		axis.text.y=element_text(size=6)
	)
for(i in 1:nrow(ZNF683_nLung_test)) {
	p10 = p10 + annotate("text", x=i, y= nLung.y, color=origin_cols["nLung"], label=ZNF683_nLung_test$label[i], size=4)
	p10 = p10 + annotate("text", x=i, y= nLN.y, color=origin_cols["nLN"], label=ZNF683_nLN_test$label[i], size=4)
}
p11 = bar_data %>%
	filter(celltype_new == "CD8-GZMK") %>%
	ggplot(aes(x=group, y=ratio)) + 
	geom_jitter(size=0.5) + 
	scale_fill_manual(values= origin_cols) + 
	geom_boxplot(aes(fill=group), outlier.colour="#FFFFFF") + 
	guides(fill="none") + 
	xlab("") + ylab("") + 
	ggtitle("CD8-GZMK") + 
	ylim(-0.01,1.05) +
	theme_classic() + 
	theme(plot.title=element_text(size=9),
		axis.title.x=element_text(size=7),
		axis.title.y=element_blank(),
		axis.text.x=element_text(size=6,angle=60,hjust=1),
		axis.text.y=element_text(size=6)
	)
for(i in 1:nrow(GZMK_nLung_test)) {
	p11 = p11 + annotate("text", x=i, y= nLung.y, color=origin_cols["nLung"], label=GZMK_nLung_test$label[i], size=4)
	p11 = p11 + annotate("text", x=i, y= nLN.y, color=origin_cols["nLN"], label=GZMK_nLN_test$label[i], size=4)
}
p12 = bar_data %>%
	filter(celltype_new == "CD4-NR4A3") %>%
	ggplot(aes(x=group, y=ratio)) + 
	geom_jitter(size=0.5) + 
	scale_fill_manual(values= origin_cols) + 
	geom_boxplot(aes(fill=group), outlier.colour="#FFFFFF") + 
	guides(fill="none") + 
	xlab("") + ylab("") + 
	ggtitle("CD4-NR4A3") + 
	ylim(-0.01,1.05) +
	theme_classic() + 
	theme(plot.title=element_text(size=9),
		axis.title.x=element_text(size=7),
		axis.title.y=element_blank(),
		axis.text.x=element_text(size=6,angle=60,hjust=1),
		axis.text.y=element_text(size=6)
	)
for(i in 1:nrow(NR4A3_nLung_test)) {
	p12 = p12 + annotate("text", x=i, y= nLung.y, color=origin_cols["nLung"], label=NR4A3_nLung_test$label[i], size=4)
	p12 = p12 + annotate("text", x=i, y= nLN.y, color=origin_cols["nLN"], label=NR4A3_nLung_test$label[i], size=4)
}
p13 = bar_data %>%
	filter(celltype_new == "CD4-FOXP3") %>%
	ggplot(aes(x=group, y=ratio)) + 
	geom_jitter(size=0.5) + 
	scale_fill_manual(values= origin_cols) + 
	geom_boxplot(aes(fill=group), outlier.colour="#FFFFFF") + 
	guides(fill="none") + 
	xlab("") + ylab("") + 
	ggtitle("CD4-FOXP3") + 
	ylim(-0.01,1.05) +
	theme_classic() + 
	theme(plot.title=element_text(size=9),
		axis.title.x=element_text(size=7),
		axis.title.y=element_blank(),
		axis.text.x=element_text(size=6,angle=60,hjust=1),
		axis.text.y=element_text(size=6)
	)
for(i in 1:nrow(FOXP3_nLung_test)) {
	p13 = p13 + annotate("text", x=i, y= nLung.y, color=origin_cols["nLung"], label=FOXP3_nLung_test$label[i], size=4)
	p13 = p13 + annotate("text", x=i, y= nLN.y, color=origin_cols["nLN"], label=FOXP3_nLung_test$label[i], size=4)
}
p14 = ggplot() + geom_blank() + xlim(0,1) + ylim(0.7,1) + theme_void() + 
	annotate("text", x=0.3, y=0.9, label="*", color=origin_cols["nLung"])  +
	annotate("text", x=0.3, y=0.9, label="p < 0.05, vs nLung", hjust=-0.1, vjust=0,size=2.5) + 
	annotate("text", x=0.3, y=0.8, label="**", color=origin_cols["nLung"])  +
	annotate("text", x=0.3, y=0.8, label="p < 0.01, vs nLung", hjust=-0.1, vjust=0,size=2.5) + 
	annotate("text", x=0.6, y=0.9, label="*", color=origin_cols["nLN"])  +
	annotate("text", x=0.6, y=0.9, label="p < 0.05, vs nLN", hjust=-0.1, vjust=0,size=2.5) + 
	annotate("text", x=0.6, y=0.8, label="**", color=origin_cols["nLN"])  +
	annotate("text", x=0.6, y=0.8, label="p < 0.01, vs nLN", hjust=-0.1, vjust=0,size=2.5) 



p1 = add_cluster_labels(p1,meta,size=2.5) + 
	guides(color="none") + 
	theme(plot.title=element_text(size=10),
		axis.text=element_text(size=7),
		axis.title=element_text(size=7)
	)
p1 = ggdraw(p1) + draw_text("A",x=0.02,y=0.98,size=18, fontface ="bold")
p2 = p2 + 
	guides(color=guide_legend(title="",ncol=1,byrow=TRUE,override.aes = list(size=0.8))) + 
	annotate("text", x=3, y=4, label="CD8+ T", hjust=0, size=2.5) + 
	annotate("text", x=8, y=5, label="CD4+ T", hjust=0, size=2.5) + 
	theme(plot.title=element_text(size=10),
		legend.position="right",
		axis.text=element_text(size=7),
		axis.title=element_text(size=7),
		legend.key.height=unit(0.1,"mm"),
		legend.text=element_text(size=5),
		legend.justification=c(0.5,0.5)
	)
p3 = p3 + 
	guides(color=guide_legend(title="",ncol=1,byrow=TRUE,override.aes = list(size=1))) + 
	theme(plot.title=element_text(size=10),
		legend.position="right",
		axis.text=element_text(size=7),
		axis.title=element_text(size=7),
		legend.key.height=unit(1,"mm"),
		legend.text=element_text(size=7),
		legend.justification=c(0.5,0.5)	
	)
p3 = ggdraw(p3) + draw_text("B",x=0.02,y=0.98,size=18, fontface ="bold")
p4 = p4 + 
	guides(color="none")
p4 = ggdraw(p4) + draw_text("C",x=0.02,y=0.98,size=18, fontface ="bold")
p5 = p5 + 
	guides(color="none") + 
	theme(
		axis.ticks.y=element_blank(),
		axis.text.y=element_blank(),
		axis.title.y=element_blank()
	)
p6 = p6 + 
	guides(color="none") + 
	theme(
		axis.ticks.y=element_blank(),
		axis.text.y=element_blank(),
		axis.title.y=element_blank()
	)
p7 = p7 + 
	theme(
		axis.ticks.y=element_blank(),
		axis.text.y=element_blank(),
		axis.title.y=element_blank(),
		legend.position="right",
		legend.direction="vertical",
		legend.key.height=unit(4,"mm"),
		legend.key.width=unit(3,"mm"),
	)
p8 = ggdraw(p8) + draw_text("D",x=0.01,y=0.98,size=18, fontface ="bold")
p10 = ggdraw(p10) + draw_text("E",x=0.02,y=0.98,size=18, fontface ="bold")

ps1 = ggdraw() + 
	draw_plot(p1, x=0.05, y=0.78, width=0.25, height=0.2) +
	draw_plot(p2, x=0.3, y=0.78, width=0.43, height=0.2) +
	draw_plot(p3, x=0.73, y=0.78, width=0.42, height=0.2) +
	draw_plot(p4, x=0.05, y=0.56, width=0.25, height=0.2) +
	draw_plot(p5, x=0.3, y=0.56, width=0.22, height=0.2) +
	draw_plot(p6, x=0.52, y=0.56, width=0.22, height=0.2) +
	draw_plot(p7, x=0.74, y=0.56, width=0.3, height=0.2) +
	draw_plot(p8, x=0.05, y=0.45, width=0.95, height=0.1) +
	draw_plot(p9, x=0.045, y=0.31, width=0.955, height=0.14) +
	draw_plot(p10, x=0.05, y=0.17, width=0.25, height=0.15) +
	draw_plot(p11, x=0.3, y=0.17, width=0.25, height=0.15) +
	draw_plot(p12, x=0.55, y=0.17, width=0.25, height=0.15) +
	draw_plot(p13, x=0.8, y=0.17, width=0.25, height=0.15) +
	draw_plot(p14, x=0.05, y=0.145, width=1, height=0.05) +
	plot_annotation(title="Tcell validation",theme=theme(title=element_text(face="bold",size=15,hjust=-0.1, vjust=-0.5)))

ggsave(ps1,filename="Tcel_nat.pdf", paper="a4",pagecentre=FALSE, width=7,height=12.5)

