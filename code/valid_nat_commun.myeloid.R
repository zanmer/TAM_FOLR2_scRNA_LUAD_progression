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
celltype_cols = c("Alveolar Mac"="#CE5D2A", "Microglia/Mac"="#fabebe", "mo-Mac"="#C54A53", "Monocytes"="#B4DCAA", 
	"DCs"="#8C50A0", "Pleural Mac"="#E3B150", "Undetermined"="#888888")

celltype_new2_cols = c("ARM-RETN"="#3cb44b", "Mono-CD14"="#7570B3", "Mono-CD16"="#F781BF", "TAM-FOLR2"="#f58231")

celltypesub_cols = c("Activated DCs"="#B6327B", "Alveolar Mac"="#CE5D2A", "CD141+ DCs"="#5C509D", 
	"CD163+CD14+ DCs"="#FAE196", "CD1c+ DCs"="#4B86B8", "CD207+CD1a+ LCs"="#7EC0A6", 
	"Microglia/Mac"="#F1B16E", "mo-Mac"="#C54A53", "Monocytes"="#B4DCAA", "DCs"="#8C50A0",
	"pDCs"="#E887BD", "Pleural Mac"="#E3B150", "Undetermined"="#888888")
celltypesubs = names(celltypesub_cols)
celltypesub_cols["Undetermined"] = "#888888"

DC.cells = c("Activated DCs", "CD141+ DCs", "CD163+CD14+ DCs", "CD1c+ DCs", "CD207+CD1a+ LCs", "pDCs")
load("../data/Myeloid.dotplot.Rdata")
marker.data = data.matrix(marker.data)

meta = read.table("../data/NC.Myeloid.meta_data.txt",head=T,sep="\t")
meta$celltype = meta$Cell_subtype
meta$celltype[meta$Cell_subtype %in% DC.cells] = "DCs"
meta$seurat_clusters = meta$RNA_snn_res.0.8
meta$CD68 = marker.data["CD68",meta$barcode]
meta$MSR1 = marker.data["MSR1",meta$barcode]
meta$MRC1 = marker.data["MRC1",meta$barcode]
meta$C1QA = marker.data["C1QA",meta$barcode]
meta$APOE = marker.data["APOE",meta$barcode]
meta$FABP4 = marker.data["FABP4",meta$barcode]
meta$MCEMP1 = marker.data["MCEMP1",meta$barcode]
meta$MARCO = marker.data["MARCO",meta$barcode]
meta$RETN = marker.data["RETN",meta$barcode]
meta$FOLR2 = marker.data["FOLR2",meta$barcode]
meta$VCAN = marker.data["VCAN",meta$barcode]
meta$CD14 = marker.data["CD14",meta$barcode]
meta$FCGR3A = marker.data["FCGR3A",meta$barcode]
meta$CD1C = marker.data["CD1C",meta$barcode]

cell_type_defined = c(
  "0" = "ARM-RETN",
  "1" = "DCs",
  "2" = "TAM-FOLR2",
  "3" = "TAM-FOLR2",
  "4" = "Mono-CD14",
  "5" = "TAM",
  "6" = "ARM-RETN",
  "7" = "TAM",
  "8" = "TAM",
  "9" = "TRM",
  "10" = "ARM",
  "11" = "Undetermined",
  "12" = "Mono-CD16",
  "13" = "ARM",
  "14" = "ARM-RETN",
  "15" = "ARM",
  "16" = "TAM",
  "17" = "DCs",
  "18" = "TAM-FOLR2",
  "19" = "ARM",
  "20" = "TAM",
  "21" = "ARM",
  "22" = "Undetermined",
  "23" = "ARM",
  "24" = "Undetermined",
  "25" = "ARM-RETN",
  "26" = "ARM-RETN",
  "27" = "Undetermined"
)
meta$celltype_this = cell_type_defined[as.character(meta$seurat_clusters)]
meta$celltype_new = meta$celltype
meta$celltype_new[meta$celltype_this=="TAM-FOLR2"] = "TAM-FOLR2"
meta$celltype_new[meta$celltype_this=="ARM-RETN"] = "ARM-RETN"
meta$celltype_new[meta$celltype_this=="Mono-CD14"] = "Mono-CD14"
meta$celltype_new[meta$celltype_this=="Mono-CD16"] = "Mono-CD16"
meta$celltype_new2 = meta$celltype_new
meta$celltype_new2[!meta$celltype_new2 %in% c("TAM-FOLR2","ARM-RETN","Mono-CD14","Mono-CD16")] = NA


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
	guides(color=guide_legend(title="",nrow=5,byrow=TRUE,override.aes = list(size=1))) + 
	scale_color_manual(values= celltype_cols) + 
	xlab("UMAP 1") + ylab("UMAP 2") + 
	ggtitle("Cell types") + 
	theme_classic() + 
	theme(plot.title=element_text(hjust=0.5),
		legend.position=c(0.05,1),
		legend.background=element_blank(),
		legend.key.height=unit(2,"mm"),
		legend.key.width=unit(1,"mm"),
		legend.justification=c(0,1))
#p2 = add_cluster_labels(p2,meta,label="celltype",size=5,stat="median")	
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
p4 = meta %>%
	ggplot(aes(x=umap.1, y=umap.2, color=RETN)) + 
	geom_point_rast(size=0.1, scale=rast.scale, raster.dpi = getOption("ggrastr.default.dpi", 500)) + 
	scale_color_viridis(option="A") + 
	guides(color=guide_colorbar(title="")) + 
	xlab("UMAP 1") + ylab("UMAP 2") + 
	ggtitle("RETN (ARM-RETN)") + 
	theme_marker
p5 = meta %>%
	ggplot(aes(x=umap.1, y=umap.2, color= FOLR2)) + 
	geom_point_rast(size=0.1, scale=rast.scale, raster.dpi = getOption("ggrastr.default.dpi", 500)) + 
	scale_color_viridis(option="A") + 
	guides(color=guide_colorbar(title="")) + 
	xlab("UMAP 1") + ylab("UMAP 2") + 
	ggtitle("FOLR2 (TAM-FOLR2)") + 
	theme_marker
p6 = meta %>%
	ggplot(aes(x=umap.1, y=umap.2, color= CD14)) + 
	geom_point_rast(size=0.1, scale=rast.scale, raster.dpi = getOption("ggrastr.default.dpi", 500)) + 
	scale_color_viridis(option="A") + 
	guides(color=guide_colorbar(title="")) + 
	xlab("UMAP 1") + ylab("UMAP 2") + 
	ggtitle("CD14 (Mono-CD14)") + 
	theme_marker
p7 = meta %>%
	ggplot(aes(x=umap.1, y=umap.2, color= FCGR3A)) + 
	geom_point_rast(size=0.1, scale=rast.scale, raster.dpi = getOption("ggrastr.default.dpi", 500)) + 
	scale_color_viridis(option="A") + 
	guides(color=guide_colorbar(title="")) + 
	xlab("UMAP 1") + ylab("UMAP 2") + 
	ggtitle("FCGR3A (Mono-CD16)") + 
	theme_marker


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

RETN_nLung_test = test_celltype(bar_data, levels(bar_data$group), "nLung", celltype="ARM-RETN")
RETN_nLN_test = test_celltype(bar_data, levels(bar_data$group), "nLN", celltype="ARM-RETN")
FOLR2_nLung_test = test_celltype(bar_data, levels(bar_data$group), "nLung", celltype="TAM-FOLR2")
FOLR2_nLN_test = test_celltype(bar_data, levels(bar_data$group), "nLN", celltype="TAM-FOLR2")
CD14_nLung_test = test_celltype(bar_data, levels(bar_data$group), "nLung", celltype="Mono-CD14")
CD14_nLN_test = test_celltype(bar_data, levels(bar_data$group), "nLN", celltype="Mono-CD14")
CD16_nLung_test = test_celltype(bar_data, levels(bar_data$group), "nLung", celltype="Mono-CD16")
CD16_nLN_test = test_celltype(bar_data, levels(bar_data$group), "nLN", celltype="Mono-CD16")
nLung.y = 1.02
nLN.y = 0.95

celltype_new_cols = c(celltype_new2_cols, celltype_cols)
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
	filter(celltype_new == "ARM-RETN") %>%
	ggplot(aes(x=group, y=ratio)) + 
	geom_jitter(size=0.5) + 
	scale_fill_manual(values= origin_cols) + 
	geom_boxplot(aes(fill=group), outlier.size=0) + 
	guides(fill="none") + 
	xlab("") + ylab("Proportion") + 
	ggtitle("ARM-RETN") + 
	theme_classic() + 
	ylim(-0.01,1.05) +
	theme(plot.title=element_text(size=9),
		axis.title=element_text(size=7),
		axis.text.x=element_text(size=6,angle=60,hjust=1),
		axis.text.y=element_text(size=6)
	)
for(i in 1:nrow(RETN_nLung_test)) {
	p10 = p10 + annotate("text", x=i, y= nLung.y, color=origin_cols["nLung"], label=RETN_nLung_test$label[i], size=4)
	p10 = p10 + annotate("text", x=i, y= nLN.y, color=origin_cols["nLN"], label=RETN_nLN_test$label[i], size=4)
}
p11 = bar_data %>%
	filter(celltype_new == "TAM-FOLR2") %>%
	ggplot(aes(x=group, y=ratio)) + 
	geom_jitter(size=0.5) + 
	scale_fill_manual(values= origin_cols) + 
	geom_boxplot(aes(fill=group), outlier.size=0) + 
	guides(fill="none") + 
	xlab("") + ylab("") + 
	ggtitle("TAM-FOLR2") + 
	ylim(-0.01,1.05) +
	theme_classic() + 
	theme(plot.title=element_text(size=9),
		axis.title.x=element_text(size=7),
		axis.title.y=element_blank(),
		axis.text.x=element_text(size=6,angle=60,hjust=1),
		axis.text.y=element_text(size=6)
	)
for(i in 1:nrow(FOLR2_nLung_test)) {
	p11 = p11 + annotate("text", x=i, y= nLung.y, color=origin_cols["nLung"], label=FOLR2_nLung_test$label[i], size=4)
	p11 = p11 + annotate("text", x=i, y= nLN.y, color=origin_cols["nLN"], label=FOLR2_nLN_test$label[i], size=4)
}
p12 = bar_data %>%
	filter(celltype_new == "Mono-CD14") %>%
	ggplot(aes(x=group, y=ratio)) + 
	geom_jitter(size=0.5) + 
	scale_fill_manual(values= origin_cols) + 
	geom_boxplot(aes(fill=group), outlier.size=0) + 
	guides(fill="none") + 
	xlab("") + ylab("") + 
	ggtitle("Mono-CD14") + 
	ylim(-0.01,1.05) +
	theme_classic() + 
	theme(plot.title=element_text(size=9),
		axis.title.x=element_text(size=7),
		axis.title.y=element_blank(),
		axis.text.x=element_text(size=6,angle=60,hjust=1),
		axis.text.y=element_text(size=6)
	)
for(i in 1:nrow(CD14_nLung_test)) {
	p12 = p12 + annotate("text", x=i, y= nLung.y, color=origin_cols["nLung"], label=CD14_nLung_test$label[i], size=4)
	p12 = p12 + annotate("text", x=i, y= nLN.y, color=origin_cols["nLN"], label=CD14_nLung_test$label[i], size=4)
}
p13 = bar_data %>%
	filter(celltype_new == "Mono-CD16") %>%
	ggplot(aes(x=group, y=ratio)) + 
	geom_jitter(size=0.5) + 
	scale_fill_manual(values= origin_cols) + 
	geom_boxplot(aes(fill=group), outlier.size=0) + 
	guides(fill="none") + 
	xlab("") + ylab("") + 
	ggtitle("Mono-CD16") + 
	ylim(-0.01,1.05) +
	theme_classic() + 
	theme(plot.title=element_text(size=9),
		axis.title.x=element_text(size=7),
		axis.title.y=element_blank(),
		axis.text.x=element_text(size=6,angle=60,hjust=1),
		axis.text.y=element_text(size=6)
	)
for(i in 1:nrow(CD16_nLung_test)) {
	p13 = p13 + annotate("text", x=i, y= nLung.y, color=origin_cols["nLung"], label=CD16_nLung_test$label[i], size=4)
	p13 = p13 + annotate("text", x=i, y= nLN.y, color=origin_cols["nLN"], label=CD16_nLung_test$label[i], size=4)
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
	guides(color=guide_legend(title="",ncol=1,byrow=TRUE,override.aes = list(size=1))) + 
	theme(plot.title=element_text(size=10),
		legend.position="right",
		axis.text=element_text(size=7),
		axis.title=element_text(size=7),
		legend.key.height=unit(1,"mm"),
		legend.text=element_text(size=7),
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
	annotate("segment",x=5,xend=10,y=0,yend=0, color=celltype_cols["Monocytes"], lty="22") + 
	annotate("segment",x=5,xend=10,y=4,yend=4, color=celltype_cols["Monocytes"], lty="22") + 
	annotate("segment",x=5,xend=5,y=0,yend=4, color=celltype_cols["Monocytes"], lty="22") + 
	annotate("segment",x=10,xend=10,y=0,yend=4, color=celltype_cols["Monocytes"], lty="22") +
 	annotate("text",x=10,y=1,label="Monocytes",hjust=-0.1,size=2) +
	theme(
		axis.ticks.y=element_blank(),
		axis.text.y=element_blank(),
		axis.title.y=element_blank()
	)
p7 = p7 + 
	annotate("segment",x=5,xend=10,y=0,yend=0, color=celltype_cols["Monocytes"], lty="22") + 
	annotate("segment",x=5,xend=10,y=4,yend=4, color=celltype_cols["Monocytes"], lty="22") + 
	annotate("segment",x=5,xend=5,y=0,yend=4, color=celltype_cols["Monocytes"], lty="22") + 
	annotate("segment",x=10,xend=10,y=0,yend=4, color=celltype_cols["Monocytes"], lty="22") +
 	annotate("text",x=10,y=1,label="Monocytes",hjust=-0.07,size=2) +
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
	draw_plot(p2, x=0.3, y=0.78, width=0.4, height=0.2) +
	draw_plot(p3, x=0.7, y=0.78, width=0.4, height=0.2) +
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
	plot_annotation(title="Myeloid validation",theme=theme(title=element_text(face="bold",size=15,hjust=-0.1, vjust=-0.5)))
ggsave(ps1,filename="Myeloid.valid_nat_commun.pdf", paper="a4",pagecentre=FALSE, width=7,height=12.5)




