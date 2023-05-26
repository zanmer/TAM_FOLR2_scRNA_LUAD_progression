library(ggplot2)
library(reshape2)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)

## top10 differentail expressed genes of AT1/AT2 cells in single cell datasets
AT1.top10.genes = c('AGER', 'RTKN2', 'CAV1', 'EMP2', 'S100A4', 'MYL9', 'SPOCK2', 'ADIRF', 'VEGFA', 'AHNAK')
AT2.top10.genes = c('NAPSA', 'NPC2', 'SFTPA1', 'C4BPA', 'SFTPB', 'SERPINA1', 'PGC', 'CTSH', 'AQP3', 'SLC34A2')
epi.signatures = list(AT1_sig = AT1.top10.genes, AT2_sig = AT2.top10.genes)

## TCGA LUAD bulk RNAseq gene expression data
luad.expr = read.table("../data/luad_tcga.rnaseq.txt",head=T,sep="\t",check.names=F,stringsAsFactors=F)
luad.expr = subset(luad.expr, Hugo_Symbol!="")
rownames(luad.expr) = make.unique(luad.expr$Hugo_Symbol)
luad.expr = luad.expr[,-c(1,2)]
luad.scale = t(scale(t(luad.expr)))

## run GSVA method
luad_signature = gsva(luad.scale, 
	gset.idx.list = epi.signatures,
	method="gsva")
luad_signature = data.frame(t(luad_signature))


## read the clinical data and classify TCGA LUAD into two groups (AT1 high and low or AT2 high and low, using 75% quantile)
luad_meta = read.table("../data/luad_tcga.clinical.txt",head=T,sep="\t",stringsAsFactors=F,row.names=2,na.strings="[Not Available]")
luad.patient = substr(rownames(luad_signature),1,12)
luad_os = data.frame(
	luad_signature,
	pathologic_stage=luad_meta[luad.patient,"AJCC_PATHOLOGIC_TUMOR_STAGE"], 
	os_status=luad_meta[luad.patient,"OS_STATUS"], 
	os_month=luad_meta[luad.patient,"OS_MONTHS"]
)
luad_os$status = 0
luad_os$status[luad_os$os_status=="DECEASED"] = 1
luad_os$stage = NA
luad_os$stage[luad_os$pathologic_stage %in% c("Stage I", "Stage IA", "Stage IB", "Stage II", "Stage IIA", "Stage IIB")] = "early"
luad_os$stage[luad_os$pathologic_stage %in% c("Stage IIIA", "Stage IIIB", "Stage IV")] = "advanced"

luad_os$AT1_group = "low"
luad_os$AT1_group[luad_os$AT1_sig > quantile(luad_os$AT1_sig,0.75)] = "high"
luad_os$AT2_group = "low"
luad_os$AT2_group[luad_os$AT2_sig > quantile(luad_os$AT2_sig,0.75)] = "high"

## survival analysis and perform log rank test between two groups
fit.AT1 <- survfit(Surv(os_month, status) ~ AT1_group, data = luad_os)
fit.AT2 <- survfit(Surv(os_month, status) ~ AT2_group, data = luad_os)

## survival curve plottings
ggsurvplot(fit.AT1, conf.int=T, risk.table=T, pval=T, legend.title="", xlab ='Time in month', ggtheme =theme_light())
ggsurvplot(fit.AT2, conf.int=T, risk.table=T, pval=T, legend.title="", xlab ='Time in month', ggtheme =theme_light())


