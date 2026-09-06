# Exploratory locus screen and same-model AKT2 measurements. Inputs: scripts46/47.
PROJ <- normalizePath(Sys.getenv("OVCAN_PROJ", unset=getwd()), mustWork=TRUE)
source(file.path(PROJ,"scripts/00b_figure_theme.R"))
suppressPackageStartupMessages(library(patchwork))
d <- read.csv(file.path(PROJ,"reports/molecular_extension_2026-09-06/copy_number/gene_model_evidence.csv"),
              na.strings=c("NA",""),check.names=FALSE)
genes <- c("CCNE1","AKT2","CCND1","MYC","MECOM","PRKCI","PIK3CA","KRAS","EGFR","ERBB2",
           "NF1","CDKN2A","CDKN2B","RB1","PTEN","BRCA1","BRCA2","RAD51D","CDK12","SMARCA2")
w <- d[d$wes_measured %in% c(TRUE,"True") & d$gene %in% genes,]
models <- unique(w[c("cell_line","patient_id","histotype_code")])
models <- models[order(match(models$histotype_code,c("HGS","CC","EC","LGS","MC","MMMT","SCCOHT")),
                       models$patient_id,models$cell_line),]
stopifnot(nrow(models)==23,nrow(w)==23*length(genes))
w$model <- factor(w$cell_line,levels=rev(models$cell_line))
w$gene <- factor(w$gene,levels=genes)
w$display_ratio <- pmax(-3,pmin(3,w$target_bin_median_log2))
w$dropout <- w$target_bins >= 3 & w$zero_depth_target_bins/w$target_bins >= .5
row_labels <- setNames(ifelse(models$cell_line==models$patient_id,models$cell_line,
                             paste0(models$cell_line,"  / ",models$patient_id)),models$cell_line)
a <- ggplot(w,aes(gene,model))+
  geom_tile(aes(fill=display_ratio),width=.95,height=.91)+
  geom_point(data=w[w$dropout,],aes(shape="At least half of bins have zero depth"),size=2,stroke=.6,colour=cook_ink)+
  scale_fill_gradient2(low=cook_slate,mid="white",high=cook_rust,midpoint=0,limits=c(-3,3),
                       breaks=c(-3,0,3),labels=c("≤−3","0","≥3"),name="Target-bin median\nlog2 ratio",na.value="#E2E8F0")+
  scale_shape_manual(values=4,name=NULL)+
  scale_y_discrete(labels=row_labels,expand=expansion(add=.3))+
  scale_x_discrete(expand=expansion(add=.3))+
  guides(fill=guide_colourbar(order=1,barwidth=unit(1.25,"in"),barheight=unit(.08,"in")),
         shape=guide_legend(order=2))+
  labs(title="A  Relative DNA profiles",x=NULL,y="Model / shared patient")+
  theme_ovcan(base_size=9)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=8),axis.text.y=element_text(size=8),
        axis.ticks=element_blank(),axis.line=element_blank(),axis.title.y=element_text(size=9),
        legend.position="bottom",legend.box="horizontal",legend.text=element_text(size=8),
        legend.title=element_text(size=8),legend.spacing.x=unit(5,"pt"),
        plot.title=element_text(size=11,face="bold"),plot.margin=margin(6,8,0,5))
k <- d[d$gene=="AKT2",]
kr <- k[!is.na(k$gene_segment_log2)&!is.na(k$rna_tpm),]
kp <- k[!is.na(k$rna_tpm)&!is.na(k$protein_log2_normalised),]
stopifnot(nrow(kr)==13,nrow(kp)==30)
style <- theme_ovcan(base_size=9)+theme(plot.title=element_text(size=11,face="bold"),
           plot.subtitle=element_text(size=8),plot.margin=margin(12,9,5,5))
b <- ggplot(kr,aes(gene_segment_log2,log2(rna_tpm+1)))+
  geom_point(colour=cook_slate,size=2)+
  geom_point(data=kr[kr$cell_line=="OV3331",],colour=cook_rust,size=2.5)+
  geom_text(data=kr[kr$cell_line=="OV3331",],aes(label=cell_line),hjust=1.12,vjust=-.55,size=2.8,colour=cook_rust)+
  scale_y_continuous(breaks=log2(c(64,256,1024)+1),labels=c("64","256","1,024"),expand=expansion(mult=c(.08,.2)))+
  scale_x_continuous(expand=expansion(mult=c(.08,.08)))+
  labs(title="B  AKT2 DNA and RNA",subtitle="13 models with both assays",
       x="Gene segment log2 ratio",y="RNA TPM (log scale)")+style
c <- ggplot(kp,aes(log2(rna_tpm+1),protein_log2_normalised))+
  geom_point(colour=cook_slate,size=2)+
  geom_point(data=kp[kp$cell_line=="OV3331",],colour=cook_rust,size=2.5)+
  geom_text(data=kp[kp$cell_line=="OV3331",],aes(label=cell_line),hjust=1.15,vjust=-.55,size=2.8,colour=cook_rust)+
  scale_x_continuous(breaks=log2(c(16,64,256,1024)+1),labels=c("16","64","256","1,024"),expand=expansion(mult=c(.07,.08)))+
  scale_y_continuous(expand=expansion(mult=c(.08,.2)))+
  labs(title="C  AKT2 RNA and protein",subtitle="30 models with both assays",
       x="RNA TPM (log scale)",y="Protein log2 normalised abundance")+style
p <- a/wrap_elements(full=b+c)+plot_layout(heights=c(2.15,1))
dest <- file.path(PROJ,"output/pdf/molecular_copy_number_extension.pdf")
save_fig(p,dest,w=7.6,h=8.4)
message("Exported molecular copy-number figure: 23 models / 20 selected loci; AKT2 paired assays.")
