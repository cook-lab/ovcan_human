# Exploratory model-selection figure; script 40 supplies all numerical inputs.
# Separate from the numbered manuscript figures. Narrative stays in the caption.
PROJ <- normalizePath(Sys.getenv("OVCAN_PROJ", unset=getwd()), mustWork=TRUE)
source(file.path(PROJ, "scripts/00b_figure_theme.R"))
suppressPackageStartupMessages(library(patchwork))
d <- read.csv(file.path(PROJ, "reports/clinical_classification_2026-09-06/ccne1_model_summary.csv"),
              check.names=FALSE, na.strings=c("NA", ""))
stopifnot(nrow(d)==42, sum(!is.na(d$gene_segment_log2c_auto))==23,
          sum(!is.na(d$rna_tpm))==31, sum(!is.na(d$protein_log2_normalised))==31)
hist_order <- c("HGS", "CC", "EC", "LGS", "MC", "MMMT", "SCCOHT")
d <- d[order(match(d$histotype_code,hist_order), is.na(d$gene_segment_log2c_auto),
             -d$gene_segment_log2c_auto, -d$rna_tpm, d$cell_line), ]
d$y <- rev(seq_len(nrow(d)))
guides_data <- d[d$y %% 2 == 0, ]
breaks <- d$y[c(FALSE,d$histotype_code[-1]!=d$histotype_code[-nrow(d)])]+.5
base <- function() {
  ggplot(d,aes(y=y)) +
    geom_rect(data=guides_data,aes(xmin=-Inf,xmax=Inf,ymin=y-.5,ymax=y+.5),
              inherit.aes=FALSE,fill="#F5F7F9",colour=NA) +
    geom_hline(yintercept=breaks,colour=cook_hair,linewidth=.35) +
    scale_y_continuous(limits=c(.5,42.5),expand=c(0,0),breaks=NULL) +
    theme_ovcan(base_size=9) +
    theme(axis.title.y=element_blank(),axis.line.y=element_blank(),
          plot.title=element_text(size=10,face="bold",margin=margin(b=8)),
          plot.margin=margin(5,5,5,5),legend.position="bottom")
}
labels <- base()+geom_text(aes(x=0,label=cell_line),hjust=0,size=2.85,colour=cook_ink)+
  geom_text(aes(x=1.19,label=histotype_code),hjust=0,size=2.3,colour=cook_ink_muted)+
  scale_x_continuous(limits=c(0,1.68),expand=c(0,0),breaks=NULL)+
  labs(title="Model / histotype",x=" ")+theme(axis.line.x=element_blank(),axis.ticks.x=element_blank())
cn <- base()+geom_vline(xintercept=0,colour=cook_hair,linewidth=.4)+
  geom_segment(aes(x=gene_target_median_log2c_auto,xend=gene_segment_log2c_auto,yend=y),
               colour=cook_hair,linewidth=.55,na.rm=TRUE)+
  geom_point(aes(x=gene_segment_log2c_auto,shape="Segment",colour="Segment"),size=2,na.rm=TRUE)+
  geom_point(aes(x=gene_target_median_log2c_auto,shape="Target-bin median",colour="Target-bin median"),
             size=1.7,na.rm=TRUE)+
  geom_text(data=d[is.na(d$gene_segment_log2c_auto),],aes(x=-.42,label="-"),colour="#94A3B8",size=3)+
  scale_shape_manual(values=c("Segment"=16,"Target-bin median"=0),name=NULL)+
  scale_colour_manual(values=c("Segment"=cook_rust,"Target-bin median"=cook_slate),name=NULL)+
  scale_x_continuous(limits=c(-.5,3.5),breaks=c(0,1,2,3),expand=c(0,0))+
  labs(title="A  CCNE1 DNA",x="Log2 ratio to\nautosomal baseline")
rn <- base()+geom_point(aes(x=rna_log2_tpm_plus1),colour=cook_rust,size=2,na.rm=TRUE)+
  geom_text(data=d[is.na(d$rna_tpm),],aes(x=2.4,label="-"),colour="#94A3B8",size=3)+
  scale_x_continuous(limits=c(2.2,7.6),breaks=log2(c(4,16,64,128)+1),labels=c("4","16","64","128"),expand=c(0,0))+
  labs(title="B  CCNE1 RNA",x="TPM\n(log scale)")
pr <- base()+geom_point(aes(x=protein_log2_normalised),colour=cook_slate,size=2,na.rm=TRUE)+
  geom_text(data=d[is.na(d$protein_log2_normalised),],aes(x=12.27,label="-"),colour="#94A3B8",size=3)+
  scale_x_continuous(limits=c(12.2,13.3),breaks=c(12.25,12.75,13.25),expand=c(0,0))+
  labs(title="C  Cyclin E1 protein",x="Log2 normalised abundance\n(isoDoping target)")
p <- (labels+cn+rn+pr)+plot_layout(widths=c(1.34,1.35,1.05,1.10),guides="collect") &
  theme(legend.position="bottom",legend.text=element_text(size=8))
dir.create(file.path(PROJ,"output/pdf"),recursive=TRUE,showWarnings=FALSE)
save_fig(p,file.path(PROJ,"output/pdf/ccne1_exploration.pdf"),w=7.2,h=9.3)
save_fig(p,file.path(PROJ,"output/pdf/ccne1_exploration.png"),w=7.2,h=9.3,dpi=220)
message("CCNE1 figure exported; 42 rows, 23 DNA, 31 RNA, 31 protein.")
