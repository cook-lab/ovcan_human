# Separate exploratory figures from script42. No canonical pipeline outputs.
PROJ <- normalizePath(Sys.getenv("OVCAN_PROJ", unset=getwd()), mustWork=TRUE)
source(file.path(PROJ,"scripts/00b_figure_theme.R"))
suppressPackageStartupMessages(library(patchwork))
D <- file.path(PROJ,"reports/molecular_extension_2026-09-06/expression")
OUT <- file.path(PROJ,"output/pdf")
read <- function(n) read.csv(file.path(D,n),check.names=FALSE,na.strings=c("NA",""))
panel <- read("panel.csv"); raw <- read("model_expression.csv")
con <- read("concordance.csv"); rank <- read("joint_ranking.csv")
summary <- read("target_summary.csv")
genes <- panel$gene; ng <- length(genes)
stopifnot(ng==19,nrow(raw)==42*ng,all(table(raw$gene)==42))
pos <- setNames(rev(seq_along(genes)),genes)
flag <- setNames(!is.na(summary$protein_isodoping)&tolower(summary$protein_isodoping)=="true",summary$gene)
gene_labels <- setNames(paste0(genes,ifelse(flag[genes],"*","")),genes)
groups <- panel$panel_group
group_y <- pos[genes[c(FALSE,groups[-1]!=groups[-length(groups)])]]+.5
stripes <- data.frame(y=seq_len(ng)[seq_len(ng)%%2==0])
common <- function() {
  theme_ovcan(base_size=8)+theme(plot.title=element_text(size=9,face="bold",margin=margin(b=9)),
    axis.title=element_text(size=8),axis.text=element_text(size=7.5),
    axis.line.y=element_blank(),axis.ticks.y=element_blank(),
    legend.position="bottom",legend.text=element_text(size=7),
    plot.margin=margin(5,5,5,5))
}
row_bg <- function() geom_rect(data=stripes,aes(ymin=y-.5,ymax=y+.5,xmin=-Inf,xmax=Inf),
                               inherit.aes=FALSE,fill="#F4F6F8",colour=NA)
a <- con[con$cohort %in% c("selected_patient_representatives","HGS_selected_representatives"),]
a$set <- factor(ifelse(a$cohort=="selected_patient_representatives","27 patient models","12 HGSC patient models"),
                levels=c("27 patient models","12 HGSC patient models"))
a$y <- pos[a$gene]+ifelse(a$cohort=="selected_patient_representatives",.14,-.14)
no_prot <- data.frame(gene=summary$gene[summary$n_protein_models==0])
no_prot$y <- pos[no_prot$gene]
pa <- ggplot(a,aes(spearman,y))+row_bg()+
  geom_hline(yintercept=group_y,colour=cook_hair,linewidth=.3)+
  geom_vline(xintercept=0,colour=cook_hair,linewidth=.4)+
  geom_segment(data=a[a$cohort=="selected_patient_representatives",],
               aes(x=bootstrap_lo,xend=bootstrap_hi,yend=y),colour="#D7B5A7",linewidth=.7,na.rm=TRUE)+
  geom_point(aes(colour=set,shape=set),size=1.8,na.rm=TRUE)+
  geom_text(data=no_prot,aes(x=.28,y=y,label="RNA only"),inherit.aes=FALSE,
            size=2.5,colour=cook_ink_muted)+
  scale_colour_manual(values=c("27 patient models"=cook_rust,"12 HGSC patient models"=cook_slate),name=NULL)+
  scale_shape_manual(values=c("27 patient models"=16,"12 HGSC patient models"=18),name=NULL)+
  scale_y_continuous(limits=c(.5,ng+.5),breaks=pos[genes],labels=gene_labels[genes],expand=c(0,0))+
  scale_x_continuous(limits=c(-1.02,1.02),breaks=c(-1,-.5,0,.5,1),expand=c(0,0))+
  labs(title="A  RNA-protein agreement",x="Spearman correlation",y=NULL)+common()+
  theme(legend.direction="vertical",legend.box="vertical",legend.key.height=unit(9,"pt"))

primary <- rank[rank$cohort=="selected_patient_representatives",]
primary <- primary[order(match(primary$gene,genes),primary$joint_mean_rank,primary$cell_line),]
lead <- primary[!duplicated(primary$gene),]
pm <- rank[rank$cohort=="patient_means_of_matched_models",c("gene","patient_id","joint_mean_rank")]
names(pm)[3] <- "patient_mean_rank"
lead <- merge(lead,pm,by=c("gene","patient_id"),all.x=TRUE,sort=FALSE)
lead$y <- pos[lead$gene]
max_rank <- max(10,ceiling(max(lead$patient_mean_rank,lead$rank_worst_across_methods,na.rm=TRUE)))
points <- rbind(data.frame(y=lead$y,rank=lead$joint_mean_rank,set="Selected model"),
                data.frame(y=lead$y,rank=lead$patient_mean_rank,set="Patient mean"))
points$set <- factor(points$set,levels=c("Selected model","Patient mean"))
pb <- ggplot()+row_bg()+geom_hline(yintercept=group_y,colour=cook_hair,linewidth=.3)+
  geom_vline(xintercept=1,colour=cook_hair,linewidth=.4)+
  geom_segment(data=lead,aes(x=rank_best_across_methods,xend=rank_worst_across_methods,y=y,yend=y),
               colour=cook_slate,linewidth=1)+
  geom_point(data=points,aes(rank,y,shape=set,colour=set),size=2)+
  geom_text(data=lead,aes(x=-9,y=y,label=cell_line),hjust=0,size=2.7,colour=cook_ink)+
  geom_text(data=no_prot,aes(x=-9,y=y,label="No paired protein"),hjust=0,size=2.5,colour=cook_ink_muted)+
  scale_colour_manual(values=c("Selected model"=cook_rust,"Patient mean"=cook_slate),name=NULL)+
  scale_shape_manual(values=c("Selected model"=16,"Patient mean"=5),name=NULL)+
  scale_x_continuous(limits=c(-9.5,max_rank+.5),breaks=c(1,3,5,7,9),expand=c(0,0))+
  scale_y_continuous(limits=c(.5,ng+.5),breaks=NULL,expand=c(0,0))+
  labs(title="B  Leading model rank sensitivity",x="Joint rank among 27 patients",y=NULL)+common()+
  theme(axis.line.x=element_blank(),legend.direction="vertical",legend.key.height=unit(9,"pt"))
fig1 <- (pa+pb)+plot_layout(widths=c(1.02,1))
save_fig(fig1,file.path(OUT,"expression_concordance_rank_sensitivity.pdf"),7.2,6.8)
save_fig(fig1,file.path(OUT,"expression_concordance_rank_sensitivity.png"),7.2,6.8,dpi=220)

# Full model-level view: each assay is ranked among its own measured models.
# This display uses HIGH abundance at the high end for every target, including
# MMR/SWI-SNF, unlike the explicitly lower-directed candidate score above.
m <- unique(raw[,c("cell_line","histotype_code")])
hist <- c("HGS","CC","EC","LGS","MC","MMMT","SCCOHT")
m <- m[order(match(m$histotype_code,hist),m$cell_line),]
m$y <- rev(seq_len(nrow(m)))
raw$y <- m$y[match(raw$cell_line,m$cell_line)]
raw$x <- match(raw$gene,genes)
bound_y <- m$y[c(FALSE,m$histotype_code[-1]!=m$histotype_code[-nrow(m)])]+.5
bound_x <- which(c(FALSE,groups[-1]!=groups[-length(groups)]))-.5
heat <- function(assay,title,labels=FALSE) {
  q <- raw
  q$value <- 100*q[[if(assay=="RNA") "rna_abundance_percentile" else "protein_abundance_percentile"]]
  ggplot(q,aes(x,y))+geom_tile(aes(fill=value),width=.96,height=.94,colour=NA)+
    geom_point(data=q[is.na(q$value),],shape=4,size=.85,stroke=.22,colour="#94A3B8")+
    geom_hline(yintercept=bound_y,colour=cook_hair,linewidth=.3)+
    geom_vline(xintercept=bound_x,colour=cook_hair,linewidth=.4)+
    scale_fill_gradientn(colours=c("#EEF2F6","#CBD5E1","#F7BEA2",cook_rust),
                         limits=c(0,100),breaks=c(0,25,50,75,100),na.value="white",
                         name="Within-target abundance percentile")+
    scale_x_continuous(limits=c(.5,ng+.5),breaks=seq_along(genes),labels=gene_labels[genes],expand=c(0,0))+
    scale_y_continuous(limits=c(.5,nrow(m)+.5),breaks=if(labels)m$y else NULL,
      labels=if(labels)paste0(m$cell_line,"  ",m$histotype_code) else NULL,expand=c(0,0))+
    labs(title=title,x=NULL,y=NULL)+common()+theme(axis.line=element_blank(),axis.ticks=element_blank(),
      axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=7.2),axis.text.y=element_text(size=7.1),
      legend.title=element_text(size=8),legend.position="bottom")+
    guides(fill=guide_colourbar(barwidth=unit(2.6,"in"),barheight=unit(.09,"in"),title.position="top"))
}
fig2 <- (heat("RNA","A  RNA (31 models)",TRUE)+heat("Protein","B  Protein (31 models)",FALSE))+
  plot_layout(widths=c(1,1),guides="collect") & theme(legend.position="bottom")
save_fig(fig2,file.path(OUT,"expression_target_model_atlas.pdf"),7.2,9.6)
save_fig(fig2,file.path(OUT,"expression_target_model_atlas.png"),7.2,9.6,dpi=220)
message("Expression extension: 2 PDFs/PNGs, 19 targets, 42 model rows; legend text is external.")
