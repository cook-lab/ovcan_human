# Figure 1: workflow, histotype counts and model-level assay coverage.
source("scripts/00_setup.R"); source("scripts/00b_figure_theme.R")
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
fam <- read_csv(file.path(META,"line_family_map.csv"),show_col_types=FALSE) %>%
  mutate(subtype=factor(subtype,levels=c("HGS","CC","MC","EC","MMMT","SCCOHT","LGS")),
         all_three=has_rna & has_prot & has_wes_cnv)
stopifnot(nrow(fam)==42,n_distinct(fam$patient_id)==34,sum(fam$all_three)==13,
          identical(fam$has_wes_cnv,fam$has_wes_maf))
check_palette_keys(fam$subtype,subtype_colours,"histotype")
n_assay <- colSums(as.data.frame(fam[c("has_rna","has_prot","has_wes_cnv")]))
boxes <- tribble(~x,~y,~w,~h,~label,~fill,~border,
  10,50,9,16,"Ovarian\ntumours","#F1F5F9",cook_slate,
  33,50,10,23,sprintf("%d cell lines\n%d patients",nrow(fam),n_distinct(fam$patient_id)),"#FEF3EE",cook_rust,
  63,82,13,12,sprintf("RNA-seq\nn = %d",n_assay[1]),"#FEF3EE",cook_rust,
  63,50,13,12,sprintf("Proteomics\nn = %d",n_assay[2]),"#F0FDFA",cook_teal,
  63,18,13,12,sprintf("Exome\nn = %d",n_assay[3]),"#F1F5F9",cook_slate,
  91,50,8,23,"Data\nresource","#F1F5F9",cook_slate)
edges <- tribble(~x,~y,~xend,~yend,19,50,23,50,43,56,50,82,43,50,50,50,43,44,50,18,
                 76,82,83,56,76,50,83,50,76,18,83,44)
pA <- ggplot()+
  geom_segment(data=edges,aes(x,y,xend=xend,yend=yend),colour=cook_slate,linewidth=.35,
               arrow=arrow(length=unit(3,"pt"),type="closed"))+
  geom_rect(data=boxes,aes(xmin=x-w,xmax=x+w,ymin=y-h,ymax=y+h),
            fill=boxes$fill,colour=boxes$border,linewidth=.45)+
  geom_text(data=boxes,aes(x,y,label=label),size=2.5,lineheight=1,colour=cook_ink,fontface="bold")+
  scale_x_continuous(limits=c(0,100),expand=c(0,0))+
  scale_y_continuous(limits=c(1,99),expand=c(0,0))+
  theme_void(base_family=FIG_FONT)+theme(plot.margin=margin(3,10,5,3))
# Count columns distinguish cell lines from independent patients.
counts <- fam %>% group_by(subtype) %>%
  summarise(models=n(),patients=n_distinct(patient_id),.groups="drop") %>% mutate(y=8-as.integer(subtype))
stopifnot(sum(counts$models)==42,sum(counts$patients)==34)
pB <- ggplot(counts)+
  geom_segment(aes(x=4.5,xend=4.5+models/3,y=y,yend=y),colour=cook_grey,linewidth=4,lineend="butt")+
  geom_point(aes(x=3.3,y=y,colour=subtype),size=2.1)+scale_colour_subtype(guide="none")+
  geom_text(aes(x=0,y=y,label=subtype),hjust=0,size=2.5,colour=cook_ink)+
  geom_text(aes(x=13.7,y=y,label=models),size=2.6,colour=cook_ink)+
  geom_text(aes(x=18.5,y=y,label=patients),size=2.6,colour=cook_ink)+
  annotate("text",x=c(13.7,18.5),y=8.05,label=c("Models","Patients"),size=2.5,colour=cook_ink,fontface="bold")+
  scale_x_continuous(limits=c(0,21),expand=c(0,0))+
  scale_y_continuous(limits=c(.45,8.6),expand=c(0,0))+
  theme_void(base_family=FIG_FONT)+theme(plot.margin=margin(3,3,5,3))
# Split coverage across two blocks to preserve readable model labels.
ordered <- fam %>% mutate(fs=ifelse(is_multiline_family,paste0("A_",family),paste0("Z_",cell_line))) %>%
  arrange(subtype,fs,desc(patient_representative),cell_line)
make_coverage <- function(d,nonhgs=FALSE) {
  d <- d %>% mutate(g=as.integer(factor(subtype,levels=unique(subtype))),
                    y=25-row_number()-(g-1)*ifelse(nonhgs,.7,0))
  long <- d %>% select(cell_line,y,has_rna,has_prot,has_wes_cnv) %>%
    pivot_longer(starts_with("has_"),names_to="assay",values_to="present") %>%
    mutate(x=c(has_rna=9,has_prot=11.5,has_wes_cnv=14)[assay])
  brk <- d %>% filter(is_multiline_family) %>% group_by(family) %>%
    summarise(lo=min(y)-.36,hi=max(y)+.36,mid=mean(y),n=n(),.groups="drop")
  stopifnot(all(abs(brk$hi-brk$lo-.72-(brk$n-1))<1e-6))
  ggplot()+
    geom_segment(data=d,aes(x=0,xend=14.6,y=y-.5,yend=y-.5),linewidth=.15,colour="#F1F5F9")+
    geom_text(data=d,aes(x=0,y=y,label=cell_line),hjust=0,size=2.55,colour=cook_ink)+
    geom_tile(data=d,aes(x=5.25,y=y,fill=subtype),width=.33,height=.78)+scale_fill_subtype(guide="none")+
    geom_segment(data=brk,aes(x=5.85,xend=5.85,y=lo,yend=hi),colour=cook_slate,linewidth=.3)+
    geom_segment(data=brk,aes(x=5.7,xend=5.85,y=lo,yend=lo),colour=cook_slate,linewidth=.3)+
    geom_segment(data=brk,aes(x=5.7,xend=5.85,y=hi,yend=hi),colour=cook_slate,linewidth=.3)+
    geom_text(data=brk,aes(x=6.05,y=mid,label=family),hjust=0,size=2.3,colour=cook_slate)+
    geom_tile(data=long,aes(x,y),fill=ifelse(long$present,cook_slate,"white"),
              width=1.85,height=.74,colour=cook_hair,linewidth=.2)+
    annotate("text",x=c(0,9,11.5,14),y=25.5,label=c("Cell line","RNA","Protein","Exome"),
             hjust=c(0,.5,.5,.5),size=2.55,fontface="bold",colour=cook_ink)+
    scale_x_continuous(limits=c(0,15.1),expand=c(0,0))+
    scale_y_continuous(limits=c(.35,26.1),expand=c(0,0))+
    theme_void(base_family=FIG_FONT)+theme(plot.margin=margin(2,6,1,3))
}
pC1 <- make_coverage(filter(ordered,subtype=="HGS"))
pC2 <- make_coverage(filter(ordered,subtype!="HGS"),TRUE)
key <- ggplot()+
  annotate("rect",xmin=c(0,2.4),xmax=c(.28,2.68),ymin=.35,ymax=.7,
           fill=c(cook_slate,"white"),colour=cook_hair,linewidth=.25)+
  annotate("text",x=c(.42,2.82,5.6),y=.525,hjust=0,size=2.5,colour=cook_ink,
           label=c("Profiled","Not profiled",sprintf("All three assays: %d models",sum(fam$all_three))))+
  xlim(0,10.5)+ylim(0,1)+theme_void(base_family=FIG_FONT)+theme(plot.margin=margin(0,3,0,3))
fig1 <- ((pA+labs(tag="A")) | (pB+labs(tag="B"))) /
  ((pC1+labs(tag="C")) | pC2) / key + plot_layout(heights=c(1.38,3.85,.27)) &
  theme(plot.tag=element_text(size=11,face="bold",family=FIG_FONT,colour=cook_ink),plot.tag.position="topleft")
save_fig(fig1,file.path(MSFIG,"fig1.pdf"),w=W2,h=5.55)
save_fig(fig1,file.path(MSFIG,"fig1.png"),w=W2,h=5.55)
message("Figure 1: 42 models, 34 patients; RNA/protein/exome = ",paste(n_assay,collapse="/"),"; all three = 13")
