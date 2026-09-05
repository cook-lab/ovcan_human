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
stopifnot(identical(unname(as.integer(n_assay)),c(31L,31L,23L)))
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
  geom_segment(aes(x=4.95,xend=4.95+models/3,y=y,yend=y),colour=cook_grey,linewidth=4,lineend="butt")+
  geom_point(aes(x=3.8,y=y,colour=subtype),size=2.1)+scale_colour_subtype(guide="none")+
  geom_text(aes(x=0,y=y,label=subtype),hjust=0,size=2.5,colour=cook_ink)+
  geom_text(aes(x=13.7,y=y,label=models),size=2.6,colour=cook_ink)+
  geom_text(aes(x=18.5,y=y,label=patients),size=2.6,colour=cook_ink)+
  annotate("text",x=c(13.7,18.5),y=8.05,label=c("Models","Patients"),size=2.5,colour=cook_ink,fontface="bold")+
  scale_x_continuous(limits=c(0,21),expand=c(0,0))+
  scale_y_continuous(limits=c(.45,8.6),expand=c(0,0))+
  theme_void(base_family=FIG_FONT)+theme(plot.margin=margin(3,3,5,3))
# Coverage uses one fixed-aspect canvas: square assay cells stay square when
# exported. Two adjacent blocks retain the HGSC/other split and patient brackets.
# A/B occupy the other figure column, so compact cells also reduce page area.
ordered <- fam %>% mutate(fs=ifelse(is_multiline_family,paste0("A_",family),paste0("Z_",cell_line))) %>%
  arrange(subtype,fs,desc(patient_representative),cell_line)
coverage <- ordered %>% mutate(block=ifelse(subtype=="HGS",1L,2L)) %>%
  group_by(block) %>%
  mutate(g=as.integer(factor(subtype,levels=unique(subtype))),
         y=25-row_number()-(g-1)*ifelse(block==2L,.7,0),
         x0=ifelse(block==1L,0,13.6),
         assay_x0=ifelse(block==1L,8.85,6.65)) %>% ungroup()
coverage_long <- coverage %>% select(cell_line,y,x0,assay_x0,has_rna,has_prot,has_wes_cnv) %>%
  pivot_longer(starts_with("has_"),names_to="assay",values_to="present") %>%
  mutate(x=x0+assay_x0+c(has_rna=0,has_prot=1.15,has_wes_cnv=2.3)[assay])
stopifnot(n_distinct(coverage_long$cell_line)==42,nrow(coverage_long)==126,
          sum(coverage_long$present)==sum(n_assay))
brk <- coverage %>% filter(is_multiline_family) %>% group_by(family) %>%
  summarise(lo=min(y)-.36,hi=max(y)+.36,mid=mean(y),n=n(),.groups="drop")
stopifnot(all(abs(brk$hi-brk$lo-.72-(brk$n-1))<1e-6))
headers <- coverage %>% distinct(x0,assay_x0) %>%
  tidyr::crossing(assay=c("RNA","Protein","Exome")) %>%
  mutate(x=x0+assay_x0+c(RNA=0,Protein=1.15,Exome=2.3)[assay])
pC <- ggplot()+
  geom_segment(data=coverage,aes(x=x0,xend=x0+assay_x0+2.72,y=y-.5,yend=y-.5),
               linewidth=.15,colour="#F1F5F9")+
  geom_text(data=coverage,aes(x=x0,y=y,label=cell_line),hjust=0,size=2.55,colour=cook_ink)+
  geom_tile(data=coverage,aes(x=x0+5.5,y=y,fill=subtype),width=.33,height=.78)+
  scale_fill_subtype(guide="none")+
  geom_segment(data=brk,aes(x=6.03,xend=6.03,y=lo,yend=hi),colour=cook_slate,linewidth=.3)+
  geom_segment(data=brk,aes(x=5.88,xend=6.03,y=lo,yend=lo),colour=cook_slate,linewidth=.3)+
  geom_segment(data=brk,aes(x=5.88,xend=6.03,y=hi,yend=hi),colour=cook_slate,linewidth=.3)+
  geom_text(data=brk,aes(x=6.18,y=mid,label=family),hjust=0,size=2.3,colour=cook_slate)+
  geom_tile(data=coverage_long,aes(x,y),fill=ifelse(coverage_long$present,cook_slate,"white"),
            width=.74,height=.74,colour=cook_hair,linewidth=.2)+
  annotate("text",x=c(0,13.6),y=25.4,label="Cell line",hjust=0,
           size=2.55,fontface="bold",colour=cook_ink)+
  geom_text(data=headers,aes(x=x,y=24.8,label=assay),angle=90,hjust=0,vjust=.5,
            size=2.55,fontface="bold",colour=cook_ink)+
  scale_x_continuous(limits=c(0,23.15),expand=c(0,0))+
  scale_y_continuous(limits=c(.35,27.7),expand=c(0,0))+
  coord_fixed(ratio=1)+
  theme_void(base_family=FIG_FONT)+theme(plot.margin=margin(3,3,3,7))
# A compact key uses the same square mark as the table; the shared histotype
# colours are identified directly in panel B.
key <- ggplot()+
  annotate("point",x=c(.2,3.4),y=1.18,shape=22,size=2.5,stroke=.25,
           fill=c(cook_slate,"white"),colour=cook_hair)+
  annotate("text",x=c(.52,3.72),y=1.18,hjust=0,size=2.5,colour=cook_ink,
           label=c("Profiled","Not profiled"))+
  annotate("text",x=.02,y=.25,hjust=0,size=2.5,colour=cook_ink,
           label=sprintf("All three assays: %d models",sum(fam$all_three)))+
  xlim(0,8)+ylim(0,1.7)+theme_void(base_family=FIG_FONT)+theme(plot.margin=margin(0,3,0,3))
left <- (pA+labs(tag="A")) / (pB+labs(tag="B")) / key +
  plot_layout(heights=c(1.35,1.7,.43))
fig1 <- (left | (pC+labs(tag="C"))) + plot_layout(widths=c(3.2,3.8)) &
  theme(plot.tag=element_text(size=11,face="bold",family=FIG_FONT,colour=cook_ink),plot.tag.position="topleft")
save_fig(fig1,file.path(MSFIG,"fig1.pdf"),w=W2,h=4)
save_fig(fig1,file.path(MSFIG,"fig1.png"),w=W2,h=4)
message("Figure 1: 42 models, 34 patients; RNA/protein/exome = ",paste(n_assay,collapse="/"),"; all three = 13")
