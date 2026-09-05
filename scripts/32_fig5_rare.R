# Figure 5: exploratory mutation spectra and model-specific molecular evidence.
# Figure redesign 2026-09-05. Narrative/methods are in legends_genomics.md.
source("scripts/00_setup.R"); source("scripts/00b_figure_theme.R")
suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
ASSET <- file.path(PROJ,"reports","assets")
th <- function() theme_ovcan(8.5) + theme(axis.text=element_text(size=7.5),
  legend.text=element_text(size=7.5),legend.title=element_text(size=7.5),
  plot.margin=margin(5,5,5,5),plot.tag=element_text(size=11,face="bold",colour=cook_ink))
read_out <- function(x) read_csv(file.path(OUT,x),show_col_types=FALSE)
loads <- read_out("wes_mutation_load.csv")
msi <- read_out("wes_msi_mmr.csv")
sbs <- read_out("wes_sbs_context.csv")
cosines <- read_out("wes_sbs_cosine.csv")
sensitivity <- read_out("wes_signature_target_sensitivity.csv")
markers <- read_out("auth_mucinous_marker_ranks.csv")
swi <- read_out("auth_swisnf_long.csv")
swip <- read_out("auth_swisnf_panel.csv")
HYP <- loads$cell_line[loads$is_hypermutator]
stopifnot(length(HYP)==1,nrow(loads)==23)

# A: a vertical model list avoids 23 rotated labels and the outlier-driven empty field.
lo <- loads %>% arrange(n_coding) %>% mutate(cell_line=factor(cell_line,levels=cell_line))
ifr <- c(floor(min(lo$indel_frac)*50)/50,ceiling(max(lo$indel_frac)*50)/50)
pA <- ggplot(lo,aes(n_coding,cell_line)) +
  geom_segment(aes(x=100,xend=n_coding,yend=cell_line),colour=cook_hair,linewidth=.45) +
  geom_vline(xintercept=median(lo$n_coding),colour=cook_slate,linetype=2,linewidth=.4) +
  geom_point(aes(fill=indel_frac),shape=21,size=2.35,stroke=.25,colour=cook_slate) +
  scale_x_log10(breaks=c(100,500,1500),labels=c("100","500","1,500"),
               limits=c(95,1650),expand=c(0,0)) +
  scale_y_discrete(expand=expansion(add=.55)) +
  scale_fill_gradientn(colours=cook_diverging[4:7],name="Indel fraction",limits=ifr,breaks=ifr,
    labels=sprintf("%.2f",ifr),guide=guide_colourbar(direction="horizontal",title.position="top",
      barwidth=unit(27,"mm"),barheight=unit(2.5,"mm"))) +
  labs(x="Coding candidates",y=NULL,tag="A") + th() +
  theme(legend.position="bottom",legend.justification="left",legend.margin=margin(0,0,0,0),
        axis.ticks.y=element_blank(),axis.line.y=element_blank())

# B: all 96 channels, grouped by substitution class; one subdued fill is sufficient.
cls <- c("C>A","C>G","C>T","T>A","T>C","T>G")
b <- sbs %>% transmute(context,n=.data[[HYP]]) %>%
  mutate(class=factor(sub("^.\\[(.*)\\].$","\\1",context),levels=cls),
         five=substr(context,1,1),three=substr(context,nchar(context),nchar(context)),frac=n/sum(n)) %>%
  arrange(class,five,three) %>% mutate(idx=row_number())
NSNV <- sum(b$n)
stopifnot(NSNV==msi$n_snv_used[msi$cell_line==HYP])
pB <- ggplot(b,aes(idx,frac)) + geom_col(fill=cook_slate,width=.88) +
  facet_grid(~class,scales="free_x",space="free_x") +
  scale_y_continuous(breaks=c(0,.05,.1),limits=c(0,.105),expand=c(0,0)) +
  labs(x="Trinucleotide context",y="SNV fraction",title=sprintf("%s (%s SNVs)",HYP,scales::comma(NSNV)),tag="B") +
  th() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
    strip.text=element_text(size=7.5),panel.spacing=unit(2,"pt"),
    plot.title=element_text(size=8.5),plot.margin=margin(5,5,3,5))

# C: strongest spectrum similarity in each prespecified signature group.
grp <- c("MMR-d / MSI","clock-like (age)","POLE proofreading")
c <- cosines %>% filter(cell_line==HYP,group %in% grp) %>% group_by(group) %>%
  slice_max(cosine,n=1,with_ties=FALSE) %>% ungroup() %>%
  mutate(group=factor(group,levels=rev(grp)),
    label=paste0(signature," (",recode(as.character(group),`MMR-d / MSI`="MMR",`clock-like (age)`="clock",`POLE proofreading`="POLE"),")"))
pC <- ggplot(c,aes(cosine,group)) +
  geom_segment(aes(x=0,xend=cosine,yend=group),colour=cook_hair,linewidth=.6) +
  geom_point(colour=cook_slate,size=2.3) +
  geom_text(aes(label=sprintf("%.3f",cosine)),nudge_x=.035,hjust=0,size=2.6,colour=cook_ink) +
  scale_y_discrete(labels=setNames(c$label,c$group),expand=expansion(add=.55)) +
  scale_x_continuous(limits=c(0,1.06),breaks=c(0,.5,1),expand=c(0,0)) +
  labs(x="Cosine similarity",y=NULL,tag="C") + th() +
  theme(axis.ticks.y=element_blank(),axis.line.y=element_blank(),plot.margin=margin(2,5,3,5))

# D: same MMR-associated estimand across three refit settings.
modes <- c("target-restricted variants; genome reference",
  "target-restricted variants; target-opportunity-adjusted reference",
  "target-restricted variants; target-adjusted restricted reference")
d <- sensitivity %>% filter(cell_line==HYP) %>% mutate(mode=factor(mode,levels=rev(modes)))
stopifnot(nrow(d)==3)
pD <- ggplot(d,aes(mmr_relative_exposure,mode)) +
  geom_segment(aes(x=0,xend=mmr_relative_exposure,yend=mode),colour=cook_hair,linewidth=.6) +
  geom_point(colour=cook_rust,size=2.3) +
  geom_text(aes(label=sprintf("%.3f",mmr_relative_exposure)),nudge_x=.035,hjust=0,size=2.6,colour=cook_ink) +
  scale_y_discrete(labels=setNames(c("Genome / full","Target / full","Target / restricted"),modes),
                   expand=expansion(add=.55)) +
  scale_x_continuous(limits=c(0,1.06),breaks=c(0,.5,1),expand=c(0,0)) +
  labs(x="MMR-associated fitted fraction",y=NULL,tag="D") + th() +
  theme(axis.ticks.y=element_blank(),axis.line.y=element_blank(),plot.margin=margin(2,5,5,5))

# Shared display range for the two expression panels, with explicit endpoint labels.
zscale <- function() scale_fill_cook_div(midpoint=0,name="Within-assay z score",limits=c(-3,3),
  oob=scales::squish,na.value=cook_grey,breaks=c(-3,0,3),labels=c("≤ -3","0","≥ 3"),
  guide=guide_colourbar(direction="horizontal",title.position="top",barwidth=unit(27,"mm"),barheight=unit(2.5,"mm")))
# E: marker groups are explained in the external legend, removing redundant side labels.
genes <- c("KRT7","PAX8","WT1","SATB2","CDX2","KRT20","MUC2","MUC5AC","TFF1","TFF3")
lines <- sort(unique(markers$cell_line[markers$is_mucinous_line]))
e <- markers %>% filter(is_mucinous_line,symbol %in% genes) %>%
  select(cell_line,symbol,assay,z) %>% complete(cell_line=lines,symbol=genes,assay=c("RNA","protein")) %>%
  mutate(symbol=factor(symbol,levels=rev(genes)),cell_line=factor(cell_line,levels=lines),
         assay=factor(assay,levels=c("RNA","protein"),labels=c("RNA","Protein")))
pE <- ggplot(e,aes(cell_line,symbol,fill=z)) +
  geom_tile(colour="white",linewidth=.55) +
  geom_point(data=filter(e,is.na(z)),shape=4,size=1.7,stroke=.35,colour=cook_slate) +
  facet_grid(~assay) + zscale() + scale_x_discrete(position="top",expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) + labs(x=NULL,y=NULL,tag="E") + th() +
  theme(axis.text.x.top=element_text(angle=45,hjust=0,size=7.5),axis.text.y=element_text(size=7.8,face="italic"),
    strip.text=element_text(size=8.3),panel.spacing=unit(4,"pt"),axis.line=element_blank(),
    axis.ticks=element_blank(),legend.position="bottom")

# F: assay ranks stay visible; compact WES labels preserve priority and exclusions.
sg <- c("SMARCA4","SMARCA2","ARID1A")
show <- swip %>% filter(swisnf_deficient | !is.na(swisnf_tier3_only_calls)) %>%
  arrange(desc(swisnf_deficient),subtype,cell_line)
f <- swi %>% filter(cell_line %in% show$cell_line,gene %in% sg) %>%
  mutate(cell_line=factor(cell_line,levels=rev(show$cell_line)),gene=factor(gene,levels=sg),
    assay=factor(assay,levels=c("RNA","protein")),x=2*(as.integer(gene)-1)+as.integer(assay),y=as.integer(cell_line))
wt <- swi %>% filter(cell_line %in% show$cell_line,gene %in% sg,!is.na(wes_class),assay=="RNA") %>%
  transmute(cell_line,gene,lab=sprintf("%s T%d%s",gene,wes_tier_best,
    ifelse(wes_class=="truncating" & !wes_trunc_admitted,"*",""))) %>%
  group_by(cell_line) %>% summarise(lab=paste(lab,collapse="\n"),.groups="drop") %>%
  mutate(y=match(cell_line,levels(f$cell_line)))
NY <- nlevels(f$cell_line)
heads <- tibble(gene=sg,x=c(1.5,3.5,5.5))
pF <- ggplot(f,aes(x,y,fill=z)) +
  geom_tile(width=.94,height=.9,colour="white",linewidth=.45) +
  geom_text(aes(label=rank_low_is_lowest,colour=abs(z)>1.8),size=2.6) +
  scale_colour_manual(values=c(`TRUE`="white",`FALSE`=cook_ink),guide="none") + zscale() +
  geom_text(data=heads,aes(x,y=NY+1.5,label=gene),inherit.aes=FALSE,
            size=2.6,fontface="italic",colour=cook_ink) +
  geom_segment(data=heads,aes(x=x-.94,xend=x+.94,y=NY+1.18,yend=NY+1.18),inherit.aes=FALSE,
               colour=cook_hair,linewidth=.35) +
  annotate("text",x=1:6,y=NY+.8,label=rep(c("RNA","Prot"),3),size=2.45,colour=cook_ink) +
  geom_text(data=wt,aes(x=6.9,y,label=lab),inherit.aes=FALSE,size=2.5,hjust=0,lineheight=1.1,colour=cook_ink) +
  annotate("text",x=6.9,y=NY+1.5,label="WES tier",size=2.6,hjust=0,colour=cook_ink) +
  scale_y_continuous(breaks=1:NY,labels=levels(f$cell_line),expand=c(0,0)) +
  scale_x_continuous(breaks=NULL,expand=c(0,0)) +
  coord_cartesian(xlim=c(.5,9.4),ylim=c(.45,NY+1.85),clip="off") +
  labs(x=NULL,y=NULL,tag="F") + th() +
  theme(axis.text.y=element_text(size=7.8),axis.line=element_blank(),axis.ticks=element_blank(),
        legend.position="bottom")

right_top <- (pB / pC / pD) + plot_layout(heights=c(1.25,1,1))
top <- (pA | right_top) + plot_layout(widths=c(2.5,4.7))
bottom <- (pE | pF) + plot_layout(widths=c(2.9,4.3),guides="collect") &
  theme(legend.position="bottom",legend.justification="left")
fig5 <- (top / bottom) + plot_layout(heights=c(3.85,2.7)) & theme(text=element_text(family=FIG_FONT))
save_fig(fig5,file.path(MSFIG,"fig5.pdf"),w=W2,h=6.9)
save_fig(fig5,file.path(MSFIG,"fig5.png"),w=W2,h=6.9)
save_fig(top,file.path(ASSET,"f_wes_hypermutation.png"),w=W2,h=3.85)
save_fig(pE,file.path(ASSET,"f_auth_mucinous.png"),w=3.6,h=3.0)
message("Figure 5 complete: six panels; data unchanged; interpretive notes external.")
