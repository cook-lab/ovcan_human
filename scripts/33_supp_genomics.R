# Supplementary Figures S5 and S6; figure redesign 2026-09-05.
# Full external legends: reports/figure_redesign_2026-09-05/legends_genomics.md.
source("scripts/00_setup.R"); source("scripts/00b_figure_theme.R")
suppressPackageStartupMessages({library(tidyverse); library(data.table); library(GenomicRanges)
  library(ComplexHeatmap); library(circlize); library(grid); library(patchwork)})
ht_opt_cook()
ASSET <- file.path(PROJ,"reports","assets")
# Font measurement uses metric-compatible Helvetica; exported glyphs use Arial.
if (!FIG_FONT %in% names(grDevices::pdfFonts()))
  do.call(grDevices::pdfFonts,setNames(list(grDevices::pdfFonts()$Helvetica),FIG_FONT))
if (!FIG_FONT %in% names(grDevices::postscriptFonts()))
  do.call(grDevices::postscriptFonts,setNames(list(grDevices::postscriptFonts()$Helvetica),FIG_FONT))
read_out <- function(x) read_csv(file.path(OUT,x),show_col_types=FALSE)

# S5: chromosome-wide overview; locus callouts and repeated family paragraphs removed.
seg <- as.data.table(read_out("wes_cnv_segments.csv"))
fam <- read_csv(file.path(META,"line_family_map.csv"),show_col_types=FALSE)
fga <- read_out("wes_cnv_fga.csv")
AUT <- paste0("chr",1:22)
HG38 <- c(chr1=248956422,chr2=242193529,chr3=198295559,chr4=190214555,chr5=181538259,
  chr6=170805979,chr7=159345973,chr8=145138636,chr9=138394717,chr10=133797422,
  chr11=135086622,chr12=133275309,chr13=114364328,chr14=107043718,chr15=101991189,
  chr16=90338345,chr17=83257441,chr18=80373285,chr19=58617616,chr20=64444167,
  chr21=46709983,chr22=50818468)
seg <- seg[chromosome %in% AUT]
meta <- seg[,.(subtype=subtype[1]),by=cell_line] %>%
  left_join(fam %>% select(cell_line,patient_id,family,is_multiline_family,has_wes_maf),by="cell_line") %>%
  left_join(fga %>% select(cell_line,fga=fga_auto_0.2),by="cell_line") %>%
  mutate(subtype=factor(subtype,levels=c("HGS","CC","EC","LGS","MC"))) %>%
  arrange(subtype,patient_id,cell_line)
stopifnot(nrow(meta)==23,all(meta$has_wes_maf))
bins <- tileGenome(HG38[AUT],tilewidth=1e7,cut.last.tile.in.chrom=TRUE)
chr <- factor(as.character(seqnames(bins)),levels=AUT,labels=1:22)
bin_sample <- function(s) {
 d <- seg[cell_line==s]
 gr <- GRanges(as.character(d$chromosome),IRanges(d$start+1L,d$end))
 h <- findOverlaps(bins,gr)
 w <- width(pintersect(bins[queryHits(h)],gr[subjectHits(h)]))
 q <- queryHits(h); v <- d$log2c_auto[subjectHits(h)]
 num <- tapply(v*w,q,sum); den <- tapply(w,q,sum)
 z <- rep(NA_real_,length(bins));z[as.integer(names(num))]<-num/den;z
}
m <- t(vapply(meta$cell_line,bin_sample,numeric(length(bins))));rownames(m)<-meta$cell_line
fp <- c(family_colours,Other=cook_grey)
fg <- ifelse(meta$is_multiline_family,meta$family,"Other")
check_palette_keys(fg,fp,"patient family")
tp <- droplevels(meta$subtype)
legtext <- gpar(fontsize=7.5,fontfamily=FIG_FONT,col=cook_ink)
legtitle <- gpar(fontsize=8,fontfamily=FIG_FONT,col=cook_ink)
la <- rowAnnotation(Patient=fg,col=list(Patient=fp),simple_anno_size=unit(2.7,"mm"),
 show_annotation_name=FALSE,
 annotation_legend_param=list(Patient=list(title="Patient family",nrow=1,at=names(fp),labels=names(fp),
   title_gp=legtitle,labels_gp=legtext,grid_height=unit(3,"mm"),grid_width=unit(3,"mm"))))
ra <- rowAnnotation(FGA=anno_barplot(meta$fga,baseline=0,bar_width=.75,
 gp=gpar(fill=cook_slate,col=NA),border=FALSE,ylim=c(0,.85),
 axis_param=list(side="bottom",at=c(0,.4,.8),labels_rot=0,gp=legtext),width=unit(15,"mm")),
 annotation_name_gp=legtitle,annotation_name_rot=0,annotation_name_side="top")
ht <- Heatmap(m,name="Log2 copy ratio",show_heatmap_legend=FALSE,col=cook_div_colfun(1.5),na_col=cook_grey,
 cluster_rows=FALSE,cluster_columns=FALSE,row_split=tp,row_title=c("HGS","CC","EC","LGS","MC"),row_title_rot=0,
 row_title_gp=gpar(fontsize=8,fontfamily=FIG_FONT,col=cook_ink),row_gap=unit(1.8,"mm"),
 column_split=chr,cluster_column_slices=FALSE,column_title_rot=90,
 column_title_gp=gpar(fontsize=7.5,fontfamily=FIG_FONT,col=cook_ink),column_gap=unit(.35,"mm"),
 border=TRUE,border_gp=gpar(col=cook_hair,lwd=.4),show_column_names=FALSE,
 row_names_side="left",row_names_gp=gpar(fontsize=7.6,fontfamily=FIG_FONT,col=cook_ink),
 left_annotation=la,right_annotation=ra,
 heatmap_legend_param=list(direction="horizontal",at=c(-1.5,0,1.5),labels=c("≤ -1.5","0","≥ 1.5"),
  title_gp=legtitle,labels_gp=legtext,legend_width=unit(35,"mm")))
missing_leg <- Legend(labels="No segment",legend_gp=gpar(fill=cook_grey,col=NA),
 labels_gp=legtext,grid_height=unit(3,"mm"),grid_width=unit(3,"mm"))
ratio_leg <- Legend(title="Log2 copy ratio",col_fun=cook_div_colfun(1.5),direction="horizontal",
 at=c(-1.5,0,1.5),labels=c("≤ -1.5","0","≥ 1.5"),legend_width=unit(35,"mm"),
 title_gp=legtitle,labels_gp=legtext)
draw_s5 <- function() draw(ht,heatmap_legend_side="bottom",annotation_legend_side="bottom",
 merge_legend=TRUE,heatmap_legend_list=list(missing_leg,ratio_leg),padding=unit(c(3,3,3,3),"mm"))
for (out in list(c(file.path(MSFIG,"figs5.pdf"),"pdf"),c(file.path(MSFIG,"figs5.png"),"png"),
                c(file.path(ASSET,"f_wes_cnv.png"),"png"))) {
 if(out[2]=="pdf") figure_pdf(out[1],width=W2,height=4.2)
 else ragg::agg_png(out[1],width=W2,height=4.2,units="in",res=400)
 draw_s5();dev.off()
}

# S6: probability matrix, margin plot and compact numbered exploratory clusters.
cons <- read_out("consensusov_calls.csv")
het <- read_out("hgs_heterogeneity.csv")
lv <- c("DIF","PRO","IMR","MES")
s <- cons %>% filter(applicable,!is.na(consensusov_call)) %>%
 left_join(het %>% select(cell_line,cluster,cluster_label),by="cell_line") %>%
 mutate(call=factor(consensusov_call,levels=lv)) %>%
 arrange(call,desc(margin_top_vs_second),cell_line) %>%
 mutate(cell_line=factor(cell_line,levels=rev(cell_line)))
stopifnot(nrow(s)==15,!anyNA(s$cluster))
clustermap <- s %>% distinct(cluster,cluster_label) %>% arrange(cluster)
write_csv(clustermap,file.path(PROJ,"reports/figure_redesign_2026-09-05/s6_cluster_key.csv"))
ylabels <- setNames(paste0(as.character(s$cell_line),ifelse(s$call_stable_across_input_sets,"","*")),s$cell_line)
pr <- s %>% select(cell_line,consensusov_call,all_of(paste0("prob_",lv))) %>%
 pivot_longer(all_of(paste0("prob_",lv)),names_to="class",values_to="probability") %>%
 mutate(class=factor(sub("prob_","",class),levels=lv),selected=as.character(class)==consensusov_call)
th <- function() theme_ovcan(8.5)+theme(axis.text=element_text(size=8),
 plot.margin=margin(5,5,5,5),plot.tag=element_text(size=11,face="bold",colour=cook_ink))
pA <- ggplot(pr,aes(class,cell_line,fill=probability)) +
 geom_tile(colour="white",linewidth=.5,width=.96,height=.92) +
 geom_tile(data=filter(pr,selected),fill=NA,colour=cook_ink,linewidth=.4,width=.96,height=.92) +
 geom_text(aes(label=sprintf("%.2f",probability),colour=probability>.58),size=2.7) +
 scale_colour_manual(values=c(`TRUE`="white",`FALSE`=cook_ink),guide="none") +
 scale_fill_gradient(low="white",high=cook_slate,limits=c(0,1),breaks=c(0,.5,1),
   name="Probability (uncalibrated)",guide=guide_colourbar(direction="horizontal",title.position="top",
     barwidth=unit(35,"mm"),barheight=unit(2.5,"mm"))) +
 scale_x_discrete(position="top",expand=c(0,0)) +scale_y_discrete(labels=ylabels,expand=expansion(add=.5)) +
 labs(x=NULL,y=NULL,tag="A")+th()+theme(axis.line=element_blank(),axis.ticks=element_blank(),
   legend.position="bottom",legend.justification="left",legend.title=element_text(size=7.5),legend.text=element_text(size=7.5))
pB <- ggplot(s,aes(margin_top_vs_second,cell_line)) +
 annotate("rect",xmin=0,xmax=.1,ymin=-Inf,ymax=Inf,fill=cook_rust,alpha=.09) +
 geom_vline(xintercept=median(s$margin_top_vs_second),colour=cook_slate,linetype=2,linewidth=.4) +
 geom_segment(aes(x=0,xend=margin_top_vs_second,yend=cell_line),colour=cook_hair,linewidth=.55) +
 geom_point(size=2.3,colour=cook_slate) +
 scale_x_continuous(limits=c(0,.75),breaks=c(0,.25,.5,.75),expand=c(0,0),position="top") +
 scale_y_discrete(expand=expansion(add=.5)) + labs(x="Top - second probability",y=NULL,tag="B") + th() +
 theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
       axis.title.x.top=element_text(size=8,margin=margin(b=5)))
pC <- ggplot(s,aes(x=1,y=cell_line,fill=cluster_label)) +
 geom_tile(width=.55,height=.92) +
 geom_text(aes(label=cluster,colour=cluster_label=="Low-signaling"),size=2.8) +
 scale_fill_stratum(guide="none") + scale_colour_manual(values=c(`TRUE`=cook_ink,`FALSE`="white"),guide="none") +
 scale_x_continuous(breaks=1,labels="RNA\ncluster",position="top",limits=c(.5,1.5),expand=c(0,0)) +
 scale_y_discrete(expand=expansion(add=.5)) + labs(x=NULL,y=NULL,tag="C") + th() +
 theme(axis.text.y=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),
       axis.text.x.top=element_text(size=8,lineheight=.9))
# Equal row ranges are essential: the margin endpoints must align to tile centres.
y_ranges <- lapply(list(pA,pB,pC), function(p) ggplot_build(p)$layout$panel_params[[1]]$y.range)
stopifnot(identical(y_ranges[[1]],y_ranges[[2]]),identical(y_ranges[[1]],y_ranges[[3]]))
figs6 <- (pA | pB | pC) + plot_layout(widths=c(3.6,2.65,.75)) & theme(text=element_text(family=FIG_FONT))
save_fig(figs6,file.path(MSFIG,"figs6.pdf"),w=W2,h=3.8)
save_fig(figs6,file.path(MSFIG,"figs6.png"),w=W2,h=3.8)
message("S5 and S6 complete; 23 CNV models and 15 HGSC RNA models.")
