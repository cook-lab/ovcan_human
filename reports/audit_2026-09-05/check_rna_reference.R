suppressPackageStartupMessages({library(data.table); library(SummarizedExperiment); library(matrixStats)})
base_dir <- 'reports/audit_2026-09-05'
read_matrix <- function(path, key = 1L) {
  d <- fread(path); m <- as.matrix(d[, -key, with = FALSE]); rownames(m) <- d[[key]]; m
}
old <- read_matrix(file.path(base_dir, 'rna_before/rna_tpm.csv'))
new <- read_matrix('output/rna_tpm.csv')
om <- fread('data/reference/tx2gene_ensembl_rel105.csv')
nm <- fread('output/tx2gene_matched.csv')
to_symbol <- function(m, map) {
  meta <- unique(map[, .(ensembl_gene_id, external_gene_name)])
  meta <- meta[!is.na(external_gene_name) & external_gene_name != '']
  stopifnot(!anyDuplicated(meta$ensembl_gene_id))
  idx <- match(rownames(m), meta$ensembl_gene_id); keep <- !is.na(idx)
  rowsum(m[keep, , drop = FALSE], meta$external_gene_name[idx[keep]])
}
os <- to_symbol(old, om); ns <- to_symbol(new, nm)
common <- intersect(rownames(os), rownames(ns)); lines <- intersect(colnames(os), colnames(ns))
profile_rho <- vapply(lines, function(l) cor(os[common,l], ns[common,l], method='spearman'), numeric(1))
oldpc <- readRDS(file.path(base_dir,'rna_before/rna_pca.rds'))
newpc <- readRDS('output/rna_pca.rds')
shared <- intersect(rownames(old), rownames(new))
summary <- data.table(metric = c('TPM_gene_rows_before','TPM_gene_rows_after','shared_gene_ids','new_gene_ids','retired_gene_ids',
                                 'symbol_rows_before','symbol_rows_after','shared_symbols','new_symbols','removed_symbols',
                                 'median_sample_profile_Spearman','min_sample_profile_Spearman','abs_PC1_score_correlation'),
                      value=c(nrow(old),nrow(new),length(shared),sum(!rownames(new)%in%rownames(old)),sum(!rownames(old)%in%rownames(new)),
                              nrow(os),nrow(ns),length(common),sum(!rownames(ns)%in%rownames(os)),sum(!rownames(os)%in%rownames(ns)),
                              median(profile_rho),min(profile_rho),abs(cor(oldpc$x[rownames(newpc$x),1],newpc$x[,1]))))
fwrite(summary,file.path(base_dir,'rna_reference_before_after.csv'))
fwrite(data.table(cell_line=lines, symbol_profile_spearman=profile_rho),file.path(base_dir,'rna_reference_profile_sensitivity.csv'))
# RNA abundance that was discarded because the old map omitted the transcript.
ss <- fread('output/rna_sample_annotation.csv')
unmapped <- nm[!ensembl_transcript_id %in% om$ensembl_transcript_id]
restored <- rbindlist(lapply(seq_len(nrow(ss)),function(i) {
  ab <- fread(ss$file[i]); x <- ab[match(unmapped$transcript_id_versioned,ab$target_id)]
  data.table(cell_line=ss$cell_line[i],gene_id=unmapped$ensembl_gene_id,symbol=unmapped$external_gene_name,
             primary_assembly=unmapped$primary_assembly,tpm_restored=x$tpm,counts_restored=x$est_counts)
}))[, .(tpm_restored=sum(tpm_restored),counts_restored=sum(counts_restored)),by=.(cell_line,gene_id,symbol,primary_assembly)]
fwrite(restored,file.path(base_dir,'rna_restored_gene_abundance.csv'))
top <- restored[,.(median_TPM_restored=median(tpm_restored),max_TPM_restored=max(tpm_restored)),by=.(gene_id,symbol,primary_assembly)][order(-median_TPM_restored)]
fwrite(top,file.path(base_dir,'rna_restored_gene_summary.csv'))
fo <- rbindlist(lapply(c('FOXL2','SMARCA4','SMARCA2','PAX8','FOLR1'),function(g) data.table(
  symbol=g,cell_line=lines,TPM_before=if(g%in%rownames(os))os[g,lines] else NA_real_,TPM_after=ns[g,lines])))
fwrite(fo,file.path(base_dir,'rna_marker_reference_delta.csv'))
# Sensitivity to primary-assembly-only features in otherwise identical VST/PCA.
v <- assay(readRDS('output/rna_vst.rds'))
primary_genes <- unique(nm[primary_assembly==TRUE,ensembl_gene_id])
vp <- v[rownames(v)%in%primary_genes,,drop=FALSE]
pca_primary <- prcomp(t(vp[head(order(rowVars(vp),decreasing=TRUE),2000),]),center=TRUE,scale.=FALSE)
primary <- data.table(PC=paste0('PC',1:3),
    all_loci_var_pct=100*newpc$sdev[1:3]^2/sum(newpc$sdev^2),
    primary_loci_var_pct=100*pca_primary$sdev[1:3]^2/sum(pca_primary$sdev^2),
    absolute_score_correlation=vapply(1:3,function(i)abs(cor(newpc$x[,i],pca_primary$x[,i])),numeric(1)),
    scope='sensitivity using existing VST; reselect top2000 variable primary-assembly genes; not raw-read re-quantification')
fwrite(primary,file.path(base_dir,'rna_primary_assembly_pca_sensitivity.csv'))
print(summary);print(top[1:10]);print(primary)
# Independent implementation of representative-level correlation tests: compare
# deposited rho/p with stats::cor.test on the underlying matched abundance pairs.
pm <- read_matrix('output/prot_abundance_matrix.csv')
fm <- fread('metadata/line_family_map.csv')
rep_lines <- intersect(fm[patient_representative==TRUE,cell_line],intersect(colnames(ns),colnames(pm)))
ct <- fread('output/integ_rnaprot_patientrep_cor.csv')
set.seed(90205); genes <- unique(c(ct[order(spearman)][1:10,gene],ct[order(-spearman)][1:10,gene],sample(ct$gene,100)))
checks <- rbindlist(lapply(genes,function(g){
  ok <- is.finite(ns[g,rep_lines]) & is.finite(pm[g,rep_lines]); cl <- rep_lines[ok]
  z <- suppressWarnings(cor.test(ns[g,cl],pm[g,cl],method='spearman',exact=FALSE))
  oldrow <- ct[gene==g]
  data.table(gene=g,n=length(cl),rho_error=abs(unname(z$estimate)-oldrow$spearman),p_error=abs(z$p.value-oldrow$p_value))
}))
stopifnot(max(checks$rho_error)<1e-12,max(checks$p_error)<1e-12)
fwrite(checks,file.path(base_dir,'rna_patient_correlation_verification.csv'))
print(data.table(n_genes_checked=nrow(checks),max_rho_error=max(checks$rho_error),max_p_error=max(checks$p_error)))
