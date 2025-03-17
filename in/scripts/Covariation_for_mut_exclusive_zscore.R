### mutually exclusive subunits zscores ### 
library(data.table);library(stringr); library(ggplot2);library(WGCNA) 
#### functions used #### 
zscoring = function(x){
  (x-mean(x,na.rm =T ))/sd(x,na.rm=T)
}
f_BIC<- function( input_data, value_name,min_n ){
  missingness_tmp <- base::crossprod(!is.na(input_data))>min_n # with microproteins being very rare it's common to only have 1 point in common produce fake correlations
  #using suggested bicor
  tmp <- bicor( input_data , use = "pairwise.complete.obs", nThreads = 6)                     # Robust correlation using all pairwise complete observations
  tmp[missingness_tmp == FALSE] <- NaN
  tmp <- as.data.table( reshape2::melt( as.matrix( tmp  )))                                   # Turn distance matrix into a pair-wise data.table
  tmp <- tmp[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), value ) ]   # Rename and change to character
  tmp <- tmp[ Protein_1 > Protein_2 ]                                                         # Remove redundant pairs by keeping only A > B, removing A == B and B < A pairs
  names(tmp)[3] <- value_name                                                                 # Assign new value name
  return(tmp)
}
convert_to_gene <- function(x){
  gene_names = c()
  for(i in x){
    Uniprots = strsplit(i,';') |> unlist() |> str_remove('-.$')
    Genes = Uniprot_mapping[Type == 'Gene_Name' & Protein %in% Uniprots, ID] |> unique() |> paste(collapse = ';')
    gene_names[i] = Genes
  }
  gene_names
}
calc_pairwise_overlaps <- function(sets) {
  # Ensure that all sets are unique character vectors
  sets_are_vectors <- vapply(sets, is.vector, logical(1))
  if (any(!sets_are_vectors)) {
    stop("Sets must be vectors")
  }
  sets_are_atomic <- vapply(sets, is.atomic, logical(1))
  if (any(!sets_are_atomic)) {
    stop("Sets must be atomic vectors, i.e. not lists")
  }
  sets <- lapply(sets, as.character)
  is_unique <- function(x) length(unique(x)) == length(x)
  sets_are_unique <- vapply(sets, is_unique, logical(1))
  if (any(!sets_are_unique)) {
    stop("Sets must be unique, i.e. no duplicated elements")
  }
  
  n_sets <- length(sets)
  set_names <- names(sets)
  n_overlaps <- choose(n = n_sets, k = 2)
  
  vec_name1 <- character(length = n_overlaps)
  vec_name2 <- character(length = n_overlaps)
  vec_num_shared <- integer(length = n_overlaps)
  vec_overlap <- numeric(length = n_overlaps)
  vec_jaccard <- numeric(length = n_overlaps)
  overlaps_index <- 1
  
  for (i in seq_len(n_sets - 1)) {
    name1 <- set_names[i]
    set1 <- sets[[i]]
    for (j in seq(i + 1, n_sets)) {
      name2 <- set_names[j]
      set2 <- sets[[j]]
      
      set_intersect <- set1[match(set2, set1, 0L)]
      set_union <- .Internal(unique(c(set1, set2), incomparables = FALSE,
                                    fromLast = FALSE, nmax = NA))
      num_shared <- length(set_intersect)
      overlap <- num_shared / min(length(set1), length(set2))
      jaccard <- num_shared / length(set_union)
      
      vec_name1[overlaps_index] <- name1
      vec_name2[overlaps_index] <- name2
      vec_num_shared[overlaps_index] <- num_shared
      vec_overlap[overlaps_index] <- overlap
      vec_jaccard[overlaps_index] <- jaccard
      
      overlaps_index <- overlaps_index + 1
    }
  }
  
  result <- data.frame(name1 = vec_name1,
                       name2 = vec_name2,
                       num_shared = vec_num_shared,
                       overlap = vec_overlap,
                       jaccard = vec_jaccard,
                       stringsAsFactors = FALSE)
  return(result)
}
#### Uniprot annotation file #### 
Uniprot_mapping = fread(here::here('in','datasets','HUMAN_9606_idmapping.dat.gz'), header = F)
setnames(Uniprot_mapping,c('Protein','Type','ID'))
Gene_names = Uniprot_mapping[,N_ids := .N, by =Protein][Type == 'Gene_Name']
Gene_names = Gene_names[order(ID,N_ids)][, tail(.SD, 1), by = ID]
Gene_names = Gene_names[,.(ID,Protein)][, tail(.SD, 1), by = Protein]

# paralogs = readxl::read_xlsx(here::here('in','datasets','NIHMS1737332-supplement-2.xlsx'), skip = 1) |> as.data.table()
# paralogs = paralogs[target2 != 'NA'][target1 != 'NA'][,.(target1,target2)]  |> unique()
# setnames(paralogs,c('ID.x','ID.y'))
# paralogs_rev = copy(paralogs)
# paralogs_rev[,`:=`(ID.x = ID.y,
#                    ID.y = ID.x)]
# paralogs = rbind(paralogs,paralogs_rev)
# paralogs[,is_paralog := T]
paralog_egg = fread(here::here('in','datasets','uniprotkb_taxonomy_id_9606_AND_reviewed_mapEGGNOG_2024_10_07.tsv'))
paralog_egg[,N_eggNOG:= .N, by = eggNOG]
paralog_egg= paralog_egg[between(N_eggNOG,2,10)]
paralog_pairs = data.table()
for(i in unique(paralog_egg$eggNOG)){
  print(i)
  paralog_pairs_tmp = paralog_egg[eggNOG == i,Entry] |> 
    combn(2) |> t() |> as.data.table() 
  paralog_pairs_tmp[,eggNOG:=i]
  paralog_pairs = rbind(paralog_pairs_tmp,paralog_pairs)
  }
setnames(paralog_pairs,c('Protein_1','Protein_2',''))

paralog_pairs[Protein_1<Protein_2,`:=`(Protein_1=Protein_2,
                                  Protein_2 = Protein_1)]
paralog_pairs[Protein_1<Protein_2]

### readhumap3 clusters & define clusters to use #### 
humap_3 = fread(here::here('in','datasets','Supplemental_Table_3_hu.MAP3.1_complexes_wConfidenceScores_total15326_20240922.csv'),header = T,  sep = ',')
humap_3[,N_members := stringr::str_count(uniprotACCs,' ')+1, by = clustID]
clusters_to_check = humap_3[cluster_confidence<4 & N_members>3]
prots_per_complex = clusters_to_check[,.(clustID,uniprotACCs,cluster_confidence)] |> 
  tidyr::separate_longer_delim(uniprotACCs,' ') |> as.data.table()
pairs_to_clust = data.table()
for(i in clusters_to_check$clustID){
  print(i)
  clust_tmp = combn(strsplit(clusters_to_check[clustID == i,uniprotACCs],' ') |> unlist(),2) |> 
    t() |> as.data.table()
  pairs_to_clust = rbind(pairs_to_clust,clust_tmp[,clustID := i])
                         
}
setnames(pairs_to_clust,c('V1','V2'),c('Protein_1','Protein_2'))
pairs_to_clust[Protein_1<Protein_2,`:=`(Protein_1=Protein_2,
                                  Protein_2 = Protein_1)]
pairs_to_clust[Protein_1<Protein_2]
fwrite(pairs_to_clust,here::here('out','datasets','pairs_to_clust.gz'))
pairs_to_clust = fread(here::here('out','datasets','pairs_to_clust.gz'))
### read Sam mut_excl sets #### 
mut_excl_old = fread(here::here('in','datasets','refiltered_mutually_exclusive_and_structurally_consistent_protein_pairs_03OCT2024.tsv') )
mut_excl_old =mut_excl_old[,.(protein_1,protein_2,cmmn_p_homodimer)]
mut_excl = fread(here::here('in','datasets','MutEx_pairs_refiltered_w_large_comps_removed_04OCT2024.tsv'))
new_pairs = fread(here::here('in','datasets',
                             'TableS4_structurally_consistent_mutually_exclusive_modeled_pairs_w_janes_w_burke_12MAR2024.csv'))
new_pairs = new_pairs[,.(pair_1,pair_2,interface,common_protein,protein_1,protein_2)]
setnames(new_pairs,c('pair_1','pair_2','interface'),c('V1','V2','class'))
                 # '/mnt/kustatscher/members/savvas/LFQ_DIA_SWATH/mutually_exclusive_dimer_analysis_10SEP2024.csv')
# mut_excl = fread(here::here('in','datasets','humap3_mutual_exclusion_vs_non_dimers_annotated_homodimers_30JUL2024.csv'))
mut_excl = merge(mut_excl,mut_excl_old,c('protein_1','protein_2'))
mut_excl  =mut_excl[interface_overlap=='yes',.(V1,V2,common_protein,protein_1,protein_2,cmmn_p_homodimer)]
mut_excl[,Class:='mut_excl']
mut_excl = new_pairs
setnames(mut_excl,c('protein_1','protein_2'),c('Protein_1','Protein_2'))
mut_excl[Protein_1<Protein_2,`:=`(Protein_1=Protein_2,
                                  Protein_2 = Protein_1)]
mut_excl[Protein_1<Protein_2]

### read Procan and normalise #### 
# to_load = list.files('/home/v1skourt/Downloads/drug_peptide',full.names = T)
mapping_1 = fread(here::here('in','datasets','ProCan-DepMapSanger_mapping_file_averaged.txt'))
tissue_n = mapping_1[,.N, by = Cancer_type]

ProCan = readxl::read_xlsx(here::here('in','datasets','1-s2.0-S1535610822002744-mmc3.xlsx'), sheet = 2, skip = 1)
ProCan_matrix = ProCan[,-1] |> as.matrix()
rownames(ProCan_matrix) =  ProCan[,1] |> unlist()
colnames(ProCan_matrix) = colnames(ProCan_matrix) |> str_remove(';[:print:]*$')

ProCan_matrix[1:20,] |> t() |> boxplot()
tmp_medians <- apply( ProCan_matrix , 1, median, na.rm = TRUE )  
ProCan_matrix_norm <- sweep( ProCan_matrix , 1, tmp_medians, FUN = "-" )
ProCan_matrix_norm[1:20,] |> t() |> boxplot()
ProCan_matrix_norm_piv = ProCan_matrix_norm |> t() |> as.data.frame() |> tibble::rownames_to_column('Uniprot') |> as.data.table() |> 
  melt(id.vars = 'Uniprot', variable.name = 'Project_Identifier',value.name = 'Abundance')
ProCan_matrix_norm_piv = merge(ProCan_matrix_norm_piv[!is.na(Abundance)],mapping_1[Cancer_type %in% tissue_n[N>20,Cancer_type],
                                                                .(Project_Identifier,Cancer_type)], by = 'Project_Identifier')
ProCan_matrix_norm_piv[,N_tissue := .N, by = .(Uniprot,Cancer_type )]
ProCan_matrix_norm_piv = ProCan_matrix_norm_piv[N_tissue>20]
ProCan_matrix_norm_piv_plot = ProCan_matrix_norm_piv[,.(Prot_SD = sd(Abundance,na.rm = T)),.(Cancer_type,Uniprot )
                                                     ][,median_protSD := median(Prot_SD,na.rm =T),by = Cancer_type]
  ggplot(ProCan_matrix_norm_piv_plot,aes(x = Prot_SD, y = reorder(Cancer_type,median_protSD) ))+
  geom_boxplot()+theme_bw()+
  geom_vline(xintercept = median(unique(ProCan_matrix_norm_piv_plot$median_protSD) ),colour  = 'red')+
  labs(x= 'across cell line protein SD',y = 'cancer lineage')+
  ggtitle('Proteins have similar variation excell B-Cell')
ggsave(here::here('out','plots','Procan_protein_SD.pdf'))

ProCan_matrix_norm_piv_plot_subset =ProCan_matrix_norm_piv[,head(.SD,20),by = .(Uniprot,Cancer_type) 
                                                           ][,.(Prot_SD = sd(Abundance)),.(Cancer_type,Uniprot )
                                                             ][,median_protSD := median(Prot_SD,na.rm =T),by = Cancer_type]
  ggplot(ProCan_matrix_norm_piv_plot_subset,aes(x =Prot_SD,y = reorder(Cancer_type,median_protSD)  ))+
  geom_boxplot()+theme_bw()+
  geom_vline(xintercept = median(unique(ProCan_matrix_norm_piv_plot$median_protSD) ),colour  = 'red')+
  labs(x= 'across 20 cell line protein SD',y = 'cancer lineage')+
  ggtitle('Tissues subsetted to 20')
ggsave(here::here('out','plots','Procan_protein_SD_subset.pdf'))

### Global bicor #### 
Pairwise_procan = f_BIC(ProCan_matrix_norm,'bicor_Procan_across_tissue',30)
# fwrite(Pairwise_procan,here::here('out','datasets','ProCAN_global_bicor_min_30.gz'))
Pairwise_procan = fread(here::here('out','datasets','ProCAN_global_bicor_min_30.gz'))

pairs_to_clust_mutexc = merge(pairs_to_clust,mut_excl, by = c('Protein_1','Protein_2'),all.x = T)
pairs_to_clust_mutexc[,cmmn_p_homodimer:='non_homodimer']
Bicor_pairs = merge(Pairwise_procan,pairs_to_clust_mutexc,by =c('Protein_1','Protein_2'), all.x = T)
# Bicor_pairs[,cmmn_p_homodimer:=]
Bicor_pairs[,Type:= fifelse(is.na(cmmn_p_homodimer),as.character(!is.na(clustID)),cmmn_p_homodimer)]
Bicor_pairs[Protein_1 =='Q8WVM7'& Protein_2 == 'Q8N3U4']

Bicor_pairs[,Type:= fcase(
  Type =='no','interphase overlap',
  Type =='yes','interphase overlap (homodimer)',
  Type =='TRUE','subunit pairs',
  Type =='FALSE','non-subunits pairs'
) ]
types_tmp = unique(Bicor_pairs$Type)
subunit_colour = c('interphase overlap'  ='darkblue',
                   'interphase overlap (homodimer)'= '#00A36C',
                   'subunit pairs' = '#A7C7E7',
                   'non-subunits pairs' = 'grey70')
names(subunit_colour) = types_tmp[c(3,4,2,1)]
fwrite(Bicor_pairs,here::here('out','datasets','procan_bicor_mut_exclusive.gz'))
Bicor_pairs = fread(here::here('out','datasets','procan_bicor_mut_exclusive.gz'))
  ggplot(Bicor_pairs,aes(x= bicor_Procan_across_tissue, fill =Type))+
  geom_density(alpha = 0.5)+theme_bw()+
    labs(x = 'Across ProCan Covariation', y = 'Density')+
    scale_fill_manual(values = subunit_colour)+
  ggtitle('Protein Pairs in the same complex have higher covariation',
          subtitle = 'homodimers have the highest covariation')
  ggsave(here::here('out','plots','Procan_bicor_mut_excpairs.pdf'))
### tissue specific bicor #### 
tissue_n = mapping_1[,.N, by = Cancer_type]

tissue_cor = data.table()
# 20 celllines gives us 18 tissues
tissue_n[N>20]
for(i in unique(mapping_1$Cancer_type)){
  print(i)
  cell_lines_tmp = mapping_1[Cancer_type == i,Project_Identifier]
  if(length(cell_lines_tmp)>20){
    ProCan_matrix_complex_tissue = ProCan_matrix_norm[rownames(ProCan_matrix_norm) %in% cell_lines_tmp,]
   # ProCan_matrix_complex_tissue[,'Q6STE5'] |> is.na() |> table() |> print()
    # }}
    # min 15 shared datapoints between 2 proteins for bicor
    cor_tmp =  f_BIC(ProCan_matrix_complex_tissue,'bicor',20)
    cor_tmp[,tissue:=i]
    tissue_cor = rbind(tissue_cor[!is.na(bicor)],
                       cor_tmp)
  }
}
# tissue_cor[,global_average := median(bicor,na.rm =T), by = .(Protein_1,Protein_2)]
fwrite(tissue_cor,here::here('out','datasets','procan_covariation_per_tissue.gz'))
tissue_cor = fread(here::here('out','datasets','procan_covariation_per_tissue.gz'))
tissue_cor_subset = data.table()
# 20 celllines gives us 18 tissues
tissue_n[N>20]
for(i in unique(mapping_1$Cancer_type)){
  print(i)
  cell_lines_tmp = mapping_1[Cancer_type == i,Project_Identifier]
  if(length(cell_lines_tmp)>20){
    ProCan_matrix_complex_tissue = ProCan_matrix_norm[rownames(ProCan_matrix_norm) %in% head(cell_lines_tmp,21),]
    
    # ProCan_matrix_complex_tissue[,'Q6STE5'] |> is.na() |> table() |> print()
    # }}
    # min 15 shared datapoints between 2 proteins for bicor
    cor_tmp =  f_BIC(ProCan_matrix_complex_tissue,'bicor',20)
    cor_tmp[,tissue:=i]
    tissue_cor_subset = rbind(tissue_cor_subset[!is.na(bicor)],
                       cor_tmp)
  }
}
# tissue_cor[,global_average := median(bicor,na.rm =T), by = .(Protein_1,Protein_2)]
fwrite(tissue_cor_subset,here::here('out','datasets','procan_covariation_per_tissue_subset.gz'))
tissue_cor_subset = fread(here::here('out','datasets','procan_covariation_per_tissue_subset.gz'))
setnames(tissue_cor_subset,'bicor','subset_20_bicor')
# tissue_cor = fread(here::here('out','datasets','procan_covariation_per_tissue.gz'))
tissue_cor = tissue_cor[!is.na(bicor)]
tissue_cor_subset = tissue_cor_subset[!is.na(subset_20_bicor)]
tissue_n |> ggplot(aes(x = N,y = Cancer_type))+
  geom_col()+
  labs(x= 'Number of cell lines')

subset_tissue_comparison = merge(tissue_cor_subset,tissue_cor,by =c('Protein_1','Protein_2','tissue'))
subset_tissue_comparison[,N_tissues := .N, by = c('Protein_1','Protein_2')]
subset_tissue_comparison = subset_tissue_comparison[N_tissues>15]
subset_tissue_comparison[,`:=`(zscore_subset = zscoring(subset_20_bicor),
                               zscore_all =zscoring(bicor) ),by = .(Protein_1,Protein_2)]
subset_tissue_comparison[,tissue:= factor(tissue,levels = tissue_n[N>20][order(N),Cancer_type])]
subset_tissue_comparison |> ggplot(aes(x= zscore_all,y = zscore_subset))+
  geom_hex(bins = 25)+ theme_bw()+
  facet_wrap('tissue')+
  stat_poly_eq() +
  scale_fill_continuous(type = "viridis") + 
  geom_abline(intercept = 0,slope =1 ,colour = 'red')+
  labs(x = 'zscore including all cell lines',y = 'zscore using 20 cell lines per tissue')+
  ggtitle('Subset 20 zscoring of pairwise bicor',
          'ordered by total sample size')
ggsave(here::here('out','plots','procan_bicor_subset_vs_all_zscore.pdf'),width = 13,height = 10)

subset_tissue_comparison |> ggplot(aes(x= bicor,y = subset_20_bicor))+
  geom_hex(bins = 25)+ theme_bw()+
  facet_wrap('tissue')+
  stat_poly_eq() +
  scale_fill_continuous(type = "viridis") + 
  geom_abline(intercept = 0,slope =1 ,colour = 'red')+
  labs(x = 'bicor including all cell lines',y = 'bicor using 20 cell lines per tissue')+
  ggtitle('Subset 20 zscoring of pairwise bicor',
          'ordered by total sample size')
ggsave(here::here('out','plots','procan_bicor_subset_vs_all.pdf'),width = 13,height = 10)

tissue_cor[,N_tissues := .N, by = c('Protein_1','Protein_2')]
tissue_cor_mutexcl = merge(tissue_cor[N_tissues>15],pairs_to_clust_mutexc,by =c('Protein_1','Protein_2'), all.x = T, allow.cartesian = T)
tissue_cor_mutexcl[,Type:= fifelse(is.na(cmmn_p_homodimer),as.character(!is.na(clustID)),cmmn_p_homodimer)]
tissue_cor_mutexcl[,Type:= fcase(
  Type =='no','interphase overlap',
  Type =='yes','interphase overlap (homodimer)',
  Type =='TRUE','subunit pairs',
  Type =='FALSE','non-subunits pairs'
) ]
ggplot(tissue_cor_mutexcl,aes(y = tissue, x = bicor ,fill = Type))+
  geom_boxplot()+theme_bw()+
  facet_wrap('Type')+
  labs(x ='lineage-specific covariation', y = 'cancer lineage')+   scale_fill_manual(values = subunit_colour)+
  ggtitle('ProCan lineage specific protein pairwise Covariation')
ggsave(here::here('out','plots','Procan_lineage_bicor_mut_excpairs.pdf'),width = 11,height = 9)
fwrite(tissue_cor_mutexcl,here::here('out','datasets','Procan_lineage_bicor_mut_excpairs.gz'))
### tissue specific complex bicor #### 
tissue_complex_cor = data.table()
proteins_correlations =  unique(c(tissue_cor$Protein_1,tissue_cor$Protein_2))
complex_id = 'huMAP3_00356.1'
Protein_tmp = 'P28062'

for(complex_id in clusters_to_check[, clustID]){
  print(complex_id)
  
  partners_tmp = merge(tissue_cor,pairs_to_clust[clustID == complex_id], by = c('Protein_1','Protein_2'))
  partners_tmp  =partners_tmp[!is.na(bicor)]
  if(nrow( partners_tmp)>15){
  partners_tmp[,N_pair_per_tissue:=.N, by =tissue]
 
  partners_tmp  =partners_tmp[N_pair_per_tissue>3]
  # in order to calculate zscore per pair i need to see the same pair in multiple tissues
  partners_tmp[,N_tissue_per_pair :=.N, by =.(Protein_1,Protein_2)]
  partners_tmp = partners_tmp[N_tissue_per_pair>5]
  if(nrow( partners_tmp)>15){
  # partners_tmp = partners_tmp[partners_tmp != Protein_tmp]
  # test_tissue = "Esophageal Carcinoma"
  # # # 
  # DTm <- dcast( data = rbind( tissue_cor_tmp[tissue == test_tissue, .(Protein_1, Protein_2, bicor)],                             # These steps create a "redundant" table...
  #                             tissue_cor_tmp[tissue ==test_tissue, .(Protein_1 = Protein_2, Protein_2 = Protein_1, bicor)]),    # ... containing both A <-> B and B <-> pairs
  #               Protein_1 ~ Protein_2 , value.var = "bicor")                      # And this casts them into a matrix-shaped data.table
  # DTm <- as.data.frame(DTm)                 # Turn into data.frameconv
  # rownames(DTm) <- DTm$Protein_1            # Create rownames
  # DTm$Protein_1 <- NULL                     # Drop original name column
  # DTm <- as.matrix(DTm)
  # diag(DTm) = 1
  # rownames_tmp_cor_dcast = data.table(IDs = rownames(DTm) |> convert_to_gene(),
  #                                     Uniprot = rownames(DTm))
  # rownames_tmp_cor_dcast[,isoform:= str_match(Uniprot,'-.') ]
  # rownames_tmp_cor_dcast[!is.na(isoform),IDs:=paste0(IDs,isoform,collapse = ''), by = Uniprot]
  # rownames(DTm) = rownames_tmp_cor_dcast$IDs
  # # pdf(here::here('out','plots',glue::glue('Procan_tissue_bicor_{test_tissue}_{complex_id}.pdf')))
  # pheatmap::pheatmap(DTm, main = paste(test_tissue,complex_id))
  # # dev.off()
  
    members = unique(c(partners_tmp$Protein_1,partners_tmp$Protein_2))
    # this makes the tissues comparable
    partners_tmp[,z_score_tissues :=zscoring(bicor),by = .(tissue)]
    # this helps ask in which tissue does the pair have its best covariation
    partners_tmp[,z_score_pairs :=zscoring(z_score_tissues),by = .(Protein_1,Protein_2)]
    # all the pairs do worse in the tissue that does worse, of course, so need to normalise tissue first
    # partners_tmp[,z_score_only_pairs :=zscoring(bicor),by = .(Protein_1,Protein_2)]
    # partners_tmp[Protein_1 ==  Protein_tmp | 
    #                Protein_2 == Protein_tmp] |> ggplot(aes(x=tissue,y = z_score_only_pairs))+
    #   geom_boxplot()+coord_flip()
    # partners_tmp |> ggplot(aes(x=tissue,y = z_score_only_pairs))+
    #   geom_boxplot()+coord_flip()
    core_complex = partners_tmp[,.(pair_avg = mean(z_score_tissues,na.rm = T),
                    pair_SD = sd(z_score_tissues)),by = .(Protein_1,Protein_2)]
    # core_complex[,`:=`(IDx = convert_to_gene(Protein_1),
    #                    IDy = convert_to_gene(Protein_2)), by =.(Protein_1,Protein_2)]
    # core_complex |> ggplot(aes(x = pair_SD, label = paste(IDx,IDy,sep = '-'),
    #                            y = pair_avg))+
    #   geom_point()+
    #   ggrepel::geom_text_repel()+theme_bw()+
    #   ggtitle('identifying -core-complex subunits')
  if(length(members)>10){
    core_complex = core_complex[pair_avg>0.4][order(pair_SD)]
    core_complex = core_complex[,head(.SD,nrow(core_complex)/2)] }
    core_members = core_complex[,.(Protein_1,Protein_2)] |> unlist() |> unique()
    for(Protein_tmp in members){
      tissue_cor_tmp_member = partners_tmp[Protein_1 ==  Protein_tmp | 
                                              Protein_2 == Protein_tmp]
      core_members_tmp= core_members[!(core_members == Protein_tmp)]
      core_complex_cor = tissue_cor_tmp_member[Protein_1 %in% core_members_tmp|
                                                 Protein_2 %in% core_members_tmp
                                               ][,N_tissue:= .N, by = tissue
                                                 ][N_tissue>2]
      tissue_cor_tmp_member = tissue_cor_tmp_member[,N_tissue:= .N, by = tissue
                                                    ][N_tissue>3]
      # if(ncol(complex_procan)>5){(Protein_tmp)
      if(nrow( core_complex_cor)>3){
        complex_cor = tissue_cor_tmp_member[,.(complex_cor = mean(bicor,na.rm = T),
                                               # complex_zscore_tissue = mean(z_score_tissues,na.rm = T),
                                               complex_zscore_pair = mean(z_score_pairs,na.rm = T),
                                                           N_subunits = .N), by = tissue] 
        core_complex_cor = core_complex_cor[,.(complex_core_cor = mean(bicor,na.rm = T),
                                               # complex_core_zscore_tissue = mean(z_score_tissues,na.rm = T),
                                               complex_core_zscore_pair = mean(z_score_pairs,na.rm = T),
                                               N_core_subunits = .N), by = tissue] 
        complex_cor  =merge(complex_cor,core_complex_cor, by = 'tissue',all = T)
        # complex_cor |> ggplot(aes(x = complex_core_cor,label = tissue,
        #                           y = complex_core_zscore_pair))+
        #   geom_point()+
        #   ggrepel::geom_text_repel()+
        #   ggtitle('PSMB8 CORE subunits complex correlation ',
        #           subtitle = 'bicor avg (x), zscore tissue and protein pairs')
        complex_cor[,`:=`(subunit =Protein_tmp ,
                          complex = complex_id)]
        tissue_complex_cor = rbind(tissue_complex_cor,
                                   complex_cor)
      }
    }
  }
  }
}
tissue_complex_cor = tissue_complex_cor[is.finite(complex_cor)] |> unique()
tissue_complex_cor[complex =='huMAP3_00356.1'] |> 
  ggplot(aes(x = complex_zscore_pair, y = complex_core_zscore_pair))+
  geom_point()
tissue_complex_cor |> 
  ggplot(aes(x = complex_zscore_pair, y = complex_core_zscore_pair))+
  geom_hex()
fwrite(tissue_complex_cor,'per_tissue_subunit_cor.gz')
tissue_complex_cor = fread('per_tissue_subunit_cor.gz')

tissue_complex_cor = tissue_complex_cor[!is.na(bicor)]
tissue_complex_cor[,N_tissues := .N, by = c('subunit','complex')]
ggplot(tissue_complex_cor[N_tissues>15],aes(y = factor(tissue,levels = tissue_n[N>20][order(N),Cancer_type]),
                                            x = complex_core_cor))+
  geom_boxplot()+theme_bw()+
  labs(x ='lineage-specific subunit-complex covariation', y = 'cancer lineage')+ 
  ggtitle('ProCan lineage subunit to complex average Covariation')
ggsave(here::here('out','plots','Procan_lineage_bicor_subunit_complex.pdf'))

ggplot(tissue_complex_cor[N_tissues>15],aes(y = factor(tissue,levels = tissue_n[N>20][order(N),Cancer_type]),
                                            x = complex_core_zscore_pair ))+
  geom_boxplot()+theme_bw()+
  labs(x ='lineage-specific subunit-complex zscore covariation', y = 'cancer lineage')+ 
  ggtitle('ProCan lineage subunit to complex zscore Covariation')
ggsave(here::here('out','plots','Procan_lineage_bicor_subunit_zscorecomplex.pdf'))
fwrite(tissue_complex_cor,here::here('out','datasets','Procan_lineage_bicor_subunit_complex.gz'))


max_subunit_cor = tissue_complex_cor[order(complex_cor,subunit,complex)][,tail(.SD,1) ,by = .(subunit,complex)]

tissue_complex_cor = merge(tissue_complex_cor,Gene_names,by.x = 'subunit',by.y = 'Protein',all.x=T)

### by-tissue complex bicor #### 
pair_intertissue_complex_cor = data.table()
complex_tmp = 'huMAP3_00356.1'
for(complex_tmp in clusters_to_check[, clustID]){
  print(complex_tmp)
  subunit_complex_tmp = tissue_complex_cor[complex ==complex_tmp ]
  if(length(unique(subunit_complex_tmp$tissue))>4){
    # other_members = tissue_complex_cor[complex==complex_tmp ]
    # dcast(subunit_complex_tmp,tissue ~subunit,value.var = 'complex_core_zscore_pair' ) |>
    #   ggplot(aes(x = P25788,label=tissue,  y = P25786))+
    #   geom_point()+theme_bw()+
    #   labs(x = 'PSMA3',y = 'PSMA1')+
    #   ggrepel::geom_text_repel()
    # testing  = dcast(subunit_complex_tmp,tissue ~subunit,value.var = 'complex_cor' ) |> tibble::column_to_rownames('tissue') |> as.matrix()
    # pheatmap(testing,cluster_rows = F,cluster_cols = F)
    complex_cor_tmp = dcast(subunit_complex_tmp,tissue ~subunit,value.var = 'complex_core_zscore_pair' ) |> tibble::column_to_rownames('tissue') |> 
      f_BIC('tissue_cor',8)
    complex_cor_tmp[,clustID:=complex_tmp]
    pair_intertissue_complex_cor = rbind(pair_intertissue_complex_cor,
                          complex_cor_tmp[!is.na(tissue_cor)])
    
  }
}
fwrite(pair_intertissue_complex_cor,here::here('out','datasets','mutexl_pairs_intertissue_complex_cor.csv'))
pair_intertissue_complex_cor = fread(here::here('out','datasets','mutexl_pairs_intertissue_complex_cor.csv'))
intertissue_pairs = merge(pair_intertissue_complex_cor,unique(pairs_to_clust_mutexc[,.(Protein_1,Protein_2,cmmn_p_homodimer)]),by =c('Protein_1','Protein_2'), all.x = T)
intertissue_pairs[,Type:= fifelse(is.na(cmmn_p_homodimer),as.character(!is.na(clustID)),cmmn_p_homodimer)]

intertissue_pairs[,Type:= fcase(
  Type =='no','interphase overlap',
  Type =='yes','interphase overlap (homodimer)',
  Type =='TRUE','subunit pairs',
  Type =='FALSE','non-subunits pairs'
) ]
merge(intertissue_pairs,mut_excl,by =c('Protein_1','Protein_2'), all.x = T) |> 
  ggplot(aes(x= tissue_cor, fill = Type))+
  geom_density(alpha = 0.5)+theme_bw()+
  scale_fill_manual(values = subunit_colour)+
  labs(x = 'Complex Inclusion',y = 'Density')+
  ggtitle('Protein mut_excl Pairs with w/o homodimer',
          subtitle = 'high values, two proteins are correlate well with the complex in the same samples')
ggsave(here::here('out','plots','Procan_perlineage_bicor.pdf'))
### Across Tissue average bicor ####
ProCan_matrix_norm_centered = ProCan_matrix_norm

ProCan_matrix_norm_centered[,1:20] |> boxplot()
# centering each protein
tmp_medians <- apply( ProCan_matrix_norm_centered , 2, median, na.rm = TRUE )  
ProCan_matrix_norm_centered <- sweep( ProCan_matrix_norm_centered , 2,tmp_medians, FUN = "-" )
ProCan_matrix_norm_centered[,1:20] |> boxplot()
ProCan_matrix_norm_centered[1:20,] |> t() |> boxplot()

ProCan_tissue = ProCan_matrix_norm_centered |> t() |> as.data.frame() |> tibble::rownames_to_column('Uniprot') |> as.data.table() |> 
  melt(id.vars = 'Uniprot', variable.name = 'Project_Identifier',value.name = 'LogRatio')
ProCan_tissue = merge(ProCan_tissue,mapping_1[Cancer_type %in% tissue_n[N>3,Cancer_type],.(Project_Identifier,Cancer_type)], by = 'Project_Identifier')
ProCan_tissue =ProCan_tissue[,.(Prot_avg = mean(LogRatio,na.rm = T)),by = .(Uniprot,Cancer_type)]
ProCan_tissue =ProCan_tissue |> dcast(Cancer_type~Uniprot,value.var = 'Prot_avg') |> 
  tibble::column_to_rownames('Cancer_type')
Tissue_level_covar = f_BIC(ProCan_tissue,'tissue_level_covar',20)
fwrite(Tissue_level_covar,here::here('out','datasets','Tissue_level_covar.gz'))
Tissue_level_covar= fread(here::here('out','datasets','Tissue_level_covar.gz'))
tissue_level_pairs = merge(Tissue_level_covar,unique(pairs_to_clust_mutexc[,.(Protein_1,clustID,Protein_2,cmmn_p_homodimer)]),by =c('Protein_1','Protein_2'), all.x = T)
tissue_level_pairs[,Type:= fifelse(is.na(cmmn_p_homodimer),as.character(!is.na(clustID)),cmmn_p_homodimer)]

tissue_level_pairs[,Type:= fcase(
  Type =='no','interphase overlap',
  Type =='yes','interphase overlap (homodimer)',
  Type =='TRUE','subunit pairs',
  Type =='FALSE','non-subunits pairs'
) ]
tissue_level_pairs |> 
  ggplot(aes(x= tissue_level_covar, fill = Type))+
  geom_density(alpha = 0.5)+theme_bw()+
  scale_fill_manual(values = subunit_colour)+
  labs(x = 'tissue_level_covar',y = 'Density')+
  ggtitle('covaration for tissue-level')
ggsave(here::here('out','plots','Procan_perlineage_bicor.pdf'))
### Procan Imputation and score #### 
complex_tmp = 'huMAP3_00356.1'
tissue_complex_cor[subunit %in% prots_per_complex[clustID ==complex_tmp,uniprotACCs]] |> 
  View()
tissue_tmp = mapping_1[Cancer_type == 'Non-Small Cell Lung Carcinoma',Project_Identifier]
Procan_lineplot  = ProCan_matrix_norm_centered[tissue_tmp,]  
complex_presence = Procan_lineplot[,colnames(Procan_lineplot) %in% 
                                     c('P28074','P28062')]
samples_no_NA = complex_presence[!matrixStats::rowAnyNAs(complex_presence),] |> rownames()
Procan_lineplot = Procan_lineplot[samples_no_NA,] |> t() |> as.data.frame() |> tibble::rownames_to_column('Uniprot') |> as.data.table() |> 
  melt(id.vars = 'Uniprot', variable.name = 'Cell_line',value.name = 'LogRatio')
Procan_lineplot[,ID:=convert_to_gene(Uniprot),by = Uniprot]
Procan_lineplot[,category:=dplyr::case_when(
  ID =='PSMB5'~'PSMB5',
  ID =='PSMB8'~'PSMB8',
  Uniprot %in% prots_per_complex[clustID ==complex_tmp,uniprotACCs]~'Proteasome subunits',
  TRUE~'other'
  
)]
colour_pallet = c('grey95','#ffcccb','darkblue','darkred')
names(colour_pallet) = c('other','Proteasome subunits','PSMB8','PSMB5')
Procan_lineplot[] |> ggplot(aes(x = reorder(Cell_line,LogRatio), y= LogRatio,colour= category,group=Uniprot))+
  geom_line(data = Procan_lineplot[category == 'other'],alpha = 0.05)+
  geom_line(data = Procan_lineplot[category != 'other'],alpha = 0.9)+
  geom_line(data = Procan_lineplot[category %in% c('PSMB5','PSMB8')],alpha = 1)+
  # geom_line(data = Procan_lineplot[category == Protein_OI],alpha = 1)+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_colour_manual(values = colour_pallet)+
  theme(text = element_text(size = 15))+
  labs(x = 'Cancer Cell Lines',y = 'Protein Normalised Abundance')
# annotate('text',label = glue::glue('adj.pval = {fifelse(round(go_Term_pval,3)==0,"<0.001",as.character(round(go_Term_pval,4)))}'), 
  #          x =100,  y = 4, angle = 0)+
  # ylim(-5,5)+
  ggtitle(glue::glue('{Protein_OI} interactors {term} fdr {max_fdr}'))

# imputing the missing values of each protein to -4
ProCan_matrix_norm_imp = copy(ProCan_matrix_norm_centered)
ProCan_matrix_norm_imp[,1:20] |> boxplot()
ProCan_matrix_norm_imp |> dim()
n_missing = ProCan_matrix_norm_imp |> is.na() |> matrixStats::colSums2()
n_present =  nrow(ProCan_matrix_norm_imp)-n_missing
ProCan_matrix_norm_imp = ProCan_matrix_norm_imp[,n_present>49]
Proteins_present = ProCan_matrix_norm_centered |> as.data.frame() |> 
  tibble::rownames_to_column('Project_Identifier') |> 
  as.data.table() |> 
  melt(id.vars = 'Project_Identifier', value.name = 'Abundance', variable.name = 'Uniprot')
Proteins_present[,presence := fifelse(is.na(Abundance),'Imputed','Measured')]
ProCan_matrix_norm_imp[is.na(ProCan_matrix_norm_imp)] <- rnorm(sum(is.na(ProCan_matrix_norm_imp)),
                                                               mean = -4,
                                                               sd = 0.3)
ProCan_matrix_norm_imp[,1:20] |> boxplot()

complex_id = 'huMAP3_00356.1'
subunit_complex = data.table()
for(complex_id in clusters_to_check[, clustID]){
  print(complex_id)
  partners_tmp = humap_3[clustID ==complex_id, uniprotACCs] |> str_split(' ') |> unlist()
  complex_procan = ProCan_matrix_norm_imp[,colnames(ProCan_matrix_norm_imp) %in% partners_tmp]
  if(is.matrix( complex_procan)){
    if(ncol(complex_procan)>2){
      for(Protein_tmp in colnames(complex_procan)){
        print(Protein_tmp)
        position_protein = which(colnames(complex_procan)==Protein_tmp)
        min_cell_lines = complex_procan[,position_protein] |> sort() |> head(50) 
        max_cell_lines = complex_procan[,position_protein] |> sort() |> tail(50) 
        # complex_procan_plot = complex_procan |> t()
        # rownames(complex_procan_plot) = rownames(complex_procan_plot) |> convert_to_gene()
        # complex_procan_plot = complex_procan_plot |> as.data.frame() |> tibble::rownames_to_column('ID') |> as.data.table() |>
        #   melt(id.vars ='ID',variable.name ='Cell_line',value.name = 'Abundance')
        # complex_procan_plot[,Expression:= fcase(
        #   Cell_line %in% names(min_cell_lines),'PSMB5_low',
        # Cell_line %in% names(max_cell_lines),'PSMB5_high',
        # !(Cell_line %in% names(c(max_cell_lines,min_cell_lines))),'PSBM5_average')]
        # 
        # order_samples = c(paste0('PSMA',1:7),
        #                   paste0('PSMB',1:10),
        #                   paste0('PSME',1:4),
        #                   paste0('PSMG',1:4))
        # # order_samples = order_samples |> sort(decreasing = F)
        # order_samples = c('PSMB5',order_samples, colnames(complex_procan) |> convert_to_gene()) |> unique()
        # complex_procan_plot |> ggplot(aes(x =factor(ID,levels = rev(order_samples)), y = Abundance,colour =Expression ))+
        #   geom_jitter(data = complex_procan_plot[Expression =='PSBM5_average'],alpha = 0.1) +
        #   geom_jitter(data = complex_procan_plot[Expression !='PSBM5_average'],alpha = 0.8) +
        #   coord_flip()+theme_bw()+
        #   labs(y='Imputted Cell Line Abundance',x = 'Protein ID')+
        #   scale_colour_manual(values = c('PSMB5_low'='darkblue',
        #                                  PSMB5_high = 'darkred',
        #                                  PSMB5_average = 'grey50'))
        diff_prot = mean(max_cell_lines)- mean(min_cell_lines)
        
        subunit_complex_tmp = tibble::enframe((complex_procan[max_cell_lines  |> names(),-position_protein] |> colMeans())-
                                                (complex_procan[min_cell_lines |> names(),-position_protein] |> colMeans() )) |>
          as.data.table()
        subunit_complex_tmp[,`:=`(subunit = Protein_tmp)] 
        subunit_complex = rbind(subunit_complex,subunit_complex_tmp)
      }
    }
  }
}
setnames(subunit_complex,c('name','subunit','value'),c('Protein_1','Protein_2','diff_abundance'))
subunit_complex[Protein_1<Protein_2,`:=`(Protein_1=Protein_2,
                                        Protein_2 = Protein_1)]
# removing same pairs that exist in multiple complexes
subunit_complex = subunit_complex |> unique()
subunit_complex[Protein_1 =='Q8WVM7'& Protein_2 == 'Q8N3U4']
subunits  = c('Q8NFD5','O14497')
subunits  = c('P28074','P28062')

ProCan_tissues = ProCan_matrix_norm_imp[,subunits] |> as.data.frame() |> tibble::rownames_to_column('Project_Identifier') |> 
  as.data.table()
ProCan_tissues = merge(ProCan_tissues |> melt(id.vars = 'Project_Identifier',
                                              variable.name = 'Uniprot',
                                              value.name = 'Abundance'),mapping_1[,.(Cancer_type,Project_Identifier,Cell_line)])
order_types = ProCan_tissues[Uniprot == 'P28062',.(ARID1A_abundance= median(Abundance)),by = Cancer_type
][order(ARID1A_abundance),Cancer_type]
ProCan_tissues[,Cancer_type:=factor(Cancer_type,levels = order_types)]
ProCan_tissues[,Gene := fifelse(Uniprot  =='P28074', 'PSMB5','PSMB8')]
ProCan_tissues |> ggplot(aes(x = Cancer_type, y  = Abundance,colour = Gene))+
  geom_boxplot()+
  coord_flip()+theme_bw()


fwrite(subunit_complex, here::here('out','datasets','pairs_different_presence_imp.gz'))
subunit_complex= fread(here::here('out','datasets','pairs_different_presence_imp.gz'))
# subunit_complex = fread(here::here('out','datasets','pairs_different_presence_imp.gz'))
subunit_complex_pairs = merge(subunit_complex,pairs_to_clust_mutexc,by =c('Protein_1','Protein_2'), all.x = T)
subunit_complex_pairs[,Type:= fifelse(is.na(cmmn_p_homodimer),as.character(!is.na(clustID)),cmmn_p_homodimer)]

subunit_complex_pairs[,Type:= fcase(
  Type =='no','interphase overlap',
  Type =='yes','interphase overlap (homodimer)',
  Type =='TRUE','subunit pairs',
  Type =='FALSE','non-subunits pairs'
) ]
ggplot(subunit_complex_pairs,aes(x= diff_abundance, fill =Type))+
  geom_density(alpha = 0.5)+theme_bw()+
  scale_fill_manual(values =subunit_colour)+
  labs(x = 'Expression Similarity', y = 'Density')+
  ggtitle('Protein Pairs similar expression covariation')

ggsave(here::here('out','plots','Procan_expression_similarity.pdf'))

### humap3 jaccard score #### 
humap_3_pairwise =  fread(here::here('in','datasets','ComplexPortal_model_unlabeledpredictions_ALL_sorted_202305.pairsWprob'))
setnames(humap_3_pairwise,c('Protein_1','Protein_2','PPI_score'))
# threshold for positive interactors
interactors = humap_3_pairwise[PPI_score>0.8]
uniprots_interactors = humap_3_pairwise[,.(Protein_1,Protein_2)] |> unlist() |> unique()
setkeyv(interactors,c('Protein_1','Protein_2'))
interactor_ids =list()
for(i in uniprots_interactors){
  print(which(uniprots_interactors ==i))
  interactors_tmp = interactors[Protein_1 ==i | Protein_2 ==i
  ][,.(Protein_1,Protein_2)]|> unlist() |> unique()
  interactor_ids[[i]] = interactors_tmp[interactors_tmp != i]
  
}
n_interactors = interactor_ids |> lengths()

possible_complex = interactor_ids[n_interactors>3]
jaccard_prots = calc_pairwise_overlaps(possible_complex)
jaccard_prots = as.data.table(jaccard_prots)
setnames(jaccard_prots,c('name1','name2'),c('Protein_1','Protein_2'))
jaccard_prots[Protein_1<Protein_2,`:=`(Protein_1=Protein_2,
                                         Protein_2 = Protein_1)]
jaccard_prots |> nrow()
jaccard_prots |> unique() |> nrow()
jaccard_prots[Protein_1<Protein_2]
fwrite(jaccard_prots,here::here('out','datasets', 'jaccard_mutexcl.gz'))
jaccard_prots = fread(here::here('out','datasets', 'jaccard_mutexcl.gz'))
merge(jaccard_prots,pairs_to_clust,by =c('Protein_1','Protein_2'), all.x = T) |> 
  ggplot(aes(x= jaccard, fill = is.na(clustID)))+
  geom_density(alpha = 0.5)+
  ggtitle('Protein Pairs in the same complex have higher covariation')
merge(jaccard_prots,mut_excl,by =c('Protein_1','Protein_2'), all.x = F) |> 
  ggplot(aes(x= jaccard, fill = cmmn_p_homodimer))+
  geom_density(alpha = 0.5)+
  ggtitle('Protein mut_excl Pairs with w/o homodimer')
ggsave(here::here('out','plots','Procan_mutexc_jacc.pdf'))
### zscores and combining scores #### 
pairs_to_clust[Protein_1<Protein_2]
subunit_complex[Protein_1<Protein_2]
Pairwise_procan[Protein_1<Protein_2]
# testing_procan = copy(Pairwise_procan)
# testing_procan[ ,global_bicor := zscoring(bicor_Procan_across_tissue)]
jaccard_prots[Protein_1<Protein_2]
pair_intertissue_complex_cor[Protein_1<Protein_2]

    # pairs might be missing in global bicor because of not enough shared data points
merge_with_complex = merge(pairs_to_clust,
                           pair_intertissue_complex_cor[, .(Protein_1,Protein_2,tissue_cor,clustID) ],
                           all.x = T, by = c('Protein_1','Protein_2','clustID'))
combined_scores = Reduce(function(...) merge(..., all.x=T, by = c('Protein_1','Protein_2')), list(merge_with_complex,
                                                                                                  Tissue_level_covar,
                                                                                                  Pairwise_procan
                                                                                                  # ,
                                                                                                  # jaccard_prots[,.(Protein_1,Protein_2,jaccard)]
                                                                                                  ))
combined_scores |> nrow()
combined_scores  =  merge(combined_scores,
                          subunit_complex,
                       all = T, allow.cartesian = T, by = c('Protein_1','Protein_2'))  
combined_scores |> nrow()


combined_scores[,`:=`(complex_tissue_bicor=  zscoring(tissue_cor),
                      tissue_bicor =  zscoring(tissue_level_covar),
                      global_bicor = zscoring(bicor_Procan_across_tissue),
                      diff_abundance = zscoring(diff_abundance)
                      # ,
                      # pairwise_jacc = -(jaccard)
                      )]
combined_scores = combined_scores[,.(Protein_1,Protein_2,clustID,
                             complex_tissue_bicor,global_bicor,diff_abundance
                             # ,pairwise_jacc
                             )] |>
  melt(id.vars = c('Protein_1','Protein_2','clustID'), value.name = 'score',variable.name = 'type')
combined_scores = combined_scores|> unique()
combined_scores_summary = combined_scores[!is.na(score)]
combined_scores_summary = combined_scores_summary[,.(avg_score = mean(score),
                                            N_types = .N), by = .(Protein_1, Protein_2, clustID)]


# some entries like this only have bicor because they are complexes with less than 4 subunits (found in procan)
combined_scores[Protein_1=='O94966' & Protein_2=='O00418']
combined_scores_summary = combined_scores_summary[N_types>1]
combined_scores_summary = merge(combined_scores_summary,Gene_names,
                        by.x = c('Protein_1'),by.y = c('Protein'))
combined_scores_summary = merge(combined_scores_summary,Gene_names,
                        by.x = c('Protein_2'),by.y = c('Protein'))


# merge_b |> unique() |> View()
combined_scores_summary[Protein_1 %in% c('P28062','P28074') &
                  Protein_2 %in% c('P28062','P28074')]

combined_scores_summary  =merge(combined_scores_summary, max_subunit_cor,all.x = T, by.x = c('Protein_1','clustID'),
                        by.y = c('subunit','complex'),allow.cartesian = T)
combined_scores_summary  =merge(combined_scores_summary, max_subunit_cor,all.x = T, by.x = c('Protein_2','clustID'),
                        by.y = c('subunit','complex'),allow.cartesian = T)
combined_scores_summary[,min_max_cor:= min(c(complex_cor.y,complex_cor.x),na.rm =T ), by = .(Protein_1,Protein_2,clustID)]
combined_scores_summary = merge(combined_scores_summary,mut_excl,all.x = T, by = c('Protein_1','Protein_2')) 

# check the average correlation of the specific complex in the specific tissue
high_cor_complexes = tissue_complex_cor[,.(avg_complex_cor = mean(complex_cor)), by = .(complex,tissue)] 
high_cor_complexes = high_cor_complexes[avg_complex_cor>0.4,complex]
combined_scores_summary[,complex_cor:= fifelse(clustID %in% high_cor_complexes,'high_cor_complex','low_cor_complex')]


# Gene_names_complexes = merge_b[order(min_max_cor,bicor_Procan_across_tissue)][,head(.SD,1),by = .(ID.x,ID.y,Protein_2,Protein_1)]
plotting_mutexcl = merge(combined_scores_summary,jaccard_prots[,.(Protein_1,Protein_2,jaccard)],by = c('Protein_1','Protein_2'))
plotting_mutexcl |> 
  ggplot(aes(x = jaccard, fill = cmmn_p_homodimer, alpha= 0.5))+
  geom_density()

# mut_excl_old = fread(here::here('in','datasets','refiltered_mutually_exclusive_and_structurally_consistent_protein_pairs_03OCT2024.tsv') )
# mut_excl_old =mut_excl_old[,.(protein_1,protein_2,cmmn_p_homodimer)]
# mut_excl = fread(here::here('in','datasets','MutEx_pairs_refiltered_w_large_comps_removed_04OCT2024.tsv'))
# # '/mnt/kustatscher/members/savvas/LFQ_DIA_SWATH/mutually_exclusive_dimer_analysis_10SEP2024.csv')
# # mut_excl = fread(here::here('in','datasets','humap3_mutual_exclusion_vs_non_dimers_annotated_homodimers_30JUL2024.csv'))
# mut_excl = merge(mut_excl,mut_excl_old,c('protein_1','protein_2'))
# mut_excl  =mut_excl[,.(V1,V2,common_protein,protein_1,protein_2,interface_overlap,cmmn_p_homodimer)]
# mut_excl[,Class:='mut_excl']
# setnames(mut_excl,c('protein_1','protein_2'),c('Protein_1','Protein_2'))
# mut_excl[Protein_1<Protein_2,`:=`(Protein_1=Protein_2,
#                                   Protein_2 = Protein_1)]
# mut_excl[Protein_1<Protein_2]

struct_diff_plot = merge(combined_scores_summary[,.(Protein_1,Protein_2,avg_score)],
      mut_excl)


t_test = t.test(struct_diff_plot[class == 'mutually_exclusive' ,avg_score], 
                struct_diff_plot[class != 'mutually_exclusive',avg_score], 
                alternative = 'two.sided', var.equal = FALSE)
  ggplot(struct_diff_plot,aes(x = class,y = avg_score))+
  geom_boxplot()+theme_bw()+
 ggtitle('Compaaring ALL-IN score including new pairs revisions',subtitle = glue::glue('t.test pvalue<0.0001'))+
    # labs(x= 'Class', 'ALL-IN Score')
    # facet_wrap('cmmn_p_homodimer')+
  labs(y = 'Mutually exclusivity ProCAN score (ALL-IN)',x = 'Class')
ggsave(here::here('out','plots','structurally_consistent_vsmut_Excl.pdf'))

# loaidng ProHD2 RF scores and removing isoform info
ProHD2_treeclust_RF <- fread(here::here('out','datasets','RF_predictions.csv.gz'))
ProHD2_treeclust_RF = ProHD2_treeclust_RF[,.(Protein_1,Protein_2,RF_covariation_prob)]
isoform_pattern <- "^(.+)(;\\1-\\d+)+$"  # This regex patterns identify protein groups that contain different isoforms of the same protein, e.g."A6QL63;A6QL63-2;A6QL63-3;A6QL63-4;A6QL63-5"
ProHD2_treeclust_RF[ Protein_1 %like% isoform_pattern, Protein_1 := gsub(";.+", "", Protein_1) ]   # Where all proteins are isoforms, simplify to gene level
ProHD2_treeclust_RF[ Protein_2 %like% isoform_pattern, Protein_2 := gsub(";.+", "", Protein_2) ]   # Where all proteins are isoforms, simplify to gene level

merge(ProHD2_treeclust_RF,mut_excl, by = c('Protein_1','Protein_2')) |> 
  ggplot(aes(x = interface_overlap,y = RF_covariation_prob))+
  geom_boxplot()+theme_bw()+
  facet_wrap('cmmn_p_homodimer')+
  labs(y = 'ProHD2 RF score',x = 'Interphase Overlap')
ggsave(here::here('out','plots','structurally_consistent_vsmut_Excl_ProHD2.pdf'))

modelled =  merge(combined_scores_summary[,.(Protein_1,Protein_2,avg_score,ID.x,ID.y)],
                  mut_excl, all.x = T, c('Protein_1','Protein_2'))

paralogs = merge(modelled,jaccard_prots, by = c('Protein_1','Protein_2'), all.x = T)
paralogs= merge(paralogs,paralog_pairs,by = c('Protein_1','Protein_2'))     |> as.data.table()        
paralogs[,has_structure := !is.na(interface_overlap)]                 
paralogs = paralogs[,.(Protein_1,Protein_2,ID.x,ID.y,interface_overlap,cmmn_p_homodimer,
                       has_structure,jaccard,avg_score )]
fwrite(paralogs,here::here('out','datasets','paralogs_no-homo_ALL-IN.csv'))
paralogs= paralogs[order(-avg_score), head(.SD,min(.N,1)), by = .(Protein_1,Protein_2)]


ggplot(paralogs,
       aes(x= jaccard,y = avg_score, label =paste(ID.x,ID.y, sep = '-'), 
           colour = interface_overlap, alpha= !is.na(interface_overlap)  ))+
  geom_point()+
  theme_bw()+
  facet_wrap('has_structure')+
  scale_colour_manual(values = c('darkred','darkblue'))+
  scale_alpha_manual(values = rev(c(0.5,1)))+
  labs(x = 'Jaccard Interactome score', y = 'Pairwise subunit mut_excl score', 
       colour = 'Interface Overlap', alpha = 'Interface Overlap')+
  ggtitle('Mutual exclusivity score from ProCAN expression data')+
  # ggrepel::geom_label_repel(data = Gene_names_complexes[ID.y %in% c('PSMB5','PSMB8','ARPC1B','ARPC1A') &
  #                                                           ID.x %in% c('PSMB5','PSMB8','ARPC1B','ARPC1A')  ])+
  ggrepel::geom_text_repel(max.overlaps = 10)
ggsave(here::here('out','plots','mut_excpairs_pairs_paralogs_sep_structure.pdf'),width = 16,height = 12)


pairs = mut_excl_plot[][avg_score<(-1.4),.(Protein_1,Protein_2,ID.x,ID.y)] |> unique()
pairs= merge(pairs,pairs_to_clust, by = c('Protein_1','Protein_2'))
pairs = pairs[,head(.SD,1),by = .(Protein_1,Protein_2)]
pairs = paralogs[is.na(cmmn_p_homodimer) & avg_score<(-0.8),head(.SD,1),by = .(Protein_1,Protein_2) ][,.(Protein_1,Protein_2,ID.x,ID.y,clustID)]
to_check = c('ILF2','ZFR','STRBP','ILF3')
to_check_uniprots =Gene_names[ID %in% to_check]
prots_per_complex[uniprotACCs %in% to_check_uniprots$Protein][,.(.N), by = clustID]

combined_scores_summary[ID.x %in% to_check& ID.y %in%to_check ] |> View()
complex_tmp = 'huMAP3_07115.1'
combined_scores[!is.na(score)][Protein_1 %in% to_check_uniprots$Protein &
                                 Protein_2 %in% to_check_uniprots$Protein] |> View()

merge(Pairwise_procan,pairs_to_clust[clustID == complex_tmp,.(Protein_1,Protein_2)],
      by = c('Protein_1','Protein_2')) |> View()

for(i in 1:nrow(pairs)){
  subunits_tmp = merge(Pairwise_procan,pairs_to_clust[clustID == pairs[i,clustID],.(Protein_1,Protein_2)],
                       by = c('Protein_1','Protein_2'))
  subunits_tmp = rbind(subunits_tmp,
                       data.table(Protein_1 = unique(subunits_tmp[,.(Protein_1,Protein_2)] |> unlist()),
                                  Protein_2 = unique(subunits_tmp[,.(Protein_1,Protein_2)] |> unlist()),
                                  bicor_Procan_across_tissue = 1))
  subunits_tmp = subunits_tmp |>   dcast(Protein_1~Protein_2,value.var = 'bicor_Procan_across_tissue') |> 
    tibble::column_to_rownames('Protein_1') |> as.matrix()
  subunits_tmp = Matrix::forceSymmetric(subunits_tmp,uplo="L")
  subunits_tmp = subunits_tmp |> as.matrix()
  rownames(subunits_tmp) = rownames(subunits_tmp) |> convert_to_gene()
  colnames(subunits_tmp) = colnames(subunits_tmp) |> convert_to_gene()
  subunits_tmp[is.na(subunits_tmp)] <- 0
  pdf(here::here('out','plots',glue::glue('{pairs[i,ID.x]} {pairs[i,ID.y]} {pairs[i,clustID]} heatmap.pdf')))
  
  subunits_tmp |> pheatmap::pheatmap(main = glue::glue('{pairs[i,ID.x]} {pairs[i,ID.y]} {pairs[i,clustID]}'))
  dev.off()
  
  # subunit_complex = subunit_complex |> unique()
  # subunit_complex[Protein_1 =='P19022'& Protein_2 == 'P12830']
  subunits  = pairs[i,.(Protein_1,Protein_2)] |> unlist()
  
  ProCan_tissues = ProCan_matrix_norm_imp[,subunits] |> as.data.frame() |> tibble::rownames_to_column('Project_Identifier') |> 
    as.data.table()
  ProCan_tissues = merge(ProCan_tissues |> melt(id.vars = 'Project_Identifier',
                                                variable.name = 'Uniprot',
                                                value.name = 'Abundance'),mapping_1[,.(Cancer_type,Project_Identifier,Cell_line)])
  order_types = ProCan_tissues[Uniprot == subunits[1],.(abundance= median(Abundance)),by = Cancer_type
  ][order(abundance),Cancer_type]
  ProCan_tissues[,Cancer_type:=factor(Cancer_type,levels = order_types)]
  ProCan_tissues[,Gene := fifelse(Uniprot  ==subunits[1], pairs[i,ID.x], pairs[i,ID.y])]
  ProCan_tissues[,N_cell_lines := .N, by = .(Gene,Cancer_type )]
  ProCan_tissues[,Min_cell_lines := min(N_cell_lines), by = .(Cancer_type )]
  ProCan_tissues = ProCan_tissues[Min_cell_lines>10]
  ProCan_tissues[,Cancer_type:=factor(Cancer_type,levels = order_types)]
  tissues_to_plot = rbind(ProCan_tissues[,.(tissue_avg = mean(Abundance)), by  = .(Uniprot,Cancer_type)][,head(.SD,3),by= Uniprot],
        ProCan_tissues[,.(tissue_avg = mean(Abundance)), by  = .(Uniprot,Cancer_type)][,tail(.SD,3),by= Uniprot])
  
  ProCan_tissues = merge(ProCan_tissues,Proteins_present[,.(Uniprot,Project_Identifier,presence)],by = c('Uniprot','Project_Identifier'))
  ProCan_tissues = ProCan_tissues[Cancer_type %in% tissues_to_plot$Cancer_type] 
  N_detected = ProCan_tissues[,.(N_detection = .N), by = .(Gene,Cancer_type,presence,N_cell_lines)] 
  N_detected[,perc_N := (N_detection/N_cell_lines)*100]
  N_detected = rbind(N_detected,
                     data.table(Gene = 'SEPTIN10',
                                Cancer_type = c("Burkitt's Lymphoma",'Acute Myeloid Leukemia'),
                                presence = 'Measured',
                                N_cell_lines = c(12,28),
                                N_detection = 0,
                                perc_N = 0))
  order_tissue = N_detected[Gene =='SEPTIN6' & presence=='Measured',][order(perc_N),Cancer_type] |> as.character()
  N_detected[presence=='Measured'] |> ggplot(aes(x = factor(Cancer_type,levels = order_tissue), fill = Gene,y = perc_N  ))+
    geom_col(position = 'dodge')+ theme_bw()+
    scale_x_discrete(guide = guide_axis(n.dodge=2))+
    labs(x = 'Cancer Lineages',y = 'Percentage of cell lines detected')+
    scale_fill_manual(values = c('darkgreen','darkblue'))+
    geom_text(aes(y = 105, label = N_cell_lines))
  ggsave(here::here('out','plots',glue::glue('{pairs[i,ID.x]}_{pairs[i,ID.y]} missingness.pdf')))
  
  
  ProCan_tissues[presence =='Measured'] |> 
    ggplot(aes(x =  Cancer_type, y  = Abundance,colour =  Gene ))+
    geom_boxplot(outliers = F)+
    geom_point(position = position_dodge(width = 0.8))+
    coord_flip()+
    theme_bw()+
    labs(x = 'Relative Abundance (z-score)', y = 'Cancer Lineage')+
    theme(legend.position = 'bottom')+
    scale_colour_manual(values = c('darkgreen','darkblue'))+
    # scale_fill_manual(values = rev(c('white','grey50')))+
    scale_alpha_manual(values = c(0.4,0.8))
  
  ggsave(here::here('out','plots',glue::glue('{pairs[i,ID.x]}_{pairs[i,ID.y]} expression.pdf')),width = 4.4,height = 3.3)
  
  ProCan_tissues|> 
    ggplot(aes(x =  Cancer_type, y  = Gene,  alpha= presence ))+
    # geom_boxplot(outliers = F)+
   geom_jitter(data = ProCan_tissues[presence =='Measured'], height = 0.2, width = 0.1, 
               aes( colour =  Gene))+
    geom_jitter(data = ProCan_tissues[presence !='Detected'], shape = 1,  
                height = 0.2, width = 0.1, colour = 'grey50'
    )+
    # coord_flip()+
    theme_bw()+
    labs(y = 'Detected Protein', x = 'Cancer Lineage')+
    theme(legend.position = 'bottom',
          axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    # scale_fill_manual(values = c('darkgreen','darkblue'))+
    scale_alpha_manual(values = c(0.6,0.95))+
    scale_x_discrete(guide = guide_axis(n.dodge=2))+
    
    scale_colour_manual(values = c('darkgreen','darkblue'))
  
  presence_heatmap = ProCan_tissues[,.(Project_Identifier,Cancer_type,Abundance, presence,Gene)]
  cell_line_order = presence_heatmap[Gene=='SEPTIN10'][order(Cancer_type,Abundance),Project_Identifier]
  presence_heatmap[presence =='Imputed',Abundance:=NaN]
  
  annotation_row = data.frame(row.names = presence_heatmap[Gene=='SEPTIN10', Project_Identifier],
                               lineage = presence_heatmap[Gene=='SEPTIN10',Cancer_type]) 

  presence_heatmap = presence_heatmap |> dcast(Project_Identifier~Gene,value.var = 'Abundance') |> 
    tibble::column_to_rownames('Project_Identifier') |> as.matrix()
  presence_heatmap = presence_heatmap[cell_line_order,]
  presence_heatmap |> pheatmap::pheatmap(cluster_rows = F, annotation_row = annotation_row,
                                         cluster_cols = F,show_rownames = F)  
}



pairs_confidence = merge(clusters_to_check[,.(clustID,cluster_confidence)],pairs_to_clust,by = 'clustID')
new_pairs = merge(plotting_mutexcl[is.na(cmmn_p_homodimer)][jaccard>0.1][avg_score<0],pairs_confidence, by = c('Protein_1','Protein_2'))
fwrite(new_pairs[cluster_confidence==1][complex_cor == 'high_cor_complex'], here::here('out','datasets','mut_exl_pairs_no_structure.csv'))
plotting_mutexcl = plotting_mutexcl[!is.na(cmmn_p_homodimer)
                                    ][,head(.SD,1),by = .(Protein_1,Protein_2,jaccard,avg_score,cmmn_p_homodimer,ID.x,ID.y)]
  ggplot(plotting_mutexcl, aes(x= jaccard,y = avg_score,colour =cmmn_p_homodimer, label =paste(ID.x,ID.y, sep = '-')  ))+
  geom_point()+
    theme_bw()+
    facet_wrap('cmmn_p_homodimer')+
    scale_colour_manual(values = c('darkred','darkblue'))+
  labs(x = 'Jaccard Interactome score', y = 'Pairwise subunit mut_excl score')+
  ggtitle('Mut_exc score is an average of various metrics',
          subtitle = 'metrics: global covariation, lineage-covariation, presence/absence, tissue complex correlation')+
  # ggrepel::geom_label_repel(data = Gene_names_complexes[ID.y %in% c('PSMB5','PSMB8','ARPC1B','ARPC1A') &
  #                                                           ID.x %in% c('PSMB5','PSMB8','ARPC1B','ARPC1A')  ])+
  ggrepel::geom_text_repel(data = plotting_mutexcl[!between(avg_score,-1,1)
                                                   ][str_detect(ID.y ,           paste(c('PSMB','WIPF','SEPTIN','ARPC1','COPS7','AP1S','WASF','STAG','SMARCD','PSMA'),collapse = '|')) &
                                                                              str_detect(ID.x, paste(c('PSMB','WIPF','SEPTIN','ARPC1','COPS7','AP1S','WASF','STAG','SMARCD','PSMA'),collapse = '|'))])
  ggsave(here::here('out','plots','mut_excpairs_pairs_residues.pdf'),width = 11,height = 9)
  fwrite(plotting_mutexcl[!is.na(cmmn_p_homodimer),.(Protein_1,Protein_2,ID.x,ID.y,common_protein,clustID, avg_score,jaccard)
                          ][,head(.SD,1),by = .(Protein_1,Protein_2,avg_score)] |> unique()
         ,here::here('out','datasets','mut_excpairs_pairs_residues.gz'))
  
to_plot = combined_scores_summary[order(Class,avg_score,complex_cor),.(ID.x,ID.y,avg_score,complex_cor,min_max_cor,Class,cmmn_p_homodimer)][,head(.SD,1),by =.(ID.x,ID.y)]
# to_plot = 
# to_plot
unique_proteins_to_plot= c(combined_scores_summary$Protein_1,combined_scores_summary$Protein_2) |> unique()
to_plot = data.table()
for(i in unique_proteins_to_plot){
  print(i)
  to_plot = rbind(to_plot,
                  combined_scores_summary[order(cmmn_p_homodimer,avg_score,complex_cor),
                          .(Protein_1,Protein_2,ID.x,ID.y,avg_score,complex_cor,min_max_cor,Class,cmmn_p_homodimer)
                          ][Protein_1 == i | Protein_2 == i][,head(.SD,1)])
}
to_plot = to_plot[order(avg_score)][,head(.SD,1),by =.(ID.x,ID.y)]
to_plot = merge(to_plot,jaccard_prots[,.(Protein_1,Protein_2,jaccard)],by = c('Protein_1','Protein_2'))
# to_plot = combined_scores_summary[order(cmmn_p_homodimer,avg_score,complex_cor),
#                                   .(ID.x,ID.y,avg_score,complex_cor,min_max_cor,Class,cmmn_p_homodimer)
#                                   ][,head(.SD,1),by =.(ID.y)][,head(.SD,1),by =.(ID.x)]
to_plot|> ggplot(aes(x = jaccard, y = avg_score,colour =complex_cor  , label = paste(`ID.x`,`ID.y` ,sep = '-')))+
  geom_point(colour = 'grey40', alpha = 0.3)+theme_bw()+
  geom_point(data =  to_plot[!is.na(Class)], alpha = 0.7)+
  labs(x = 'JAcc score', y = 'Pairwise subunit mut_excl score')+
  ggtitle('Mut_exc score is an average of various metrics',
          subtitle = 'metrics: global pairwise, humap3 jaccard score, presence/absence, tissue complex correlation')+
  # ggrepel::geom_label_repel(data = Gene_names_complexes[ID.y %in% c('PSMB5','PSMB8','ARPC1B','ARPC1A') &
  #                                                           ID.x %in% c('PSMB5','PSMB8','ARPC1B','ARPC1A')  ])+
  ggrepel::geom_text_repel(data = to_plot[avg_score<(-1) & min_max_cor>0.2 & Class == 'mut_excl' ], size = 3,  max.overlaps = 100)+
  ggrepel::geom_text_repel(data = to_plot[avg_score<(-0.5) & min_max_cor>0.6  & is.na(Class)], size = 3,  max.overlaps = 10)+
  facet_wrap('cmmn_p_homodimer')+
  ggrepel::geom_text_repel(data = to_plot[avg_score>1.5 & min_max_cor>0.8& Class == 'mut_excl'], size = 3,  max.overlaps = 10)+
  
  ggrepel::geom_text_repel(data = to_plot[avg_score<(-0.5) & min_max_cor<0.2& Class == 'mut_excl'], size = 3,  max.overlaps = 10)+
  scale_colour_manual(values = c('darkred','darkblue','grey40'))

to_plot|> ggplot(aes(x = jaccard, y = avg_score, label = paste(`ID.x`,`ID.y` ,sep = '-')))+
  geom_point(colour = 'grey40', alpha = 0.1)+theme_bw()+
  geom_point(data =  to_plot[!is.na(Class)],  alpha = 0.7)+
  labs(x = 'Min covariation of protein pair to rest of complex', y = 'Pairwise subunit mut_excl score')+
  ggtitle('Mut_exc score is an average of various metrics',
          subtitle = 'metrics: global pairwise,  presence/absence, tissue complex correlation')+
  # ggrepel::geom_label_repel(data = Gene_names_complexes[ID.y %in% c('PSMB5','PSMB8','ARPC1B','ARPC1A') &
  #                                                           ID.x %in% c('PSMB5','PSMB8','ARPC1B','ARPC1A')  ])+
  ggrepel::geom_text_repel(data = to_plot[avg_score<(-1) & min_max_cor>0.2 & Class == 'mut_excl' ], size = 3,  max.overlaps = 100)+
  ggrepel::geom_text_repel(data = to_plot[avg_score<(-0.5) & min_max_cor>0.6  & is.na(Class)], size = 3,  max.overlaps = 10)+
  # facet_wrap('cmmn_p_homodimer')+
  ggrepel::geom_text_repel(data = to_plot[avg_score>1.5 & min_max_cor>0.8& Class == 'mut_excl'], size = 3,  max.overlaps = 10)+
  
  ggrepel::geom_text_repel(data = to_plot[avg_score<(-0.5) & min_max_cor<0.2& Class == 'mut_excl'], size = 3,  max.overlaps = 10)+
  scale_colour_manual(values = c('darkred','darkblue','grey40'))

mut_excl_plot = to_plot[!is.na(Class)][cmmn_p_homodimer =='no']
mut_excl_plot|> ggplot(aes(x = jaccard, y = avg_score, label = paste(`ID.x`,`ID.y` ,sep = '-')))+
  # geom_point(colour = 'grey40', alpha = 0.3)+theme_bw()+
  geom_point(colour = 'darkred', alpha = 0.7)+
  # facet_wrap('cmmn_p_homodimer')+
  theme_bw()+
  labs(x = 'Interactome  Jaccard scores', y = 'Subunit pair Mutual Exclusivity score')+
  ggtitle('Mutual exclusivity score from ProCAN expression data')+
  ggrepel::geom_text_repel(data = mut_excl_plot[str_detect(ID.y , 
                                                      paste(c('PSMB','SEPTIN','ARPC1','COPS7','AP1S','WASF','STAG','SMARCD','PSMA'),collapse = '|')) &
                                             str_detect(ID.x, paste(c('PSMB','SEPTIN','ARPC1','COPS7','AP1S','WASF','STAG','SMARCD','PSMA'),collapse = '|'))])+
  # ggrepel::geom_text_repel(data = to_plot[avg_score<(-1) & min_max_cor>0.2 & Class == 'mut_excl' ], size = 3,  max.overlaps = 100)+
  # ggrepel::geom_text_repel(data = to_plot[avg_score<(-0.5) & min_max_cor>0.6  & is.na(Class)], size = 3,  max.overlaps = 10)+
  # # facet_wrap('cmmn_p_homodimer')+
  # ggrepel::geom_text_repel(data = to_plot[avg_score>1.5 & min_max_cor>0.8& Class == 'mut_excl'], size = 3,  max.overlaps = 10)+
  # 
  # ggrepel::geom_text_repel(data = to_plot[avg_score<(-0.5) & min_max_cor<0.2& Class == 'mut_excl'], size = 3,  max.overlaps = 10)+
  # scale_colour_manual(values = c('darkred','darkblue','grey40'))+
  geom_hline(aes(yintercept = 0), linetype = 'dashed')
ggsave(here::here('out','plots','mut-excl_score_nohomodimers.pdf'))
fwrite(mut_excl_plot[,.(Protein_1,Protein_2,ID.x,ID.y,avg_score,jaccard)],
       here::here('out','datasets','ALL_in_mut_Excl_score.csv'))

combined_scores_summary[avg_score<(-0.2) & min_max_cor>0.4] |> View()
combined_scores_summary[avg_score<(-0.2)  & is.na(Class)] |> View()
top_pairs = combined_scores_summary[avg_score<(-0.1)  ][order(avg_score),.(ID.x,ID.y,avg_score,min_max_cor,Class,cmmn_p_homodimer)][,head(.SD,1),by =.(ID.x,ID.y)] 
  ggplot(top_pairs,aes(x = min_max_cor, y = avg_score,colour = cmmn_p_homodimer, label = paste(`ID.x`,`ID.y` ,sep = '-')))+
  geom_point(colour = 'grey40', alpha = 0.3)+theme_bw()+
  geom_point(data =  top_pairs[!is.na(Class)], alpha = 0.7)+
  labs(x = 'Min covariation of protein pair to rest of complex', y = 'Pairwise subunit mut_excl score')+
  ggtitle('Mut_exc score is an average of various metrics',
          subtitle = 'metrics: global pairwise, humap3 jaccard score, presence/absence, tissue complex correlation')+
  # ggrepel::geom_label_repel(data = Gene_names_complexes[ID.y %in% c('PSMB5','PSMB8','ARPC1B','ARPC1A') &
  #                                                           ID.x %in% c('PSMB5','PSMB8','ARPC1B','ARPC1A')  ])+
  ggrepel::geom_text_repel(data = top_pairs[avg_score<(-0.5) & min_max_cor>0.4 ], max.overlaps = 100)+
  ggrepel::geom_text_repel(data = top_pairs[avg_score>1.5 & min_max_cor>0.8& Class == 'mut_excl'], max.overlaps = 10)+
  ggrepel::geom_text_repel(data = top_pairs[avg_score<(-0.5) & min_max_cor<0.4& Class == 'mut_excl'], max.overlaps = 10)+
  scale_colour_manual(values = c('darkred','darkblue','grey40'))
  ggsave(here::here('out','plots','mut_exlusive_top_per_prot.pdf'),width=24, height = 9)
SMARCDs = c('Q96GM5','Q92925','Q6STE5')
STKS = c('Q9Y6E0','Q9P289')
ARIDs = c('Q8NFD5' ,'O14497')
example_proteins = ARIDs
combined_scores_summary[Protein_1 %in% example_proteins & Protein_2 %in% example_proteins] |> View()
examples = combined_scores[Protein_1 %in% example_proteins & Protein_2 %in% example_proteins] 
examples = merge(examples,Gene_names,
                 by.x = c('Protein_1'),by.y = c('Protein'))
examples = merge(examples,Gene_names,
                 by.x = c('Protein_2'),by.y = c('Protein'))

ggplot(examples,aes(x = type, y = score, colour = paste(ID.x,ID.y,sep = '-')))+
  geom_jitter(height = 0.0, width = 0.02)+
  facet_wrap('clustID')


### plot #### 
### complexes which are downregulated in specific tissues #### 
library('EnrichIntersect')
# cell_lines_tmp = mapping_1[Cancer_type == i,Project_Identifier]
ProCan_matrix_norm_complex = copy(ProCan_matrix_norm)
prots_NA= colMeans(is.na(ProCan_matrix_norm_complex)) 
hist(prots_NA)
ProCan_matrix_norm_complex = ProCan_matrix_norm_complex[,prots_NA<0.8]
ProCan_matrix_norm_complex |> dim()
# centering each protein
tmp_medians <- apply( ProCan_matrix_norm_complex , 2, median, na.rm = TRUE )  
ProCan_matrix_norm_complex <- sweep( ProCan_matrix_norm_complex , 2,tmp_medians, FUN = "-" )
ProCan_matrix_norm_complex[,1:20] |> boxplot()
ProCan_matrix_norm_complex |> hist(breaks = 100)
# imputing the missing values of each protein to -4
ProCan_matrix_norm_complex[is.na(ProCan_matrix_norm_complex)] <- rnorm(sum(is.na(ProCan_matrix_norm_complex)),
                                                               mean = -2,
                                                               sd = 0.3)

tmp_medians <- apply( ProCan_matrix_norm_complex , 2, median, na.rm = TRUE )  
ProCan_matrix_norm_complex <- sweep( ProCan_matrix_norm_complex , 2,tmp_medians, FUN = "-" )



complex_procan = ProCan_matrix_norm_complex |> as.data.frame() |> 
  tibble::rownames_to_column('Project_Identifier') |>
  as.data.table() |> 
  melt(id.vars = 'Project_Identifier', value.name = 'zscore_abundance',variable.name = 'Protein')
complex_procan = merge(complex_procan,  prots_per_complex[cluster_confidence==1
                                                          ][,cluster_confidence:=NULL], 
                       by.y  = 'uniprotACCs', by.x = 'Protein')
complex_procan_cell_lines = merge(complex_procan[clustID %in% complexes_with_enough_prots[N>3,clustID]],
                                  mapping_1[Cancer_type %in% tissue_n[N>5,Cancer_type],.(Project_Identifier,Cancer_type)], by = 'Project_Identifier')




complex_procan = complex_procan_cell_lines[,.(SD_complex = sd(zscore_abundance),
                                              mean = mean(zscore_abundance),
                                              N_subunits = .N),by = .(Cancer_type,clustID)] 

complex_procan_enrichment = complex_procan_cell_lines[,.(SD_complex = sd(zscore_abundance),
                                                         mean = mean(zscore_abundance),
                                                         N_subunits = .N),by = .(Cancer_type,Project_Identifier , clustID)] 

complex_procan_enrichment  = complex_procan_enrichment |>
  dcast(Project_Identifier ~ clustID ,value.var = 'mean') |> 
  tibble::column_to_rownames('Project_Identifier') |> 
  as.matrix()


custom.set = data.table(drug = mapping_1$Project_Identifier,
                        group = mapping_1$Cancer_type)
custom.set = custom.set[drug %in% rownames(complex_procan_enrichment)]

set.seed(123)
enrich <- enrichment(complex_procan_enrichment, custom.set, permute.n = 100)
values_complexes = enrich$pvalue |> as.data.frame() |> 
  tibble::rownames_to_column('Cancer_type') |> 
  as.data.table() |> 
  melt(id.vars = 'Cancer_type',value.name = 'pvalue', variable.name = 'clustID')
values_complexes[,padj:= p.adjust(pvalue,method = 'fdr')]
NES_complexes = enrich$S |> as.data.frame() |> 
  tibble::rownames_to_column('Cancer_type') |> 
  as.data.table() |> 
  melt(id.vars = 'Cancer_type',value.name = 'NES', variable.name = 'clustID')
enrich_complexes = merge(NES_complexes,values_complexes, by = c('Cancer_type','clustID'))
enrich_complexes[,pvalue_plot := pvalue+rnorm(1,0.03,0.0065), by = .(clustID,Cancer_type) ]
order_types = enrich_complexes[clustID %in% c('huMAP3_00154.1')][order(-NES),Cancer_type]
enrich_complexes[,Cancer_type := factor(Cancer_type,levels =order_types)]
enrich_complexes |> 
  ggplot(aes(x = NES, y = -log10(pvalue_plot)))+
  geom_point(aes(alpha = padj<0.01))+
  theme_bw()+
  geom_point(data = enrich_complexes[clustID %in% c('huMAP3_00115.1',
                                                    'huMAP3_00154.1',
                                                    'huMAP3_02748.1')],
             aes(colour = clustID))+
  facet_wrap('Cancer_type')+
  theme(legend.position = 'bottom')+
  labs(x = 'Normalised Enrichment Score', y = '-log10(pvalue)')
ggsave(here::here('out','plots','NES_clustID_cancertypes.pdf'),width = 25,height = 15)

  
diff_complexes =   complex_procan[mean>1 &SD_complex<1.2][,.N,by = clustID][order(N)] |> tail(3)
diff_complexes = diff_complexes[,clustID]

original_data = ProCan_matrix_norm_centered |> as.data.frame() |> 
  tibble::rownames_to_column('Project_Identifier') |>
  as.data.table() |> 
  melt(id.vars = 'Project_Identifier', value.name = 'imputation',variable.name = 'Protein')
plot_diff_complexes =  ProCan_matrix_norm_complex |> as.data.frame() |> 
  tibble::rownames_to_column('Project_Identifier') |>
  as.data.table() |> 
  melt(id.vars = 'Project_Identifier', value.name = 'abundance',variable.name = 'Protein')
plot_diff_complexes = merge(plot_diff_complexes,original_data, by = c('Project_Identifier','Protein') )
plot_diff_complexes[,imputation:= fifelse(is.na(imputation),'imputted','not_imputted')]
plot_diff_complexes = merge(plot_diff_complexes,  prots_per_complex[cluster_confidence==1
][,cluster_confidence:=NULL], by.y  = 'uniprotACCs', by.x = 'Protein')
plot_diff_complexes = merge(plot_diff_complexes[clustID %in% diff_complexes],
                       mapping_1[Cancer_type %in% tissue_n[N>5,Cancer_type],
                                 .(Project_Identifier,Cancer_type)], by = 'Project_Identifier')
plot_diff_complexes = merge(plot_diff_complexes, Gene_names, by = 'Protein')

order_samples = complex_procan[clustID == tail(diff_complexes,1)][order(mean),Cancer_type]

plot_diff_complexes[,Cancer_type:=factor(Cancer_type,levels = order_samples)]

plot_diff_complexes |> 
  ggplot(aes(y = Cancer_type, x = abundance ))+
  geom_boxplot()+
  geom_jitter( alpha = 0.1,aes(colour = imputation))+
  facet_wrap('clustID')+theme_bw()+
  labs(x = 'Protein Relative Abundance', colour = 'Imputation', y = 'Cancer Lineage')+
  ggrepel::geom_text_repel(data = plot_diff_complexes[order(abundance),head(.SD,1),by = .(ID,clustID)],
                           max.overlaps = 40,aes(y = Cancer_type, x = abundance,label = ID))

ggsave(here::here('out','plots','diff_complex_abundance.pdf'),width=9, height = 9)
fwrite(plot_diff_complexes,here::here('out','plots','diff_complex_abundance.gz'))
### complexes which are poorly covaring in tissues #### 

tissue_cor_complex_all = merge(tissue_cor,pairs_to_clust, by = c('Protein_1','Protein_2') ) 
tissue_cor_complex_all[,N_pairs:= .N,  by = .(tissue,clustID) ]
tissue_cor_complex_all = tissue_cor_complex_all[N_pairs>10]
tissue_cor_complex = tissue_cor_complex_all[,.(avg_complex = mean(bicor),
                      sd_complex= sd(bicor),
                      N_pairs = .N), by = .(tissue,clustID)]
tissue_cor_complex[,tissue_per_complex := .N, by = clustID]
tissue_cor_complex = tissue_cor_complex[tissue_per_complex>12]
tissue_cor_complex_plot= merge(tissue_cor_complex,tissue_n,by.x ='tissue',by.y= 'Cancer_type') 
  ggplot(tissue_cor_complex_plot,aes(y= tissue, x = avg_complex ))+
  geom_boxplot()+
  geom_point(data = tissue_cor_complex_plot[clustID=='huMAP3_00356.1'],alpha = 0.7,colour = '#FF0000')+
    # geom_point(data = tissue_cor_complex_plot[clustID=='huMAP3_04814.1'],alpha = 0.7,colour = 'green')+
    theme_bw()+
    ggtitle('Proteasome is well correlated in some tissues but poorly in others')
  ggsave(here::here('out','plots','diff_complex_covariation_proteasome.pdf'),width=9, height = 9)
fwrite(tissue_cor_complex_plot,here::here('out','plots','diff_complex_covariation_proteasome.gz'))
tissue_cor_complex[,z_score_complex := zscoring(avg_complex), by = clustID]

tissue_cor_complex |> ggplot(aes(x = z_score_complex, y = sd_complex, label  = clustID))+
  geom_point()+
  facet_wrap('tissue')+
  ggrepel::geom_label_repel(data = tissue_cor_complex[z_score_complex>2 & sd_complex<0.2])

clust_tmp = 'huMAP3_00356.1'
prots_per_complex[clustID==clust_tmp,uniprotACCs] |> convert_to_gene()

complex_test = ProCan_matrix_norm_complex |> as.data.frame() |> 
  tibble::rownames_to_column('Project_Identifier') |>
  as.data.table() |> 
  melt(id.vars = 'Project_Identifier', value.name = 'zscore_abundance',variable.name = 'Protein')
complex_test = merge(complex_test,  prots_per_complex[clustID ==clust_tmp][,cluster_confidence:=NULL], 
by.y  = 'uniprotACCs', by.x = 'Protein')
complex_test = merge(complex_test,
                       mapping_1[Cancer_type %in% tissue_n[N>5,Cancer_type],
                                 .(Project_Identifier,Cancer_type)], by = 'Project_Identifier')
merge(complex_test[clustID==clust_tmp],tissue_cor_complex[clustID ==clust_tmp],by.x = c('clustID','Cancer_type'),
      by.y =c('clustID','tissue') ,allow.cartesian = T) |> 
  ggplot(aes(y=Cancer_type, x = zscore_abundance,fill  =avg_complex        ))+
  geom_boxplot()+
  ggtitle(glue::glue('protein expression {clust_tmp} vs covariation (fill)'),
          subtitle = paste(prots_per_complex[clustID==clust_tmp,uniprotACCs] |> convert_to_gene() |> 
                             head(4),collapse = ' '))+
  theme_bw()+
  scale_fill_continuous(type = "viridis")
ggsave(here::here('out','plots','diff_complex_covariation_expression_proteasome.pdf'),width=9, height = 9)
fwrite(tissue_cor_complex_plot,here::here('out','plots','diff_complex_covariation_expression_proteasome.gz'))

tissue_cor_complex_all[clustID %in% clust_tmp] |> 
  ggplot(aes(x= bicor,y = tissue))+
  ggridges::geom_density_ridges()+theme_bw()+
  ggtitle(glue::glue('{clust_tmp}  covariation'),
          subtitle = paste(prots_per_complex[clustID==clust_tmp,uniprotACCs] |> convert_to_gene() |> 
                             head(4),collapse = ' '))
# gastric and neuroblastoma;Mesothelioma
clust_tmp = 'huMAP3_03221.1'
tissue_tmp = c("Ewing's Sarcoma")
neuro_heatmap = tissue_cor_complex_all[clustID %in% clust_tmp & tissue ==tissue_tmp,.(Protein_1,Protein_2,bicor)] 
neuro_heatmap = rbind(neuro_heatmap,
                      data.table(Protein_1 = unique(neuro_heatmap[,.(Protein_1,Protein_2)] |> unlist()),
                                 Protein_2 = unique(neuro_heatmap[,.(Protein_1,Protein_2)] |> unlist()),
                      bicor = 1))
neuro_heatmap = neuro_heatmap |>   dcast(Protein_1~Protein_2,value.var = 'bicor') |> 
  tibble::column_to_rownames('Protein_1') |> as.matrix()
neuro_heatmap = Matrix::forceSymmetric(neuro_heatmap,uplo="L")
neuro_heatmap = neuro_heatmap |> as.matrix()
rownames(neuro_heatmap) = rownames(neuro_heatmap) |> convert_to_gene()
neuro_heatmap |> pheatmap::pheatmap(main = tissue_tmp)
#### conditional subunits  ####
conditional_subunits = tissue_complex_cor[,N_tissue_subunit := .N,  by = .(complex,subunit)] 
conditional_subunits = conditional_subunits[N_tissue_subunit>7]
# average correlation of the subunit to the comple across tissues
conditional_subunits[, avg_cor_sub := mean(complex_cor), by = .(complex,subunit,ID)]
# average correlation of the all subunits to the complex across tissue
conditional_subunits[, avg_cor_complex := mean(complex_cor), by = .(complex)]
# average correlation of the all subunits to the complex per tissue
conditional_subunits[, avg_cor_complex_tissue := mean(complex_cor), by = .(complex,tissue)]
# zscoring to find tissues where complexes loses corr
conditional_subunits[, z_avg_cor_complex_tissue := zscoring(avg_cor_complex_tissue), by = .(complex)]
# how different is the corr of this subunit (to the complex) in the specific tissue, compared to the average of tissues
conditional_subunits[, corr_diff := complex_cor-avg_cor_sub, by = .(complex,tissue,subunit,ID)]

conditional_subunits[, subunit_to_complex_dif := avg_cor_sub-avg_cor_complex, by = .(complex,subunit,ID)]
conditional_subunits

conditional_subunits[ avg_cor_sub >0.4 & avg_cor_complex>0.4 ] |> ggplot(aes(x =  subunit_to_complex_dif, y = corr_diff, label  = ID, colour = tissue))+
  geom_point()+
  # facet_wrap('tissue')+
  labs(x = 'correlation difference of subunit to complex compared to the complex average subunits', y ='correlation difference in specific tissue vs across tissues')+
  ggrepel::geom_label_repel(data = conditional_subunits[abs(corr_diff)>0.3 & abs(subunit_to_complex_dif)<0.2 & avg_cor_sub>0.5  & avg_cor_complex>0.4  ])

complex_id = 'huMAP3_01276.1'

tissue_complex_cor[complex == complex_id] |> 
  ggplot(aes(x = ID, y = complex_cor, colour = tissue))+
  geom_point()

conditional_subunits[ avg_cor_sub >0.4 & avg_cor_complex_tissue>0.4 & abs(subunit_to_complex_dif)<0.2 ] |> 
  ggplot(aes(x =  z_avg_cor_complex_tissue, y = corr_diff, label  = ID, colour = tissue))+
  geom_point()+
  # facet_wrap('tissue')+
  labs(x = 'how well the complex covaries in the specific tissue', y ='correlation difference in specific tissue vs across tissues')+
  ggrepel::geom_label_repel(data = conditional_subunits[corr_diff<(-0.2) & z_avg_cor_complex_tissue>0.5  &avg_cor_sub>0.4  & avg_cor_complex_tissue>0.4  & abs(subunit_to_complex_dif)<0.2 ])
complex_id = 'huMAP3_08556.1'
tissue_tmp = 'Glioblastoma'
ID_tmp = 'MRPS12'
tissue_complex_cor[complex == complex_id] |> 
  ggplot(aes(x = ID, y = complex_cor, colour = tissue))+
  geom_point(data = tissue_complex_cor[complex == complex_id][tissue == tissue_tmp], size = 2)+
  geom_point(data = tissue_complex_cor[complex == complex_id][tissue!= tissue_tmp], alpha = 0.2)+
  coord_flip()+
  labs(x = 'Subunits Correlation to the complex',y = 'Subunit')
ggsave(here::here('out','plots','conditional_subunit_per_tissue.pdf'),width=9, height = 9)

original_data = ProCan_matrix_norm_centered |> as.data.frame() |> 
  tibble::rownames_to_column('Project_Identifier') |>
  as.data.table() |> 
  melt(id.vars = 'Project_Identifier', value.name = 'imputation',variable.name = 'Protein')
differential_presence =  ProCan_matrix_norm_complex |> as.data.frame() |> 
  tibble::rownames_to_column('Project_Identifier') |>
  as.data.table() |> 
  melt(id.vars = 'Project_Identifier', value.name = 'abundance',variable.name = 'Protein')
differential_presence = merge(differential_presence,original_data, by = c('Project_Identifier','Protein') )
differential_presence[,imputation:= fifelse(is.na(imputation),'imputted','not_imputted')]
differential_presence = merge(differential_presence,  prots_per_complex[cluster_confidence<5
][,cluster_confidence:=NULL], by.y  = 'uniprotACCs', by.x = 'Protein',allow.cartesian = T)
differential_presence = merge(differential_presence,
                            mapping_1[Cancer_type %in% tissue_n[N>5,Cancer_type],
                                      .(Project_Identifier,Cancer_type)], by = 'Project_Identifier')
to_plot = differential_presence[clustID == complex_id][Cancer_type %in% tissue_complex_cor[complex == complex_id,tissue]]

tissue_order = to_plot[Protein %in% to_plot[,Protein][1]][,.(avg_abundance = mean(abundance)), by= Cancer_type][order(avg_abundance),Cancer_type]
to_plot = merge(to_plot,Gene_names,by = 'Protein',all.x=T)

  ggplot(to_plot,aes(y = factor(Cancer_type,levels = tissue_order), x = abundance))+
  geom_boxplot()+
    facet_wrap('ID')+
  geom_jitter(aes(colour =  imputation), alpha =0.1)

  tissue_order = to_plot[Protein %in% 'O15235'][,.(avg_abundance = mean(abundance)), by= Cancer_type][order(avg_abundance),Cancer_type]
  to_plot = merge(to_plot,Gene_names,by = 'Protein',all.x=T)
  to_plot = merge(to_plot, tissue_complex_cor[complex == complex_id & ID == ID_tmp,.(ID,tissue,complex_core_zscore_pair)],
                  by.x = c('ID','Cancer_type'), by.y = c('ID','tissue'))
  ggplot(to_plot[Protein =='O15235'],aes(y = factor(Cancer_type,levels = tissue_order),
                                         fill = complex_core_zscore_pair, x = abundance))+
    geom_boxplot()+
    # facet_wrap('ID')+
    geom_jitter(aes(colour =  imputation), alpha =0.5)+
    ggtitle(glue::glue('{ID_tmp} expression against covariation to complex'))
  
  
  tissue_order = to_plot[ID %in%ID_tmp][,.(avg_abundance = mean(abundance)), by= Cancer_type][order(avg_abundance),Cancer_type]
  # to_plot = merge(to_plot,Gene_names,by = 'Protein',all.x=T)
  ggplot(to_plot[ID %in%ID_tmp],aes(y = factor(Cancer_type,levels = tissue_order), x = abundance))+
    geom_boxplot()+
    # facet_wrap('ID')+
    geom_jitter(aes(colour =  imputation), alpha =0.1)
# ITG
  proteins_int = humap_3[clustID == 'huMAP3_06815.1',uniprotACCs ] |> 
    str_split(' ') |> unlist()
  Pairwise_procan[Protein_1 %in% proteins_int & Protein_2 %in% proteins_int ]
  subunit_complex[Protein_1 %in% proteins_int & Protein_2 %in% proteins_int ]
  ProCan_prots = ProCan |> colnames() |> str_remove(';[:print:]*$')
  ProCan_prots[ProCan_prots %in% proteins_int]
  