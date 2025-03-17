# PROCAN - Hu.MAP3
library(data.table);library(stringr); library(ggplot2);library(WGCNA) 
#### functions used #### 
# zscoring function to make scores comparable
zscoring = function(x){
  (x-mean(x,na.rm =T ))/sd(x,na.rm=T)
}

# bicor correlation
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

# Uniprots into Gene names
convert_to_gene <- function(x){
  gene_names = c()
  for(i in x){
    Uniprots = strsplit(i,';') |> unlist() |> str_remove('-.$')
    Genes = Uniprot_mapping[Type == 'Gene_Name' & Protein %in% Uniprots, ID] |> unique() |> paste(collapse = ';')
    gene_names[i] = Genes
  }
  gene_names
}
# jaccard overlap function
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
# removing redundancies
Gene_names = Gene_names[order(ID,N_ids)][, tail(.SD, 1), by = ID]
Gene_names = Gene_names[,.(ID,Protein)][, tail(.SD, 1), by = Protein]

# reading paralogs pairs
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
# ordering and removing redundancies
paralog_pairs[Protein_1<Protein_2,`:=`(Protein_1=Protein_2,
                                       Protein_2 = Protein_1)]
paralog_pairs[Protein_1<Protein_2]

### reading humap3 clusters & define clusters to use #### 
humap_3 = fread(here::here('in','datasets','Supplemental_Table_3_hu.MAP3.1_complexes_wConfidenceScores_total15326_20240922.csv'),header = T,  sep = ',')
humap_3[,N_members := stringr::str_count(uniprotACCs,' ')+1, by = clustID]
clusters_to_check = humap_3[cluster_confidence<4 & N_members>3]

# making pairwise table
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
# ordering and removing redundacies
pairs_to_clust[Protein_1<Protein_2,`:=`(Protein_1=Protein_2,
                                        Protein_2 = Protein_1)]
pairs_to_clust[Protein_1<Protein_2]

fwrite(pairs_to_clust,here::here('out','datasets','pairs_to_clust.gz'))
pairs_to_clust = fread(here::here('out','datasets','pairs_to_clust.gz'))

### reading mutually exclusive sets from interface #### 
new_pairs = fread(here::here('in','datasets',
                             'TableS4_structurally_consistent_mutually_exclusive_modeled_pairs_w_janes_w_burke_12MAR2024.csv'))
new_pairs = new_pairs[,.(pair_1,pair_2,interface,common_protein,protein_1,protein_2)]
setnames(new_pairs,c('pair_1','pair_2','interface'),c('V1','V2','class'))

mut_excl = new_pairs
setnames(mut_excl,c('protein_1','protein_2'),c('Protein_1','Protein_2'))
mut_excl[Protein_1<Protein_2,`:=`(Protein_1=Protein_2,
                                  Protein_2 = Protein_1)]
mut_excl[Protein_1<Protein_2]

### reading Procan and normalise #### 
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

### Global bicor #### 
Pairwise_procan = f_BIC(ProCan_matrix_norm,'bicor_Procan_across_tissue',30)
# fwrite(Pairwise_procan,here::here('out','datasets','ProCAN_global_bicor_min_30.gz'))
Pairwise_procan = fread(here::here('out','datasets','ProCAN_global_bicor_min_30.gz'))
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
    cor_tmp =  f_BIC(ProCan_matrix_complex_tissue,'bicor',20)
    cor_tmp[,tissue:=i]
    tissue_cor = rbind(tissue_cor[!is.na(bicor)],
                       cor_tmp)
  }
}
fwrite(tissue_cor,here::here('out','datasets','procan_covariation_per_tissue.gz'))
tissue_cor = fread(here::here('out','datasets','procan_covariation_per_tissue.gz'))

tissue_cor[,N_tissues := .N, by = c('Protein_1','Protein_2')]

### tissue specific complex bicor #### 
tissue_complex_cor = data.table()
proteins_correlations =  unique(c(tissue_cor$Protein_1,tissue_cor$Protein_2))

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
fwrite(tissue_complex_cor,'per_tissue_subunit_cor.gz')
tissue_complex_cor = fread('per_tissue_subunit_cor.gz')
tissue_complex_cor = tissue_complex_cor[!is.na(complex_core_cor)]
tissue_complex_cor[,N_tissues := .N, by = c('subunit','complex')]
max_subunit_cor = tissue_complex_cor[order(complex_cor,subunit,complex)][,tail(.SD,1) ,by = .(subunit,complex)]
tissue_complex_cor = merge(tissue_complex_cor,Gene_names,by.x = 'subunit',by.y = 'Protein',all.x=T)

### by-tissue complex bicor #### 
pair_intertissue_complex_cor = data.table()
for(complex_tmp in clusters_to_check[, clustID]){
  print(complex_tmp)
  subunit_complex_tmp = tissue_complex_cor[complex ==complex_tmp ]
  if(length(unique(subunit_complex_tmp$tissue))>4){
    complex_cor_tmp = dcast(subunit_complex_tmp,tissue ~subunit,value.var = 'complex_core_zscore_pair' ) |> tibble::column_to_rownames('tissue') |> 
      f_BIC('tissue_cor',8)
    complex_cor_tmp[,clustID:=complex_tmp]
    pair_intertissue_complex_cor = rbind(pair_intertissue_complex_cor,
                                         complex_cor_tmp[!is.na(tissue_cor)])
    
  }
}
fwrite(pair_intertissue_complex_cor,here::here('out','datasets','mutexl_pairs_intertissue_complex_cor.csv'))
pair_intertissue_complex_cor = fread(here::here('out','datasets','mutexl_pairs_intertissue_complex_cor.csv'))

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

### Procan Imputation and score #### 
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

# complex_id = 'huMAP3_00356.1'
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
fwrite(subunit_complex, here::here('out','datasets','pairs_different_presence_imp.gz'))
subunit_complex= fread(here::here('out','datasets','pairs_different_presence_imp.gz'))

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
# combined_scores[Protein_1=='O94966' & Protein_2=='O00418']
combined_scores_summary = combined_scores_summary[N_types>1]
combined_scores_summary = merge(combined_scores_summary,Gene_names,
                                by.x = c('Protein_1'),by.y = c('Protein'))
combined_scores_summary = merge(combined_scores_summary,Gene_names,
                                by.x = c('Protein_2'),by.y = c('Protein'))
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
plotting_mutexcl = merge(combined_scores_summary,
                         jaccard_prots[,.(Protein_1,Protein_2,jaccard)],by = c('Protein_1','Protein_2'))
plotting_mutexcl |> 
  ggplot(aes(x = jaccard, fill = class, alpha= 0.5))+
  geom_density()

struct_diff_plot = merge(combined_scores_summary[,.(Protein_1,Protein_2,avg_score)],
                         mut_excl)

t_test = t.test(struct_diff_plot[class == 'mutually_exclusive' ,avg_score], 
                struct_diff_plot[class != 'mutually_exclusive',avg_score], 
                alternative = 'two.sided', var.equal = FALSE)
ggplot(struct_diff_plot,aes(x = class,y = -avg_score))+
  geom_boxplot()+theme_bw()+
  ggtitle('Comparing ALL-IN score including new pairs revisions',subtitle = glue::glue('t.test pvalue<0.0001'))+
  # labs(x= 'class', 'ALL-IN Score')
  # facet_wrap('cmmn_p_homodimer')+
  labs(y = 'Mutually exclusivity ProCAN score (ALL-IN)',x = 'class')
ggsave(here::here('out','plots','structurally_consistent_vsmut_Excl_revision.pdf'))

# modelled =  merge(combined_scores_summary[,.(Protein_1,Protein_2,avg_score,ID.x,ID.y)],
#                   mut_excl, all.x = T, c('Protein_1','Protein_2'))
# 
# paralogs = merge(modelled,jaccard_prots, by = c('Protein_1','Protein_2'), all.x = T)
# paralogs= merge(paralogs,paralog_pairs,by = c('Protein_1','Protein_2'))     |> as.data.table()        
# paralogs[,classified := !is.na(class)]                 
# paralogs = paralogs[,.(Protein_1,Protein_2,ID.x,ID.y,class,
#                        classified,jaccard,avg_score )]
# fwrite(paralogs,here::here('out','datasets','paralogs_no-homo_ALL-IN.csv'))
# paralogs= paralogs[order(-avg_score), head(.SD,min(.N,1)), by = .(Protein_1,Protein_2)]


# ggplot(paralogs,
#        aes(x= jaccard,y = avg_score, label =paste(ID.x,ID.y, sep = '-'), 
#            colour = class, alpha= !is.na(class)  ))+
#   geom_point()+
#   theme_bw()+
#   # facet_wrap('has_structure')+
#   scale_colour_manual(values = c('darkred','darkblue'))+
#   scale_alpha_manual(values = rev(c(0.5,1)))+
#   labs(x = 'Jaccard Interactome score', y = 'Pairwise subunit mut_excl score', 
#        colour = 'Interface Overlap', alpha = 'Interface Overlap')+
#   ggtitle('Mutual exclusivity score from ProCAN expression data')+
#   ggrepel::geom_text_repel(max.overlaps = 10)
# ggsave(here::here('out','plots','mut_excpairs_pairs_paralogs_sep_structure.pdf'),width = 16,height = 12)


to_plot = combined_scores_summary[order(class,avg_score,complex_cor),
                                  .(ID.x,ID.y,avg_score,complex_cor,min_max_cor,class)
                                  ][,head(.SD,1),by =.(ID.x,ID.y)]

# selecting to show most extreme mututally exclusive candidate candidate pair protein
unique_proteins_to_plot= c(combined_scores_summary$Protein_1,combined_scores_summary$Protein_2) |> unique()
to_plot = data.table()
for(i in unique_proteins_to_plot){
  print(i)
  to_plot = rbind(to_plot,
                  combined_scores_summary[order(avg_score,complex_cor),
                                          .(Protein_1,Protein_2,ID.x,ID.y,avg_score,complex_cor,min_max_cor,class)
                  ][Protein_1 == i | Protein_2 == i][,head(.SD,1)])
}
to_plot = to_plot[order(avg_score)][,head(.SD,1),by =.(ID.x,ID.y)]
to_plot = merge(to_plot,jaccard_prots[,.(Protein_1,Protein_2,jaccard)],by = c('Protein_1','Protein_2'))
mut_excl_plot = to_plot[!is.na(class)]
mut_excl_plot|> ggplot(aes(x = jaccard, y = -avg_score, label = paste(`ID.x`,`ID.y` ,sep = '-')))+
  # geom_point(colour = 'grey40', alpha = 0.3)+theme_bw()+
  geom_point(colour = 'darkred', alpha = 0.7)+
  # facet_wrap('cmmn_p_homodimer')+
  theme_bw()+
  labs(x = 'Interactome  Jaccard scores', y = 'Subunit pair Mutual Exclusivity score')+
  ggtitle('Mutual exclusivity score from ProCAN expression data (ALL-IN)')+
  ggrepel::geom_text_repel(data = mut_excl_plot[str_detect(ID.y , 
                                                           paste(c('PSMB','SEPTIN','ARPC1','COPS7','AP1S','WASF','STAG','SMARCD','PSMA'),collapse = '|')) &
                                                  str_detect(ID.x, paste(c('PSMB','SEPTIN','ARPC1','COPS7','AP1S','WASF','STAG','SMARCD','PSMA'),collapse = '|'))])+
  geom_hline(aes(yintercept = 0), linetype = 'dashed')
ggsave(here::here('out','plots','mut-excl_score_nohomodimers_revision.pdf'))
all_mut_exclusive_subunits = merge(combined_scores_summary[,.(Protein_1,Protein_2,ID.x,ID.y,avg_score,class)],
      jaccard_prots[,.(Protein_1,Protein_2,jaccard)],by = c('Protein_1','Protein_2'),all.x = T)

fwrite(all_mut_exclusive_subunits[,.(Protein_1,Protein_2,ID.x,ID.y,avg_score,jaccard)],
       here::here('out','datasets','ALL_in_score_all_pairs.csv'))

  # example PSMB5 - PSMB8
  subunits = c('P28074','P28062')
  ProCan_tissues = ProCan_matrix_norm_imp[,subunits] |> as.data.frame() |> tibble::rownames_to_column('Project_Identifier') |> 
    as.data.table()
  ProCan_tissues = merge(ProCan_tissues |> melt(id.vars = 'Project_Identifier',
                                                variable.name = 'Uniprot',
                                                value.name = 'Abundance'),
                         mapping_1[,.(Cancer_type,Project_Identifier,Cell_line)])
  order_types = ProCan_tissues[Uniprot == subunits[1],.(abundance= median(Abundance)),by = Cancer_type
  ][order(abundance),Cancer_type]
  ProCan_tissues[,Cancer_type:=factor(Cancer_type,levels = order_types)]
  ProCan_tissues[,Gene := fifelse(Uniprot  ==subunits[1], 'PSMB5', 'PSMB8')]
  ProCan_tissues[,N_cell_lines := .N, by = .(Gene,Cancer_type )]
  ProCan_tissues[,Min_cell_lines := min(N_cell_lines), by = .(Cancer_type )]
  ProCan_tissues = ProCan_tissues[Min_cell_lines>10]
  ProCan_tissues[,Cancer_type:=factor(Cancer_type,levels = order_types)]
  tissues_to_plot = rbind(ProCan_tissues[,.(tissue_avg = mean(Abundance)), by  = .(Uniprot,Cancer_type)][,head(.SD,3),by= Uniprot],
                          ProCan_tissues[,.(tissue_avg = mean(Abundance)), by  = .(Uniprot,Cancer_type)][,tail(.SD,3),by= Uniprot])
  
  ProCan_tissues = merge(ProCan_tissues,Proteins_present[,.(Uniprot,Project_Identifier,presence)],by = c('Uniprot','Project_Identifier'))
  ProCan_tissues = ProCan_tissues[Cancer_type %in% tissues_to_plot$Cancer_type] 
  # protein relative expression
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
  
  # example SEPTIN6-SEPTIN10
  subunits = c('Q14141','Q9P0V9')
  ProCan_tissues = ProCan_matrix_norm_imp[,subunits] |> as.data.frame() |> tibble::rownames_to_column('Project_Identifier') |> 
    as.data.table()
  ProCan_tissues = merge(ProCan_tissues |> melt(id.vars = 'Project_Identifier',
                                                variable.name = 'Uniprot',
                                                value.name = 'Abundance'),
                         mapping_1[,.(Cancer_type,Project_Identifier,Cell_line)])
  order_types = ProCan_tissues[Uniprot == subunits[1],.(abundance= median(Abundance)),by = Cancer_type
  ][order(abundance),Cancer_type]
  ProCan_tissues[,Cancer_type:=factor(Cancer_type,levels = order_types)]
  ProCan_tissues[,Gene := fifelse(Uniprot  ==subunits[1], 'SEPTIN6', 'SEPTIN10')]
  ProCan_tissues[,N_cell_lines := .N, by = .(Gene,Cancer_type )]
  ProCan_tissues[,Min_cell_lines := min(N_cell_lines), by = .(Cancer_type )]
  ProCan_tissues = ProCan_tissues[Min_cell_lines>10]
  ProCan_tissues[,Cancer_type:=factor(Cancer_type,levels = order_types)]
  tissues_to_plot = rbind(ProCan_tissues[,.(tissue_avg = mean(Abundance)), by  = .(Uniprot,Cancer_type)][,head(.SD,3),by= Uniprot],
                          ProCan_tissues[,.(tissue_avg = mean(Abundance)), by  = .(Uniprot,Cancer_type)][,tail(.SD,3),by= Uniprot])
  
  ProCan_tissues = merge(ProCan_tissues,Proteins_present[,.(Uniprot,Project_Identifier,presence)],by = c('Uniprot','Project_Identifier'))
  ProCan_tissues = ProCan_tissues[Cancer_type %in% tissues_to_plot$Cancer_type] 
  # protein relative expression
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
  
