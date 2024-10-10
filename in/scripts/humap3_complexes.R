# humap3 adding covariation
library(data.table);library(stringr);library(ggplot2);library(clusterProfiler)
# Uniprot identifiers
# from https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/
Uniprot_mapping = fread(here::here('in','datasets','HUMAN_9606_idmapping.dat.gz'), header = F)
setnames(Uniprot_mapping,c('Protein','Type','ID'))
Gene_names = Uniprot_mapping[,N_ids := .N, by =Protein][Type == 'Gene_Name']
Gene_names = Gene_names[order(ID,N_ids)][, tail(.SD, 1), by = ID]
Gene_names = Gene_names[,.(ID,Protein)][, tail(.SD, 1), by = Protein]

# converting vectors of Uniprots to HUGO IDs
convert_to_gene <- function(x){
  gene_names = c()
  for(i in x){
    Uniprots = strsplit(i,';') |> unlist() |> str_remove('-.$')
    Genes = Uniprot_mapping[Type == 'Gene_Name' & Protein %in% Uniprots, ID] |> unique() |> paste(collapse = ';')
    gene_names[i] = Genes
  }
  gene_names
}
#loading humap pairwise scores
humap_3_pairwise =  fread(here::here('in','datasets','ComplexPortal_model_unlabeledpredictions_ALL_sorted_202305.pairsWprob'))
setnames(humap_3_pairwise,c('Protein_1','Protein_2','PPI_score'))
humap_3_pairwise_rev = copy(humap_3_pairwise)
humap_3_to_join = rbind(humap_3_pairwise,
                        humap_3_pairwise_rev[,`:=`(Protein_1 = Protein_2,
                                                   Protein_2 = Protein_1)]) |> unique()

# loaidng ProHD2 RF scores and removing isoform info
ProHD2_treeclust_RF <- fread(here::here('out','datasets','RF_predictions.csv.gz'))
ProHD2_treeclust_RF = ProHD2_treeclust_RF[,.(Protein_1,Protein_2,RF_covariation_prob)]
isoform_pattern <- "^(.+)(;\\1-\\d+)+$"  # This regex patterns identify protein groups that contain different isoforms of the same protein, e.g."A6QL63;A6QL63-2;A6QL63-3;A6QL63-4;A6QL63-5"
ProHD2_treeclust_RF[ Protein_1 %like% isoform_pattern, Protein_1 := gsub(";.+", "", Protein_1) ]   # Where all proteins are isoforms, simplify to gene level
ProHD2_treeclust_RF[ Protein_2 %like% isoform_pattern, Protein_2 := gsub(";.+", "", Protein_2) ]   # Where all proteins are isoforms, simplify to gene level

 

ProHD2_humap = merge(humap_3_to_join,ProHD2_treeclust_RF,all.x = F)
ProHD2_humap = merge(ProHD2_humap, Gene_names,all.x = T, by.x = 'Protein_1', by.y = 'Protein')
ProHD2_humap = merge(ProHD2_humap, Gene_names,all.x = T, by.x = 'Protein_2', by.y = 'Protein')
# ProHD2_humap[RF_covariation_prob>0.95 & PPI_score<0.1] |> View()
pairs_to_show = data.table(ID.x = c('COX6C','DHX9','RPL13','NDUFS5','TUBA1C','MRPS18B','POLR3GL','MTHFD2L','SLC12A7','MRPL30','POLR2H','UQCR10'),
                           ID.y = c('NDUFS4','DHX15','EIF3F','NDUFS8','TUBAL3','MRPS23','POLR1C','MTHFD2','SLC12A4','CDH3','IDH3A','POLD2'))
ProHD2_humap |> ggplot(aes(x = RF_covariation_prob, y = PPI_score, label = paste(ID.x,ID.y, sep = '-')))+
  geom_hex(binwidth = 0.01)+ theme_bw()+
  annotate("rect", xmin=c(0), xmax=c(0.15), ymin=c(0.85) , ymax=c(1), alpha=0.4, color="darkblue",linewidth = 1.1)+
  annotate("rect", xmin=c(0), xmax=c(0.15), ymin=c(0) , ymax=c(0.15), alpha=0.4, color="darkgreen",linewidth = 1.1)+
  annotate("rect", xmin=c(0.85), xmax=c(1), ymin=c(0.85) , ymax=c(1), alpha=0.4, color="darkred",linewidth = 1.1)+
  annotate("rect", xmin=c(0.85), xmax=c(1), ymin=c(0) , ymax=c(0.15), alpha=0.4, color="darkorange",linewidth = 1.1)+
  # scale_fill_continuous(type = "viridis") + 
  scale_fill_viridis_c(trans = "log10",
                       breaks = c(1e+00,1e+01,1e+02,1e+03,1e+04,1e+05))+
  # lims(x = c(0.1,0.99),y = c(0.0001,0.9999))+
  labs(x = 'ProHD2 covariation', y = 'huMAP3.0 PPI probability')
  # annotate("rect", xmin=c(0,0), xmax=c(3,5), ymin=c(20,10) , ymax=c(30,20), alpha=0.2, color="blue", fill="blue")+
  # ggrepel::geom_text_repel(data = merge(ProHD2_humap,pairs_to_show, by = c('ID.x','ID.y')),min.segment.length = 0, colour = 'grey30')
ggsave(here::here('out','plots','humap3_prohd2_hex_plot.pdf'),width = 9, height = 8)
ggsave(here::here('out','plots','humap3_prohd2_hex_plot.png'),width = 9, height = 8,device='png', dpi=900)

fwrite(ProHD2_humap,here::here('out','datasets','ProHD2_humap_pairwise_merge.gz'))
file.size(here::here('out','datasets','ProHD2_humap_pairwise_merge.gz'))/(1024 * 1024)

# categorising interactions
categories = list(low_lowhumap = unique(unlist(ProHD2_humap[RF_covariation_prob<0.15 & PPI_score<0.15,
                                                            .(Protein_1,Protein_2)])),
                  high_lowhumap = unique(unlist(ProHD2_humap[RF_covariation_prob>0.85 & PPI_score<0.15,
                                                             .(Protein_1,Protein_2)])),
                  low_highhumap = unique(unlist(ProHD2_humap[RF_covariation_prob<0.15 & PPI_score>0.85,
                                                             .(Protein_1,Protein_2)])),
                  high_highhumap = unique(unlist(ProHD2_humap[RF_covariation_prob>0.85 & PPI_score>0.85,
                                                              .(Protein_1,Protein_2)])))
values_col = c(lowProhd2_lowhumap ='darkgreen' ,
               highProhd2_lowhumap = 'darkorange',
               lowProhd2_highhumap = 'darkblue',
               highProhd2_highhumap = 'darkred',
               random ='grey50')
# check if these interactions are validated in open cell
OpenCell <- fread(here::here('in','datasets','opencell-protein-interactions.csv'))
OpenCell[target_gene_name <interactor_gene_name   ,`:=`(target_gene_name =interactor_gene_name ,
                                                        interactor_gene_name    = target_gene_name )]
OpenCell = OpenCell[str_detect(interactor_gene_name,';',negate = T) & 
           str_detect(target_gene_name,';',negate = T)]
OpenCell_coverage = unique(c(OpenCell$target_gene_name,OpenCell$interactor_gene_name))

# the four quadrants

ProHD2_humap[, class:= fcase(
  RF_covariation_prob<0.15 & PPI_score<0.15 , 'lowProhd2_lowhumap',
  RF_covariation_prob>0.85 & PPI_score<0.15,'highProhd2_lowhumap',
  RF_covariation_prob<0.15 & PPI_score>0.85, 'lowProhd2_highhumap',
  RF_covariation_prob>0.85 & PPI_score>0.85, 'highProhd2_highhumap',
  RF_covariation_prob<2,'other')]
reorder_ProHD2_humap = copy(ProHD2_humap)

reorder_ProHD2_humap[ID.x<ID.y ,`:=`(ID.x=ID.y ,
                                     ID.y  = ID.x)]
opencell_validation = merge(reorder_ProHD2_humap,OpenCell[,.(target_gene_name,interactor_gene_name,enrichment)], 
                            by.x = c('ID.x','ID.y'), all.x = T,
      by.y = c('target_gene_name','interactor_gene_name'))
opencell_validation[,validated_opencell := fifelse(!is.na(enrichment),T,F)]
# identifying pairs that both are covered by opencell
opencell_validation[,Protein_1_in_opencell := ID.x %in% OpenCell_coverage,by = .(Protein_1)]
opencell_validation[,Protein_2_in_opencell := ID.y %in% OpenCell_coverage,by = .(Protein_2)]
opencell_validation[,pair_in_opencell := Protein_1_in_opencell & Protein_2_in_opencell]
random_pairs = copy(opencell_validation)
random_pairs= random_pairs[sample(.N,24331)][,class:='random']
opencell_validation = rbind(opencell_validation,random_pairs)
fwrite(opencell_validation,here::here('out','datasets','ProHD2_humap_opencellvalidation.gz'))
file.size(here::here('out','datasets','ProHD2_humap_opencellvalidation.gz'))/(1024 * 1024)
opencell_validation[class !='other',.(perc_opencell = mean(validated_opencell),
                       Total_Number_of_pairs = .N),by = .(class)] |> 
  ggplot(aes(x = class ,y = perc_opencell,fill = class))+
  geom_col()+ theme_bw()+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  ggtitle('Supported by OpenCell',
          subtitle = 'All pairs in hu.MAP3 ProHD2')+
  geom_text(aes(y = 0.21, label = Total_Number_of_pairs))+
  scale_fill_manual(values = values_col)+
  labs(x= 'Score category between ProHD2 and hu.MAP3', y = 'Percentage present in Opencell')
ggsave(here::here('out','plots','humap3_prohd2_opencell_validation.pdf'))

opencell_validation[pair_in_opencell==T &class !='other'][,.(perc_opencell = mean(validated_opencell),
                                          Total_Number_of_pairs = .N),by = .(class)] |> 
  ggplot(aes(x = class ,y = perc_opencell,fill = class))+
  geom_col()+ theme_bw()+
  # theme(legend.position = 'bottom')+
  # guides(fill=guide_legend(nrow=3,byrow=TRUE))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  geom_text(aes(y = 0.31, label = Total_Number_of_pairs))+
  ggtitle('Supported by OpenCell',
          subtitle = 'Only Pairs with proteins covered by Opencell')+
  scale_fill_manual(values = values_col)+
  labs(x= 'Score category between ProHD2 and hu.MAP3', y = 'Percentage present in Opencell')
ggsave(here::here('out','plots','humap3_prohd2_opencell_validation_opencell_universe.pdf'), width = 9, height = 6)

#### complex_averages ####
# reading humap3
humap_3 = fread(here::here('in','datasets','Supplemental_Table_3_hu.MAP3.1_complexes_wConfidenceScores_total15326_20240922.csv'),header = T,  sep = ',')
humap_3[,N_members := stringr::str_count(uniprotACCs,' ')+1, by = clustID]


# coverting into long pairwiseformat 
clusters_to_check = humap_3[N_members>3]
prots_per_complex = clusters_to_check[,.(clustID,uniprotACCs,cluster_confidence)] |> 
  tidyr::separate_longer_delim(uniprotACCs,' ') |> as.data.table()
genes_complex = merge(prots_per_complex,Uniprot_mapping[Type =='Gene_Name',.(Protein,ID)][,head(.SD,1),by = Protein],
                      by.x = 'uniprotACCs',by.y = 'Protein')
genes_complex[,family:=str_sub(ID,1,3),by =uniprotACCs ]
clust_fams = genes_complex[,.(N_fam = .N),by = .(family,clustID)][order(-N_fam)]
clust_fams = clust_fams[,head(.SD,1),by = clustID]

REX1B_ProHD2 = ProHD2_treeclust_RF[Protein_1 =='Q96EN9' |
                                     Protein_2 == 'Q96EN9']
REX1B_ProHD2[,interactor:= fifelse(Protein_1 =='Q96EN9',Protein_2,Protein_1)]
REX1B_ProHD2[,Protein:= str_remove_all(interactor,';[:print:]*$')]
REX1B_ProHD2[,ID:= convert_to_gene(Protein),by = 'Protein']
REX1B_ProHD2[,PSMD13 := (Protein =='Q9UNM6')]
REX1B_ProHD2 = REX1B_ProHD2[order(-RF_covariation_prob)]
which(REX1B_ProHD2$PSMD13==T)/nrow(REX1B_ProHD2)*100
pairs_to_clust = data.table()
for(i in clusters_to_check$clustID){
  print(i)
  clust_tmp = combn(strsplit(clusters_to_check[clustID == i,uniprotACCs],' ') |> unlist(),2) |> 
    t() |> as.data.table()
  level_conf = clusters_to_check[clustID == i,cluster_confidence]
  clust_tmp[,cluster_confidence:=level_conf]
  pairs_to_clust = rbind(pairs_to_clust,clust_tmp[,clustID := i])
  
}
setnames(pairs_to_clust,c('V1','V2'),c('Protein_1','Protein_2'))
pairs_to_clust[Protein_1<Protein_2,`:=`(Protein_1=Protein_2,
                                        Protein_2 = Protein_1)]
pairs_to_clust[Protein_1<Protein_2]
pairs_to_clust[,N_pairs := .N,by = clustID]
pairs_complex_ProHD2 = merge(ProHD2_treeclust_RF,pairs_to_clust, by = c('Protein_1','Protein_2'))
pairs_complex_ProHD2[,N_pairs_joined := .N,by = clustID]
pairs_complex_ProHD2[,coverage:=N_pairs_joined/N_pairs ]

# need to have seen more than 50% of relationships
ProHD2_complexes = pairs_complex_ProHD2[coverage>0.5][,.(complex_covariation  = median(RF_covariation_prob,na.rm = T),
                                                         complex_sd = sd(RF_covariation_prob,na.rm = T),
                                                         N_pairs = .N), by = .(clustID,cluster_confidence)]
# ProHD2_complexes = merge(ProHD2_complexes,humap_3[,.(clustID,uniprotACCs)], by = 'clustID')
# ProHD2_complexes[,Genes := paste(head(convert_to_gene(unlist(strsplit(uniprotACCs,' '))),2),collapse ='_'), by = clustID]
ProHD2_complexes[,confidence :=fifelse(cluster_confidence<4,as.character(cluster_confidence),'>3')]
ProHD2_complexes |> ggplot(aes(x = complex_covariation))+
  geom_histogram(colour =  'grey20', fill = 'grey80', bins = 100)+
  # geom_boxplot()+
  theme_bw()+
  # facet_wrap('cluster_confidence')+
  labs(x = 'median complex Covariation', y = 'Number of complexes', fill ='Cluster Confidence')+
  # scale_colour_manual()
  ggtitle('more than 3 members & more than 50% \nof protein pairs quantified in ProHD2')
ggsave(here::here('out','plots','humap3_prohd2_average_complex_covar_confidence.pdf'),width = 6, height = 6)
ProHD2_complexes |> ggplot(aes(x = complex_covariation, y = factor(as.character(confidence), levels = rev(c('1','2','3','>3'))),
                           fill =  factor(as.character(confidence), levels = rev(c('1','2','3','>3')))))+
  # geom_histogram(colour =  'grey20', bins = 100)+
  geom_boxplot(outliers = F)+
  theme_bw()+
  scale_fill_manual(values = rev(c('grey35','grey60','grey75','grey90')))+
  # facet_wrap('cluster_confidence')+
  labs(x = 'median complex Covariation', y = 'hu.MAP3 Complex Confidence', fill ='Cluster Confidence')+
  # scale_colour_manual()
  ggtitle('more than 3 members & more than 50% \nof protein pairs quantified in ProHD2')
ggsave(here::here('out','plots','humap3_prohd2_average_complex_covar_confidence_boxplot.pdf'),width = 6, height = 6)

fwrite(ProHD2_complexes,here::here('out','datasets','ProHD2_complex_covariation.gz'))
file.size(here::here('out','datasets','ProHD2_complex_covariation.gz'))/(1024 * 1024)

# 
# ProHD2_complexes |> ggplot(aes(x = complex_covariation,y = complex_sd, label = Genes))+
#   geom_point(fill = 'grey80',colour =  'grey20')+theme_bw()+
#   ggrepel::geom_text_repel(data = ProHD2_complexes[complex_covariation>0.8],max.overlaps = 4,size =3)+
#   # theme(element_text(size =5))+
#   labs(x = 'median complex Covariation', y = 'Complex covariation SD')
# ggsave(here::here('out','plots','humap3_prohd2_average_complex_SD.pdf'),width = 6, height = 6)

per_complex = merge(pairs_complex_ProHD2,humap_3[,.(clustID,N_members)], by = 'clustID') 
per_complex = merge(per_complex, Gene_names,all.x = T, by.x = 'Protein_1', by.y = 'Protein')
per_complex = merge(per_complex, Gene_names,all.x = T, by.x = 'Protein_2', by.y = 'Protein')
per_complex[,Type:= factor(fcase(
  clustID %in% c('huMAP3_02362.1','huMAP3_00161.1','huMAP3_00707.1 '),'low covariation',
  clustID %in% c('huMAP3_10595.1','huMAP3_02439.1','huMAP3_09906.1'),'high covariation',
 
  clustID %in% c('huMAP3_06369.1','huMAP3_09363.1','huMAP3_02637.1'),'medium covariation'
),levels= c('low covariation','medium covariation','high covariation'))]
per_complex[,complex_covariation:= mean(RF_covariation_prob), by = clustID]
to_show_complexes = c('huMAP3_00707.1','huMAP3_06369.1', 'huMAP3_09363.1', 'huMAP3_10595.1', 'huMAP3_02439.1')
# per_complex[between(N_members,10,100) & between(complex_covariation,0.4,0.6)][order(complex_covariation,-complex_sd)] |> tail(1000) |> View()
per_complex[clustID %in% to_show_complexes] |> 
  ggplot(aes(x = reorder(clustID,complex_covariation), y= RF_covariation_prob,label = ID.y))+
  geom_boxplot(outliers = F, fill = 'grey80')+
  geom_jitter(alpha = 0.4, width = 0.1, height = 0)+theme_bw()+
  # scale_fill_manual(values = c('high covariation' = 'grey40',
  #                              'medium covariation'= 'grey60',
  #                              'low covariation' = 'grey99'))+
  labs(x = 'huMAP3 complex ID',y = 'subunit pair covariation')+
  theme(legend.position = 'top')+
  # coord_flip()+
  ylim(c(-0.0001,1))+
  ggrepel::geom_text_repel(data =  per_complex[clustID %in% to_show_complexes
                                                             ][,head(.SD,1),by = .(ID.y,clustID)],size = 3.5 ,max.overlaps = 10)+
  # facet_wrap('Type',scales = 'free_x')+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))
ggsave(here::here('out','plots','humap3_per_complex_covar.pdf'),width = 8, height = 12)

BBS_proteins = unique(unlist(per_complex[str_detect(ID.x,'^BBS') & str_detect(ID.y,'^BBS') ][!is.na(Type),.(Protein_1,Protein_2)]))
BBS_RFs= data.table()
for(BBS in BBS_proteins  ){
 BBS_tmp= ProHD2_treeclust_RF[Protein_1 == BBS | Protein_2 ==BBS] 
 BBS_tmp[,Uniprot :=BBS]
 BBS_members = BBS_proteins[BBS_proteins != BBS]
 BBS_tmp[,BBS_member := fifelse(Protein_1 %in% BBS_members | Protein_2 %in% BBS_members,T,F)]
 BBS_RFs = rbind(BBS_RFs,BBS_tmp)
}
BBS_RFs[,IDs := convert_to_gene(Uniprot), by = Uniprot]
BBS_RFs |> ggplot(aes(x = IDs, y = RF_covariation_prob))+
  geom_boxplot(alpha = 0.1)+
  geom_point(data =BBS_RFs[BBS_member == T], colour = 'red' )+
  ggtitle('members of the same complex')
version_proHD2 = 'ProteomeHD_protein_razor_quant_requant_raw_corr_-1.gz'
proHD2 = fread(here::here('in','datasets',version_proHD2))
proHD2[,Uniprot := matching_PGs]
proHD2[Uniprot %like% isoform_pattern, Uniprot := gsub(";.+", "", matching_PGs) ]   # Where all proteins are isoforms, simplify to gene level
proHD2_melt  =proHD2[Uniprot %in% per_complex[!is.na(Type),Protein_1]] |> 
  melt(id.vars = c('Uniprot','matching_PGs' ))
proHD2_melt[,IDs := convert_to_gene(Uniprot), by = Uniprot]
proHD2_melt[str_detect(variable,'^N_')] |> 
  ggplot(aes(y = IDs,x = value))+
  geom_boxplot()+
  ggtitle('N_psms per subunit complex')

proHD2_melt[str_detect(variable,'^Log')] |> 
  ggplot(aes(y = IDs,x = value))+
  geom_boxplot()+
  geom_point(alpha = 0.1)+
  ggtitle('SILAC ratio per subunit complex')

proHD2_melt[str_detect(variable,'^prot_')] |> 
  ggplot(aes(y = IDs,x = value))+
  geom_boxplot()+
  # geom_point(alpha = 0.1)+
  ggtitle('ion agreement per subunit complex')



#proteasomal substructure
complex_id = 'huMAP3_09656.1'
proteasome_subunits = humap_3[clustID == complex_id,uniprotACCs] |> str_split(' ') |> unlist()
ProHD2_treeclust_RF[str_detect(Protein_1,'Q8TAA3')]
proteasome_covar = ProHD2_treeclust_RF[Protein_1 %in% proteasome_subunits &
                                         Protein_2 %in% proteasome_subunits]
proteasome_covar_rev = copy(proteasome_covar)
proteasome_covar = rbind(proteasome_covar,
                         proteasome_covar_rev[,`:=`(Protein_1=Protein_2,
                                                    Protein_2 = Protein_1)])
proteasome_covar = proteasome_covar|> dcast(Protein_1~Protein_2,value.var = 'RF_covariation_prob') |>
  tibble:::column_to_rownames('Protein_1') |> as.matrix()
diag(proteasome_covar)= 1
colnames(proteasome_covar) = colnames(proteasome_covar) |> convert_to_gene()
rownames(proteasome_covar) = rownames(proteasome_covar) |> convert_to_gene()
annot_row = data.table(ID = rownames(proteasome_covar))
annot_row[,subunit:=fcase(
  str_detect(ID,'PSMA'),'PSMA',
  str_detect(ID,'PSMB'),'PSMB',
  str_detect(ID,'PSMC'),'PSMC',
  str_detect(ID,'PSMD'),'PSMD',
  str_detect(ID,'PSME'),'PSME',
  str_detect(ID,'PSMF'),'PSMF',
  str_detect(ID,'PSM',negate = T),'other'
)] 
annot_row = annot_row |> tibble::column_to_rownames('ID')
proteasome_order = c('PSMA1-7
PSMB1-7
PSMA8
PSMB8, PSMB9, PSMB10
PSMC1-6, PSMD1-8, PSMD10-14, ADRM1
PSMG1-4 POMP
PSME1, PSME2, PSME3, PSME3IP1, PSME4
AKIRIN1
AKIRIN2
BAG1
CCDC74A
CCDC74B
CCDC92
FBXO7
PAAF1
PSMF1
REX1BD
RNF181
STKLD1
TMEM31')
proteasome_order  = proteasome_order |> str_split('\n') |>
  unlist()|> str_split(' ') |> unlist() |> str_remove_all(',')
proteasome_order_new = c()
for(i in proteasome_order){
  if(str_detect(i,'-')){
    PSMx = str_sub(i,1,4)
    subs = str_split(i |> str_remove(PSMx),'-') |> unlist() |> 
      as.numeric()
    proteasome_order_new = c(proteasome_order_new,paste0(PSMx,subs[1]:subs[2]))
  }else{
    proteasome_order_new = c(proteasome_order_new,i)
  }
}
proteasome_order_new = proteasome_order_new[proteasome_order_new %in% rownames(annot_row) ]
proteasome_order_new = c(proteasome_order_new,rownames(annot_row)[!(rownames(annot_row) %in% proteasome_order_new)])
pdf(here::here('out','plots','humap_proteasome_intracomplex_no_rowclust.pdf'),width = 11,height = 10)
pheatmap_plot = proteasome_covar[proteasome_order_new,proteasome_order_new] 
diag(pheatmap_plot) <- NA
pheatmap::pheatmap(pheatmap_plot,
                   # annotation_row = annot_row,
                   # annotation_col = annot_row,
                   main = 'huMAP3_09656.1',
                   breaks = seq(0,1.04,0.02),
                   cluster_rows = F, 
                   cluster_cols = F,
                   color =  colorRampPalette(c("#4575b4", "#fefebd", "#d9352a"))(53))
dev.off()

pdf(here::here('out','plots','humap_proteasome_intracomplex_no_rowclust.pdf'),width = 11,height = 10)
pheatmap::pheatmap(proteasome_covar[proteasome_order_new,sort(rownames(annot_row))],
                   # annotation_row = annot_row,
                   annotation_col = annot_row,
                   main = 'huMAP3_09656.1',
                   cluster_rows = F, cluster_cols = T)
# color =  colorRampPalette(c("white", "#f69697", '#F46375', "darkred"))(20))
dev.off()

pdf(here::here('out','plots','humap_proteasome_intracomplex.pdf'),width = 11,height = 10)
pheatmap::pheatmap(proteasome_covar[sort(rownames(annot_row)),sort(rownames(annot_row))],
                   annotation_row = annot_row,
                   # annotation_col = annot_row,
                   main = 'huMAP3_09656.1',
                   cluster_rows = T,
                   cluster_cols = T)
# color =  colorRampPalette(c("white", "#f69697", '#F46375', "darkred"))(20))
dev.off()
fwrite(proteasome_covar |> as.data.frame() |> tibble::rownames_to_column('ID'),
       here::here('out','datasets','proteasome_heatmap.csv'))
file.size(here::here('out','datasets','proteasome_heatmap.csv'))/(1024 * 1024)
##### global map ####
complexes_that_covary = ProHD2_complexes[complex_covariation >0.45,clustID ]
prots_in_complexes  = humap_3[clustID %in% complexes_that_covary & N_members>10,uniprotACCs] |> 
  str_split(' ',simplify = T) |> as.vector() |> unique()
# high_confidences_prots = psm_per_prot[log2(N_psm)>12.5]

ProHD2_complex_members = ProHD2_treeclust_RF[Protein_1 %in% prots_in_complexes &
    Protein_2 %in% prots_in_complexes ]

ProHD2_complex_members = rbind(ProHD2_complex_members[,.(Protein_1,Protein_2,RF_covariation_prob)],
                               ProHD2_complex_members[,.(Protein_1,Protein_2,RF_covariation_prob)
                               ][,`:=`(Protein_1 = Protein_2,
                                       Protein_2 = Protein_1)])
ProHD2_complex_members = ProHD2_complex_members |> 
  dcast(Protein_1 ~ Protein_2, value.var = 'RF_covariation_prob') |> 
  tibble::column_to_rownames('Protein_1') |> as.matrix()
diag(ProHD2_complex_members) = 1
colnames(ProHD2_complex_members) = colnames(ProHD2_complex_members) |> convert_to_gene()
# ProHD2_complex_members[is.na(ProHD2_complex_members)] <- 1
annot_col = data.table(ID = colnames(ProHD2_complex_members) )
annot_col = annot_col[,complex := str_sub(ID,1,3)
][,N_members:= .N, by = complex
][,complex:= fifelse(N_members>6,complex,NA_character_)][,N_members:=NULL] |> 
  tibble::column_to_rownames('ID')
annot_row = data.table(uniprotACCs = rownames(ProHD2_complex_members))
annot_row = merge(annot_row,prots_per_complex, by = 'uniprotACCs')
complexes = data.table(clustID = c('huMAP3_09656.1','huMAP3_09483.1','huMAP3_09307.1',
                                   'huMAP3_14018.1','huMAP3_11832.1','huMAP3_11246.1','huMAP3_11357.1'),
                       Annotation = c('Proteasome','Mitochondrial Ribosome','Mitochondrial ATPases',
                                      'Mitochondrial NDUFs','Minichromosome maintenance complex','Eukaryotic translation initiation','Coatomer complex'))
complexes[,Annotation := paste(Annotation,clustID, sep = '-'), by = clustID]
annot_row = merge(annot_row,complexes, by= 'clustID', all.x = T)
annot_row= annot_row[order(Annotation)][,head(.SD,1), by =uniprotACCs]
annot_row[,N_prots_complex:=.N, by= clustID]
annot_row[is.na(Annotation),Annotation:='other']
annot_row = annot_row[,.(uniprotACCs,Annotation)] |> tibble::column_to_rownames('uniprotACCs')
mycolors = c('#ff9f54','darkred','#ff9289','darkblue','darkgreen','#b6c300','#00d97e','white')

names(mycolors) = annot_row$Annotation |> unique()
mycolors <- list(Annotation = mycolors)
# annot_row = data.table(ID = row.names(ProHD2_complex_members) )
# annot_row = annot_row[,complex := str_sub(ID,1,3)
# ][,N_members:= .N, by = complex
# ][,complex:= fifelse(N_members>5,complex,NA_character_)][,N_members:=NULL] |> 
#   tibble::column_to_rownames('ID')

pdf(here::here('out','plots','humap3_inter_complex.pdf'),width =16,height = 12)

inter_map = pheatmap::pheatmap(ProHD2_complex_members,
                   annotation_row = annot_row  ,
                   annotation_colors = mycolors,
                   # annotation_col = annot_col, 
                   show_rownames = F, show_colnames = F)
dev.off()

ord_names = colnames(ProHD2_complex_members)[inter_map$tree_row$order]
longData <- melt(ProHD2_complex_members)
longData$Var1 <- factor(longData$Var1, levels=rev(names(ord_names)))
longData$Var2 <- factor(longData$Var2, levels=ord_names)
longData = merge(longData,annot_row |> tibble::rownames_to_column('Var1'), by = 'Var1') |> 
  as.data.table()
library(ggh4x)
# colour_annot = names(mycolors$Annotation)
# names(colour_annot) = mycolors$Annotation
longData[,Annotation:= fifelse(Annotation =='other',NA_character_,Annotation)]
ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) +
  # geom_rug()+
  geom_rug(data = longData[,head(.SD),by = .(Var1)], length = unit(0.02, "npc"), sides=  'l', aes(colour = Annotation)) +
  # scale_x_discrete()+
  # scale_x_continuous(expand = c(0.1, 0.1))+
  scale_fill_gradient2(low="#4979b6", mid ='#fafdc7',midpoint = 0.5, high="#d9352a") +
    scale_colour_manual(values= mycolors$Annotation)+
  # labs(x="", y="") +
  # scale_y_dendrogram(hclust = inter_map$tree_row,
  #                    guide = guide_dendro(n.dodge = 2))+
  scale_x_dendrogram(hclust = inter_map$tree_col,position = 'top',
                     guide = guide_dendro(n.dodge = 2),
                     expand = expansion(mult = c(0.027, 0)))+
  theme_bw()+
  # scale_y_reverse()+
  theme(axis.title.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(fill = 'Covariation Score')
  # scale_y_discrete(limits=rev)
ggsave(here::here('out','plots','inter_complex_geomraster.pdf'),width = 11,height=8)

pheamtppheatmap::pheatmap(ProHD2_complex_members[ord,ord],
                   cluster_rows = F,cluster_cols = F,
                   show_rownames = F, show_colnames = F,
                   annotation_row = annot_row  ,
                   annotation_colors = mycolors
)



png(here::here('out','plots','humap3_inter_complex.png'),width =12,height = 8,units= 'in',res = 600)

pheatmap::pheatmap(ProHD2_complex_members,
                   annotation_row = annot_row  ,
                   annotation_colors = mycolors,
                   # annotation_col = annot_col, 
                   show_rownames = F, show_colnames = F)
dev.off()

# data from lineplots

# loading proHD2 to get the canonical prot ratio
version_proHD2 = 'ProteomeHD_protein_razor_quant_requant_raw_corr_-1.gz'
proHD2 = fread(here::here('in','datasets',version_proHD2))
proHD2 =   proHD2[,.SD,.SDcols = names(proHD2) %like% 'matching|^Log']
proHD2_cor = proHD2[,-1] |> as.matrix()
rownames(proHD2_cor) = proHD2$matching_PGs
tmp_medians <- apply( proHD2_cor, 2, median, na.rm = TRUE )  
proHD2_cor <- sweep( proHD2_cor, 2, tmp_medians, FUN = "-" )
complexes_tmp = c('huMAP3_09656.1','huMAP3_09307.1','huMAP3_11832.1')
Interactors_1 = humap_3[clustID == complexes_tmp[1], uniprotACCs] |> unlist() |> strsplit(' ') |> unlist()
Interactors_2 = humap_3[clustID == complexes_tmp[2], uniprotACCs] |> unlist() |> strsplit(' ') |> unlist()
Interactors_3 = humap_3[clustID == complexes_tmp[3], uniprotACCs] |> unlist() |> strsplit(' ') |> unlist()
interactos_present = proHD2_cor[str_detect(rownames(proHD2_cor), 
                                           paste(c(Interactors_1,Interactors_2,Interactors_3),collapse = '|')),] |> 
  is.finite() |> matrixStats::colSums2()
interactos_present = interactos_present>=quantile(interactos_present,probs = seq(0, 1, 0.1))[10]
proHD2_cor_OI_int = proHD2_cor[,interactos_present]
proHD2_cor_OI_int  = proHD2_cor_OI_int[(proHD2_cor_OI_int |> 
                                                is.finite() |> matrixStats::rowSums2())>ncol(proHD2_cor_OI_int)/1.01,]
proHD2_cor_OI_int |> dim()
# tmp_medians <- apply( proHD2_cor_OI_int, 2, median, na.rm = TRUE )  
# proHD2_cor_OI_int <- sweep( proHD2_cor_OI_int, 2, tmp_medians, FUN = "-" )
# proHD2_cor_OI_int |> boxplot()
proHD2_cor_OI_int = proHD2_cor_OI_int |> as.data.table(keep.rownames = 'Uniprot')
# ,[,c(1,270:273,280:295,245:255,325:332,435:440)]
proHD2_cor_OI_int_piv  = proHD2_cor_OI_int[,c(1,70:80,130:150,200:220)] |> 
  melt(id.vars = 'Uniprot', variable.name = 'experiment',value.name = 'LogRatio') |> 
  as.data.table()
proHD2_cor_OI_int_piv[,category:=dplyr::case_when(
  # str_detect(Uniprot,Protein_OI)~Protein_OI,
  str_detect(Uniprot,paste(Interactors_1,collapse = '|'))~complexes_tmp[1],
  str_detect(Uniprot,paste(Interactors_2,collapse = '|'))~complexes_tmp[2],
  str_detect(Uniprot,paste(Interactors_3,collapse = '|'))~complexes_tmp[3],
  TRUE~'other'
  
)]
colour_pallet = c('#00d97e','#00008b','#ff9289')
names(colour_pallet) = c(complexes_tmp)
proHD2_cor_OI_int_piv |> 
  ggplot(aes(x = experiment, y= LogRatio,colour= category,group=Uniprot))+
  geom_line(data = proHD2_cor_OI_int_piv[category != 'other'][category == 'huMAP3_09307.1'],alpha = 0.5)+
  geom_line(data = proHD2_cor_OI_int_piv[category != 'other'][category == 'huMAP3_09656.1'],alpha = 0.5)+
  geom_line(data = proHD2_cor_OI_int_piv[category != 'other'][category == 'huMAP3_11832.1'],alpha = 0.5)+
  # geom_line(data = proHD2_cor_OI_int_piv[category == Protein_OI],alpha = 1)+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(), legend.position = 'bottom',
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylim(-2.5,2.5)+
  scale_colour_manual(values = colour_pallet)+
  annotate('text',label = glue::glue('{complexes_tmp[1]} = Proteaseome \n{complexes_tmp[2]} = ATPsynthase \n{complexes_tmp[3]} = MCM complex'),
           x =10,  y = 2, angle = 0, size = 2.5)+
  ggtitle('Proteasome covaries better with MCM than ATPsynthase')
ggsave(here::here('out','plots','Complex_level_lineplot.pdf'),width = 6,height=6)
