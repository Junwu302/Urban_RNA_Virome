# Co-occurence
load('../Urban_RNA_Virus_Data/RvOTU_abd.RData')
load('../Urban_RNA_Virus_Data/meta_info.RData')
vOTUs = rownames(RvOTU_abd)
all_samples = gsub('\\.','-',gsub('^X','',colnames(RvOTU_abd)))

# High-frequency RvOTUs
envs = table(meta_info$Type[meta_info$Type != ''])
envs = envs[order(envs, decreasing = T)]
envs = names(envs)[envs >= 30]
game_samples = meta_info$Sample_id[meta_info$DataSet == 'gameanimal']

HighFreq_RvOTU = data.frame()
sample_num = c()
for(env in envs){
  print(env)
  sample_id = meta_info$Sample_id[meta_info$Type == env]
  df = RvOTU_abd[,all_samples %in% sample_id]
  if(ncol(df)>30){
    sample_num[env] = ncol(df)
    freq = apply(df, 1, function(x){
      sum(x >= 1)
    })/ncol(df)
    votus = names(freq)[freq >= 0.1]
    HighFreq_RvOTU = rbind(HighFreq_RvOTU, data.frame(Env = env, RvOTU = votus, stringsAsFactors = F))
  }
}
envs = unique(HighFreq_RvOTU$Env)
sample_num= sample_num[envs]
Cooccurrence_df = data.frame(t(combn(envs,2)),stringsAsFactors = F)
colnames(Cooccurrence_df) = c('Env1','Env2')
Cooccurrence_df$Num1 = NA
Cooccurrence_df$Num2= NA
Cooccurrence_df$Common= NA
Cooccurrence_df$JI = NA
Cooccurrence_df$pval= NA
Cooccurrence_df$JI = NA
Cooccurrence_df$pval= NA
all_highrvotu = length(unique(HighFreq_RvOTU$RvOTU))

for(i in 1:nrow(Cooccurrence_df)){
  a = HighFreq_RvOTU$RvOTU[HighFreq_RvOTU$Env==Cooccurrence_df$Env1[i]]
  b = HighFreq_RvOTU$RvOTU[HighFreq_RvOTU$Env==Cooccurrence_df$Env2[i]]
  c = intersect(a, b)
  k = c(length(a), length(b))
  k = k[order(k)]
  p = phyper(length(c)-1, k[2], all_highrvotu-k[1], k[1], lower.tail = F)
  Cooccurrence_df$Num1[i] = length(a)
  Cooccurrence_df$Num2[i]= length(b)
  Cooccurrence_df$Common[i]= length(c)
  Cooccurrence_df$JI[i] = length(c)/(sum(k)-length(c))
  Cooccurrence_df$pval[i] = p
}
Cooccurrence_df$qval = p.adjust(Cooccurrence_df$pval)

Cooccurrence_df$Cooccurence = 'No'
Cooccurrence_df$Cooccurence[Cooccurrence_df$qval < 0.001 & Cooccurrence_df$JI >= 0.1] = 'Yes'
write.table(Cooccurrence_df, file = '../Urban_RNA_Virus_Figs/RvOTU_Cooccurence.tsv',
            sep = '\t', row.names = F, col.names = T, quote = F)
Node_attri = data.frame(Envs = c(Cooccurrence_df$Env1, Cooccurrence_df$Env2), 
                        Num = c(Cooccurrence_df$Num1, Cooccurrence_df$Num2), stringsAsFactors = F)
Node_attri = Node_attri[!duplicated(Node_attri),]
write.table(Node_attri, file = '../Urban_RNA_Virus_Figs/RvOTU_Cooccurence_Node_attri.tsv',
            sep = '\t', row.names = F, col.names = T, quote = F)
# network is visualized using Cytoscape

