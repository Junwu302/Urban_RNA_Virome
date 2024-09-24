# Figure 1B
# The overlap relationships were calculated as follow and 
# the visualized manually.
library(VennDiagram) 
load('../data/RvOTU_Venn.RData')
venn.diagram(
  x = RvOTU_Venn,
  category.names =names(RvOTU_Venn),
  col=c('#D9D9D9','#7FB27F','#F8AAAC','#7030A0','#ACCDE9'),
  fill = c('#D9D9D9','#7FB27F','#F8AAAC','#7030A0','#ACCDE9'),
  alpha = 0.3,
  filename = 'venn.png',
  output=TRUE
)


# Figure 1C
library(vegan)
library(plyr)
library(ggplot2)
library(ggalluvial)
library(ggthemes)
library(eoffice)
load('../data/RvOTU_Accum.RData')
load('../data/env_cols.RData')
df_all = data.frame(x = RvOTU_Accum$All$sites,
                    y = RvOTU_Accum$All$richness,
                    sd = RvOTU_Accum$All$sd, stringsAsFactors = F)
df_all$lower = df_all$y - df_all$sd
df_all$upper = df_all$y + df_all$sd
df_all$lower[df_all$lower < 0] = 0
g = ggplot(df_all, aes(x, y)) + 
  geom_ribbon(aes(ymin = lower, ymax=upper, x = x), fill = 'gray', alpha = 0.8) +
  geom_line(color ='black') + xlab('Number of samples')+ 
  ylab('Number of RvOTUs') + theme_bw()

# Figure 1D
df_envs = data.frame()
for(env in names(RvOTU_Accum)[-1]){
  df = data.frame(x = RvOTU_Accum[[env]]$sites,
                  y = RvOTU_Accum[[env]]$richness,
                  sd = RvOTU_Accum[[env]]$sd, stringsAsFactors = F)
  df$lower = df$y - df$sd
  df$upper = df$y + df$sd
  df$lower[df$lower < 0] = 0
  df$Env = env
  df_envs = rbind(df_envs, df)
}
df_envs = df_envs[df_envs$Env %in% c('Station','Store','Bank','Hospital','Soil',
                                     'Biofilm','Sediment','Wastewater','Game animal'),]
g = ggplot(df_envs, aes(x, y)) + geom_ribbon(aes(ymin = lower, ymax=upper, x = x,fill = Env),alpha = 0.3) +
  geom_line(aes(color=Env),size = 1) + xlab('Number of samples')+ 
  ylab('Number of RvOTUs') + scale_color_manual(values = env_cols) + 
  scale_fill_manual(values = env_cols) +
  theme_bw()

# Figure 1E
load('../data/RvOTU_abd.RData')
load('../data/meta_info.RData')
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
# visualized using cytoscape

