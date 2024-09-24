library(plyr)
library(ggplot2)
library(ggsci)
library(vegan)
library(reshape2)
library(pheatmap)
library(ggplotify)

load('../data/meta_info.RData')
load('../data/RCR95_Phylum.RData')
load('../data/RCR95_Abd.RData')
load('../data/geNomad_res.RData')
load('../phyla_cols.RData')
load('../data/env_cols.RData')

RCR95_Abd$Phylum[is.na(RCR95_Abd$Phylum) | !RCR95_Abd$Phylum %in% names(phyla_cols)] = 'Unclassified'

# Phylum level

Phylum_Abundance = ddply(RCR95_Abd[,-1], .variables = 'Phylum', .fun = function(df){
  return(colSums(df[,-1]))
})
Phylum = Phylum_Abundance$Phylum
Phylum_Abundance = data.frame(t(Phylum_Abundance[,-1]))
colnames(Phylum_Abundance) = Phylum
Phylum_Abundance$Sample_id = gsub('^X','', gsub('\\.','-',rownames(Phylum_Abundance)))
Phylum_Abundance = merge(meta_info[,-c(3,7,8)],Phylum_Abundance, by = 'Sample_id')


## Diversity

Phylum_Shannon = data.frame('Sample_id'=Phylum_Abundance$Sample_id,
                            'Shannon_Index' = apply(Phylum_Abundance[,-c(1:5)], 1, diversity))
Phylum_Shannon = merge(meta_info, Phylum_Shannon)

RCR95_Shannon = data.frame('Sample_id'=gsub('^X','', gsub('\\.','-',colnames(RCR95_Abd)[-c(1:2)])),
                           'Shannon_Index' = apply(RCR95_Abd[,-c(1:2)], 2, diversity))
RCR95_Shannon = merge(meta_info, RCR95_Shannon)

# Figure 3A
envs = c('Station','Store','Bank','Hospital','Street',
         'Soil','Grassland', 'Biofilm','Sediment','Freshwater',
         'Greenhouse','Wastewater','Game animal')
df = RCR95_Shannon[RCR95_Shannon$Type %in% envs,]
df = merge(df, data.frame(table(df$Type)), by.x = 'Type', by.y = 'Var1')
df$Type = apply(df[,c('Type','Freq')], 1, function(x){
  paste0(c(x[1],' (n=',trimws(as.character(x[2])),')'), collapse = '')
})
m = ddply(df, .variables = 'Type', function(x){median(x$Shannon_Index)})
df$Type = factor(df$Type, levels =m$Type[order(m$V1, decreasing = T)])

g = ggplot(data = df, mapping = aes(x = Type, y = Shannon_Index)) +
  geom_boxplot(fill = '#75bde0', outlier.color = '#ff7500', outlier.size = .5) + 
  xlab('Environmental Type') + ylab('Shannon Index') + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90,hjust = 1))


## Figure 3B
df = Phylum_Abundance[Phylum_Abundance$Type %in% envs,-c(1,2,3,5)]
df[,-1] = t(apply(df[,-1], 1, function(x){x/sum(x)}))
df = ddply(df[!is.nan(df[,2]),], .variables = 'Type', .fun = function(x){
  colMeans(x[,-1])
})
df[,-1] = t(apply(df[,-1], 1, function(x){x/sum(x)}))
df = melt(df,id.vars = 1)
df$value = 100*df$value
df$variable = factor(df$variable, levels = names(phyla_cols))
m = df[df$variable == 'Unclassified',]
m = m[order(m$value),]
df$Type = factor(df$Type, levels = m$Type)
m = df[df$variable == 'Kitrinoviricota',]
df$Type = factor(df$Type, levels = m$Type[order(m$value)])
df$variable = factor(df$variable, levels = rev(names(phyla_cols)))

g = ggplot(df,aes(x=Type,y=value,fill=variable))+
  geom_bar(stat = 'identity',position = 'stack')+
  xlab('Environments') + ylab('Relative abundance (%)') + 
  labs (fill="Phyla")+ scale_fill_manual(values = phyla_cols) +
  coord_flip() + theme_bw() + scale_y_continuous(position = 'right')+
  theme(legend.position = 'bottom',
        axis.text = element_text(color='black'))
topptx(g, '../Urban_RNA_Virus_Figs/Viral_Env_Composition.pptx', width = 8,height = 6)

# Figure 3C
load('../data/RCR95_NMDS.RData')
load('../data/RCR95_Anosim.RData')
envs = c('Station','Store','Bank','Hospital','Street',
         'Soil','Grassland', 'Biofilm','Sediment','Freshwater',
         'Greenhouse','Wastewater','Game animal')
df = data.frame(id = 1:length(RCR95_NMDS$species),
                Sample_id = RCR95_NMDS$species, stringsAsFactors = F)
df = merge(df, meta_info[,c(1,4,5)])

RCR95_NMDS$stress
df_points = as.data.frame(RCR95_NMDS$points)
df_points$Sample_id = RCR95_NMDS$species
df_points = merge(df_points, meta_info[,c(1,4,5)])
df_points = df_points[df_points$Type %in% envs,]
m = data.frame(table(df_points$Type))
m = m[order(m$Freq, decreasing = T),]
df_points$Type = factor(df_points$Type, levels = m$Var1)
g = ggplot(data = df_points, mapping = aes(x = MDS1, y = MDS2, color = Type)) +
  geom_point(size = 2) + 
  stat_ellipse(aes(fill=Type),geom="polygon",level=0.95,alpha=0.15)+
  scale_fill_manual(values = env_cols[levels(df_points$Type)])+
  scale_color_manual(values = env_cols[levels(df_points$Type)]) +
  theme_bw()

# Figure 3D 
load('../data/RCR95_PairwiseANOSIM.RData')
df = data.frame()
for(i in 1:length(RCR95_PairwiseANOSIM)){
  pair = unlist(strsplit(names(RCR95_PairwiseANOSIM)[i],'_'))[-2]
  R = RCR95_PairwiseANOSIM[[i]]$statistic
  pval = RCR95_PairwiseANOSIM[[i]]$signif
  df = rbind(df, data.frame(Env1 = pair[1], Env2=pair[2], R = R, pval = pval))
}
df$R[df$R<0] = 0
envs=unique(c(df$Env1, df$Env2))
mat = matrix(NA, nrow = length(envs), ncol = length(envs))
colnames(mat) = rownames(mat) = envs
for(i in 1:nrow(df)){
  env1 = df$Env1[i]
  env2 = df$Env2[i]
  mat[env1, env2] = df$R[i]#as.numeric(df$pval[i] < 0.005)
  mat[env2, env1] = df$R[i]#as.numeric(df$pval[i] < 0.005)
}
diag(mat) = 0
g= as.ggplot(pheatmap(mat,na_col='#DDDDDD',display_numbers = TRUE))


# Figure 3E
df = data.frame(Type = RCR95_Anosim$class.vec, Rank = RCR95_Anosim$dis.rank,
                stringsAsFactors = F)
m = ddply(df, .variables = 'Type', function(x){median(x$Rank)})
m = m[order(m$V1, decreasing = T),]
df$Type = factor(df$Type, levels = m$Type)
g = ggplot(data = df, mapping = aes(x = Type, y = Rank, fill = Type)) +
  geom_boxplot(outlier.color = 'gray', outlier.size = .5) +
  scale_fill_manual(values = c(Between='#f0f0f4',env_cols[levels(df$Type)[-1]]))+
  xlab('') + ylab('Anosim Distance Rank') + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90,hjust = 1))


## Figure 3F
library(ggridges)
envs = c('Station','Store','Bank','Hospital','Street',
         'Soil','Grassland', 'Biofilm','Sediment','Freshwater',
         'Greenhouse','Wastewater','Game animal')
df = Phylum_Abundance[Phylum_Abundance$Type %in% envs,c('Type','urv.p.001','urv.p.002')]
ind = ddply(df, .variables = 'Type', .fun = function(x){
  miss1 = sum(x$urv.p.001 ==0)/nrow(x)
  miss2 = sum(x$urv.p.002 ==0)/nrow(x)
  return(c(miss1, miss2))
})
ind = ind$Type[!(ind$V1 > 0.8 & ind$V2 > 0.8)]
df = df[df$Type %in% ind,]
pval_df = ddply(df, .variables = 'Type', .fun = function(x){
  res = t.test(x[,2], x[3])
  return(res$p.value)
})

df = melt(df)
g = ggplot(df, aes(x=log10(value+1e-5), y=Type, color=variable, point_color=variable, fill=variable)) +
  geom_density_ridges(jittered_points=TRUE, scale = .95, rel_min_height = .01,
                      point_shape = "|", point_size = 3, size = 0.25,
                      position = position_points_jitter(height = 0))+
  scale_y_discrete(expand = c(.01, 0), name = 'Environment') +
  scale_x_continuous(expand = c(0, 0), name = "log10(FPKM + 1E-5)") +
  scale_fill_manual(values = c("#D55E0050", "#0072B250"), labels = c("urv.p.001", "urv.p.002")) +
  scale_color_manual(values = c("#D55E0050", "#0072B250"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E0050", "#0072B250"), guide = "none") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E0050", "#0072B250"),
                                                 color = NA, point_color = NA)))+
  theme_ridges(center = TRUE)


