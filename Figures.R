rm(list=ls())
options(stringsAsFactors = F)

##### codes about Figure 4 can be found under the subdirectory ./DirectionalSelection.
##### Figures of phylogenic trees are visualized using iTOL webtool.


##### Figure 1A
### Figure 1A
# Figure 1A
library(ggplot2)
library(Rmisc)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)

load('../data/mapdata.Rdata')

k=1.6
g = ggplot(data = world) +
  geom_sf(aes(fill = continent,colour = continent))+
  scale_fill_manual(values = c('#e0e0e0','#e0e0e0','#abcec0','#6a9a85','#fbcdaf','#cbd6b8','white','#f59a9b'))+
  scale_color_manual(values = c('#e0e0e0','white','#e2eee9','#a9d6c2','#fef3ec','#f1f0e4','white','#fbdfe0'))+
  theme(legend.position = "none")+
  geom_point(data = loc_shad,aes(x = lng,y = lat),shape = 21, color = '#636e6d',fill = '#636e6d', stroke = 1.5/k,size = 2.5/k,alpha=0.8)+
  #
  geom_point(data = loc,aes(x = lng,y = lat),shape = 21, color = '#fdecee',fill = '#f26975', stroke = 1.5/k,size = 2.5/k,alpha=0.9)+
  # scale_color_identity()+
  # scale_size_identity()+
  #
  theme_bw()+
  theme(
    panel.background = element_rect(fill = "#fcfcfd" ),
    legend.position = 'none',
    legend.text = element_text(colour = 'black',size = 16),
    # axis.text.x = element_text(angle=45, hjust=1, vjust=1),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank(),
    # panel.grid = element_blank(),
    rect = element_rect(fill = "transparent")
  )
g
# Figure 1A
### Figure 1A
##### Figure 1A

###### Figure 1C
### Figure 1C
# Figure 1C
library(vegan)
library(plyr)
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
g
# Figure 1C
### Figure 1C
##### Figure 1C

###### Figure 1D
### Figure 1D
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
  geom_line(aes(color=Env),linewidth = 1) + xlab('Number of samples')+ 
  ylab('Number of RvOTUs') + scale_color_manual(values = env_cols) + 
  scale_fill_manual(values = env_cols) +
  theme_bw()
g
# Figure 1D
### Figure 1D
##### Figure 1D

###### Figure 1G
### Figure 1G
# Figure 1G
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

# write.table(Cooccurrence_df, file = '../Urban_RNA_Virus_Figs/RvOTU_Cooccurence.tsv',
# sep = '\t', row.names = F, col.names = T, quote = F)
Node_attri = data.frame(Envs = c(Cooccurrence_df$Env1, Cooccurrence_df$Env2), 
                        Num = c(Cooccurrence_df$Num1, Cooccurrence_df$Num2), stringsAsFactors = F)
Node_attri = Node_attri[!duplicated(Node_attri),]
# write.table(Node_attri, file = '../Urban_RNA_Virus_Figs/RvOTU_Cooccurence_Node_attri.tsv',
# sep = '\t', row.names = F, col.names = T, quote = F)
# visualized using cytoscape
# Figure 1G
### Figure 1G
##### Figure 1G

###### Figure 2B
### Figure 2B
# Figure 2B
library(sf)
library(reshape2)
library(ggridges)
load('../data/RCR95_Abd.RData')
load('../data/meta_info.RData')
load('../data/mapdata.Rdata')
load('../data/RCR95_cluster.Rdata')
Phylum_Abundance = ddply(RCR95_Abd[,-1], .variables = 'Phylum', .fun = function(df){
  return(colSums(df[,-1]))
})
Phylum = Phylum_Abundance$Phylum
Phylum_Abundance = data.frame(t(Phylum_Abundance[,-1]))
colnames(Phylum_Abundance) = Phylum
Phylum_Abundance$Sample_id = gsub('^X','', gsub('\\.','-',rownames(Phylum_Abundance)))
Phylum_Abundance = merge(meta_info[,-c(3,7,8)],Phylum_Abundance, by = 'Sample_id')
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
g
# Figure 2B
### Figure 2B
##### Figure 2B

###### Figure 2C
### Figure 2C
# Figure 2C
rownames(loc) = loc$City
rownames(meta_info) = meta_info$Sample_id
colnames(RCR95_Abd)[3:ncol(RCR95_Abd)] = gsub('^X','', gsub('\\.','-',colnames(RCR95_Abd)[3:ncol(RCR95_Abd)]))
rownames(RCR95_Abd) = RCR95_Abd$RCR95
np1 = RCR95_Abd[which(RCR95_Abd$Phylum == 'urv.p.001'),3:ncol(RCR95_Abd)]
np2 = RCR95_Abd[which(RCR95_Abd$Phylum == 'urv.p.002'),3:ncol(RCR95_Abd)]

np1[np1 <= 10] = 0
s_ = apply(np1,2,sum)
np1 = np1[,s_ != 0]
np2[np2 <= 10] = 0
s_ = apply(np2,2,sum)
np2 = np2[,s_ != 0]

m1 = meta_info[colnames(np1),]
m2 = meta_info[colnames(np2),]

id1 = clu$V2[which(clu$V1 %in% RCR95_Abd$RCR95[which(RCR95_Abd$Phylum == 'urv.p.001')])]
id2 = clu$V2[which(clu$V1 %in% RCR95_Abd$RCR95[which(RCR95_Abd$Phylum == 'urv.p.002')])]

tf = function(x){
  return(unlist(strsplit(x,'[|]'))[1])
}
m1 = meta_info[unique(sapply(id1,tf)),]
m2 = meta_info[unique(sapply(id2,tf)),] 

c = unique(m1$City)  

df1 = loc[intersect(loc$City,unique(m1$City)),c('City','lng','lat')]
df1$num = table(m1$City)[df1$City]
df1$type = 'p001'
df2 = loc[intersect(loc$City,unique(m2$City)),c('City','lng','lat')]
df2$num = table(m2$City)[df2$City]
df2$type = 'p002'
df = rbind(df1,df2)
#######
world <- ne_countries(scale = "medium", returnclass = "sf")
projection <- "+proj=ortho +lat_0=85 +lon_0=5"

world_proj <- world %>% 
  st_transform(crs = projection)

loc_proj <- df %>% 
  st_as_sf(coords = c("lng", "lat"), crs = 4326) %>% 
  st_transform(crs = projection) %>% 
  cbind(st_coordinates(.))

earth_radius <- 6371000
circle <- st_point(x = c(0, 0)) %>%
  st_buffer(dist = earth_radius) %>%
  st_sfc(crs = projection)
g = ggplot() +
  geom_sf(data = circle, fill = "#A6CAE0", color = NA) +
  geom_sf(data = world_proj, fill = "#C5D4D0", color = "#C5D4D0", size = 0.025) +
  geom_point(data = loc_proj, aes(x = X, y = Y, color = type), size = 1) +
  theme_void() +  
  theme(panel.background = element_rect(fill = "transparent"))
g
# Figure 2C
### Figure 2C
##### Figure 2C

##### Figure 3E
### Figure 3E
# Figure 3E
load('../data/sam_count.Rdata')
sam_count$type = factor(sam_count$type,levels = c('pos','dup','pisu.u.c.1'))
g = ggplot(sam_count,aes(y=logfold,x=type,color=type))+
  geom_point(position = position_jitter(width = 0.3))+
  scale_color_manual(values = c('#2E75B6','#C00000','#2eb19f'))+
  theme_bw()+
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),  # 面板背景透明
    plot.background = element_rect(fill = "transparent", colour = NA),   # 图形整体背景透明
    #legend.title = element_blank(),
    #legend.position = 'none',
    # legend.text = element_text(colour = 'black',size = 16),
    # axis.text = element_blank(),
    #axis.text.x = element_blank(),
    axis.title = element_blank(),
    #panel.border = element_blank(),
    #panel.grid = element_blank(),
    # rect = element_rect(fill = "transparent"),
  )
g
# Figure 3E
### Figure 3E
##### Figure 3E

###### Figure 5A
### Figure 5A
# Figure 5A
df = Phylum_Abundance[Phylum_Abundance$Type %in% envs,-c(1,2,3,5)]
df[,-1] = t(apply(df[,-1], 1, function(x){x/sum(x)}))
df = ddply(df[!is.nan(df[,2]),], .variables = 'Type', .fun = function(x){
  colMeans(x[,-1])
})
df[,-1] = t(apply(df[,-1], 1, function(x){x/sum(x)}))
df = melt(df,id.vars = 1)
df$value = 100*df$value
df$variable[which(!(df$variable %in% names(phyla_cols)))] = 'Unclassified'
df$variable = factor(df$variable, levels = names(phyla_cols))
m = df[df$variable == 'Unclassified',]
m = m[order(m$value),]
df$Type = factor(df$Type, levels = m$Type)
m = df[df$variable == 'Kitrinoviricota',]
df$Type = factor(df$Type, levels = m$Type[order(m$value)])
df$variable = factor(df$variable, levels = rev(names(phyla_cols)))
df$variable[which(!(df$variable %in% names(phyla_cols)))] = 'Unclassified'

g = ggplot(df,aes(x=Type,y=value,fill=variable))+
  geom_bar(stat = 'identity',position = 'stack')+
  xlab('Environments') + ylab('Relative abundance (%)') + 
  labs (fill="Phyla")+ scale_fill_manual(values = phyla_cols) +
  coord_flip() + theme_bw() + scale_y_continuous(position = 'right')+
  theme(legend.position = 'bottom',
        axis.text = element_text(color='black'))
g
# Figure 5A
### Figure 5A
##### Figure 5A

###### Figure 5B
### Figure 5B
# Figure 5B
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
g
# Figure 5B
### Figure 5B
##### Figure 5B

###### Figure 5C
### Figure 5C
# Figure 5C
df = data.frame(Type = RCR95_Anosim$class.vec, Rank = RCR95_Anosim$dis.rank,
                stringsAsFactors = F)
m = ddply(df, .variables = 'Type', function(x){median(x$Rank)})
m = m[order(m$V1, decreasing = T),]
df$Type = factor(df$Type, levels = m$Type)
g = ggplot(data = df, mapping = aes(x = Type, y = Rank ,fill = Type)) +
  geom_boxplot(outlier.color = 'gray', outlier.size = .5) +
  scale_fill_manual(values = c('#63c2d2',rep('#d6dbe0',13)))+
  scale_y_continuous(breaks=seq(0, 1250000, by=500000))+
  xlab('') + ylab('Anosim Distance Rank') + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90,hjust = 1))
g
# Figure 5C
### Figure 5C
##### Figure 5C

###### Figure 5D
### Figure 5D
# Figure 5D
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
g = as.ggplot(pheatmap(mat,na_col='#DDDDDD',display_numbers = TRUE))
# visualized using cytoscape
# Figure 5D
### Figure 5D
##### Figure 5D

###### Figure 6A
### Figure 6A
# Figure 6A
library(plyr)
load('../data/phyla_cols.RData')

host_cols = c(Bacteria = "#DC0000FF",Fungi="#177cb0",Plants="#00A087FF",Vertebrates="#a78e44",
              Invertebrates="#F39B7FFF",Algae="#B07AA1",Others='#C7C7C7FF')

load('../data/RCR95.RData')
load('../data/RCR95_Phylum.RData')
load('../data/RdRP_Host.RData')
RdRP_Host = RdRP_Host[RdRP_Host$RdRP %in% RCR95$RCR95,]
RdRP_Host = RdRP_Host[RdRP_Host$Identity >= 50,]
RdRP_Host = merge(RdRP_Host, RCR95_Phylum[,c('RCR95','Phylum')], 
                  by.x = 'RdRP', by.y = 'RCR95', all.x = T)
RdRP_Host$Phylum.x[RdRP_Host$Phylum.x == ''] = RdRP_Host$Phylum.y[RdRP_Host$Phylum.x == '']
RdRP_Host = RdRP_Host[,-16]
colnames(RdRP_Host)[5] = 'Phylum'
RdRP_Host$Phylum[RdRP_Host$Phylum == '' | is.na(RdRP_Host$Phylum)] = 'Unclassified'
host_list = unique(RdRP_Host$Host.source.revised)

RdRP_Host = ddply(RdRP_Host, .variables = 'RdRP', function(df){
  phylum = table(df$Phylum)
  phylum = phylum[order(phylum, decreasing = T)]
  phylum = names(phylum)[1]
  host = paste0(unique(df$Host.source.revised), collapse = ';')
  virion_host = paste0(unique(df$virion_host), collapse = ';')
  return(c(Phylum = phylum, Host = host, Virion=virion_host))
})
#host_list = c("Plants","Fungi","Others","Marine","Algae","Bacteria", "Invertebrates","Protist","Vertebrates")

host_freq = unlist(lapply(host_list, function(host, RdRP_Host){
  if(host != 'vertebrates'){
    sum(grepl(host, RdRP_Host$Host))
  }else{
    sum(grepl(host, RdRP_Host$Host) & !grepl(paste('in',host, sep = ''), RdRP_Host$Host))
  }
}, RdRP_Host))
host_freq = data.frame(Host=host_list, Freq = host_freq/nrow(RdRP_Host), stringsAsFactors = F)
host_freq = host_freq[order(host_freq$Freq, decreasing = T),]

host = strsplit(RdRP_Host$Host, split = ';')
n = unlist(lapply(host, length))
df = data.frame(Phylum = rep(RdRP_Host$Phylum, n), Host = unlist(host),stringsAsFactors = F)
df = data.frame(table(df))
df = df[df$Freq > 0,]
n1 = ddply(df, .variables = 'Phylum', function(x){sum(x$Freq)})
n1 = n1[order(n1$V1, decreasing = T),]
n2 = ddply(df, .variables = 'Host', function(x){sum(x$Freq)})
n2 = n2[order(n2$V1, decreasing = T),]
n2 = n2[c(1:2,4:7,3),]
df$Phylum = factor(df$Phylum,levels = n1$Phylum)
df$Host = factor(df$Host,levels = n2$Host)

library(ggalluvial)
library(eoffice)
g = ggplot(df, aes(y= Freq, axis1 = Phylum, axis2 = Host)) +
  geom_alluvium(width= 1/12) +
  geom_stratum(width = 1/12, fill = c(phyla_cols[rev(levels(df$Phylum))],
                                      host_cols[rev(levels(df$Host))])) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size = 3) +
  scale_x_discrete(limits= c("Phylum", "Host"), expand = c(.05, .05)) +
  ylab("Number viruses in RCR95 set") + theme_bw()
g
# Figure 6A
### Figure 6A
##### Figure 6A

###### Figure 6B
### Figure 6B
# Figure 6B
library(plyr)
load('../data/RCR95_Phylum.RData')
load('../data/IMGspacer_match.RData')
load('../data/IMGspacers.RData')

RCR95_Phylum$RCR95 = gsub('\\.[0-9]*$','',gsub('_[0-9]*$','',RCR95_Phylum$RCR95))
RCR95_Phylum = RCR95_Phylum[!duplicated(RCR95_Phylum$RCR95),]

img_spacer_num = table(IMGspacers$gtdb_lineage)
img_spacer_num = img_spacer_num[img_spacer_num > 2]
IMGspacers = IMGspacers[IMGspacers$gtdb_lineage %in% names(img_spacer_num),]
IMGspacer_match = IMGspacer_match[IMGspacer_match$mismatch <= 1 & 
                                    IMGspacer_match$length > 20 &
                                    IMGspacer_match$gapopen <= 1,]
IMGspacer_match$qacc = gsub('^spacer_','', IMGspacer_match$qacc)
#IMGspacer_match = IMGspacer_match[IMGspacer_match$sacc %in% RCR95_Phylum$RCR95,]

IMGspacer_match = merge(RCR95_Phylum[!duplicated(RCR95_Phylum$Contig),c('RCR95','Phylum')], IMGspacer_match[,1:2], 
                        by.x = 'RCR95', by.y = 'sacc')

IMGspacer_match = merge(IMGspacer_match,IMGspacers, by.x = 'qacc',by.y = 'spacer_id')

colnames(IMGspacer_match) = c('Spacer_id','RCR95','Virus_Phylum','Spacer_lineage')
IMGspacer_match$Virus_Phylum[is.na(IMGspacer_match$Virus_Phylum)] = 'Unknown'
#IMGspacer_match = IMGspacer_match[!duplicated(IMGspacer_match[,c('RCR95','Spacer_id')]),]

IMGspacer_match$Spacer_Kingdom = unlist(lapply(strsplit(IMGspacer_match$Spacer_lineage, split = ';'), function(x){
  x = gsub('_[A-Z]$','',gsub('^d__','',grep('^d__', x, value = T)))
  if(length(x) == 0){x = ''}
  return(x)
}))
IMGspacer_match$Spacer_Phylum = unlist(lapply(strsplit(IMGspacer_match$Spacer_lineage, split = ';'), function(x){
  x = gsub('_[A-Z]$','',gsub('^p__','',grep('^p__', x, value = T)))
  if(length(x) == 0){x = ''}
  return(x)
}))
IMGspacer_match$Spacer_Class = unlist(lapply(strsplit(IMGspacer_match$Spacer_lineage, split = ';'), function(x){
  x = gsub('_[A-Z]$','',gsub('^c__','',grep('^c__', x, value = T)))
  if(length(x) == 0){x = ''}
  return(x)
}))
IMGspacer_match$Spacer_Genus = unlist(lapply(strsplit(IMGspacer_match$Spacer_lineage, split = ';'), function(x){
  x = gsub('_[A-Z]$','',gsub('^g__','',grep('^g__', x, value = T)))
  if(length(x) == 0){x = ''}
  return(x)
}))
IMGspacer_match$Spacer_Species = unlist(lapply(strsplit(IMGspacer_match$Spacer_lineage, split = ';'), function(x){
  x = gsub('_[A-Z]$','',gsub('^s__','',grep('^s__', x, value = T)))
  if(length(x) == 0){x = ''}
  return(x)
}))
IMGspacer_match = IMGspacer_match[!duplicated(IMGspacer_match[,c('RCR95','Spacer_id')]),]
nrow(IMGspacer_match)
length(unique(IMGspacer_match$RCR95))
length(unique(IMGspacer_match$Spacer_id))
length(unique(IMGspacer_match$Spacer_Species[IMGspacer_match$Spacer_Kingdom=='Bacteria']))
length(unique(IMGspacer_match$Spacer_Species[IMGspacer_match$Spacer_Kingdom=='Archaea']))

RCR95_host_Freq = ddply(IMGspacer_match, .variables = 'RCR95', function(df){
  length(unique(df$Spacer_Species))
})
sum(RCR95_host_Freq$V1==1)/nrow(RCR95_host_Freq)
sum(RCR95_host_Freq$V1 <= 3)/nrow(RCR95_host_Freq)
sum(RCR95_host_Freq$V1 > 10)/nrow(RCR95_host_Freq)


library(circlize)
library(tidygraph)
df = IMGspacer_match[IMGspacer_match$Spacer_Species != '',]
df = df[!duplicated(df[,c('RCR95','Spacer_Species')]),]
Virus_Phylum = unique(df$Virus_Phylum)
Spacer_Phylum = unique(df$Spacer_Phylum)
mat = matrix(0, nrow = length(Spacer_Phylum), ncol= length(Virus_Phylum))
rownames(mat) = Spacer_Phylum
colnames(mat) = Virus_Phylum
for(v in Virus_Phylum){
  x = table(df$Spacer_Phylum[df$Virus_Phylum == v])
  x = x[names(x) %in% rownames(mat)]
  mat[names(x),v] = as.numeric(x)
}
mat = mat[order(rowSums(mat), decreasing = T),]
mat = rbind(mat[rowSums(mat)>=100,], colSums(mat[rowSums(mat)<100,]))
rownames(mat)[nrow(mat)] = 'Others'
mat = mat[,order(colSums(mat))]
mat = mat[,c(3,1,2,4:7)]
grid.col = rep('grey',ncol(mat)+nrow(mat))
names(grid.col) = c(rownames(mat),rev(colnames(mat)))
grid.col[names(phyla_cols)] = phyla_cols

circos.par(start.degree = 270)
chordDiagram(t(mat), annotationTrack = "grid", preAllocateTracks = 1, grid.col = grid.col)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", 
              cex = .5,niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
circos.clear()
# Figure 6B
### Figure 6B
##### Figure 6B

###### Figure 6C
### Figure 6C
# Figure 6C
ESKAPE = c('Enterenecus faecium','Staphylococcus aureus','Klebsiella pneumoniae','Acinetobacter baumannii',
           'Pseudomonas aeruginosa','Enterobacter')
df = IMGspacer_match[IMGspacer_match$Spacer_Species %in% ESKAPE |
                       grepl('Enterobacter',IMGspacer_match$Spacer_Species),]
df$Spacer_Species[grepl('Enterobacter',df$Spacer_Species)] = 'Enterobacter spp.'
df = df[,c('RCR95','Virus_Phylum','Spacer_Species')]
df = df[!duplicated(df[,c(1,3)]),]
edges = df[,c(1,3)]
colnames(edges) = c('from','to')
nodes = data.frame(name = unique(union(edges$from, edges$to)))
nodes$type = 'Bacteria'
nodes$type[grepl('_',nodes$name)] = 'Virus'
nodes = merge(nodes, RCR95_Phylum[,c('RCR95','Phylum')], by.x = 'name', by.y = 'RCR95', all.x =T)
nodes$Phylum[is.na(nodes$Phylum) & nodes$type == 'Virus'] = 'Unknown'
nodes$Phylum[is.na(nodes$Phylum)] = nodes$name[is.na(nodes$Phylum)]
table(edges$to)
write.table(edges, file = '../Urban_RNA_Virus_Figs/IMG_ESKAPE_edges.tsv', row.names = F, sep = '\t',quote = F)
write.table(nodes, file = '../Urban_RNA_Virus_Figs/IMG_ESKAPE_nodes.tsv', row.names = F, sep = '\t',quote = F)
# visualized using cytoscape
# Figure 6C
### Figure 6C
##### Figure 6C
