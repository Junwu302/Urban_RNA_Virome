library(reshape2)
library(eoffice)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(plyr)
load('../data/env_cols.RData')
load('../data/Kraken_Res.RData')
load('../data/meta_info.RData')
load('../data/phyla_cols.RData')
# Human pathogenic RNA viruses https://doi.org/10.1073/pnas.2121335119 Table S1
human_viruses = read.table('../Urban_RNA_Virus_Data/Human pathogenic RNA viruses.tsv',
                           sep = '\t',header = T, stringsAsFactors = F)
human_viruses$Family = paste('f__',human_viruses$Family,sep = "")

ind = unlist(sapply(human_viruses$Kraken_Name, function(x, ids){
  grep(paste(x,'$',sep=''), ids)
}, colnames(Kraken_Res)))
Kraken_Res = Kraken_Res[,c(1,ind)]
colnames(Kraken_Res) = gsub('.*\\|','', colnames(Kraken_Res))
Kraken_Res[is.na(Kraken_Res)] = 0
for(i in 2:ncol(Kraken_Res)){
  Kraken_Res[,i] = Kraken_Res[,i]/Kraken_Res[,1]
}
Kraken_Res = Kraken_Res[,-1]
Kraken_Res$Sample_id = rownames(Kraken_Res)
Kraken_Res = Kraken_Res[,c(26,1:25)]
\
Kraken_Res = merge(meta_info[,c('Sample_id','Country','Geo_loc','Type')], Kraken_Res)
Env_num = data.frame(table(Kraken_Res$Type))
envs = c('Station','Store','Bank','Hospital','Street',
         'Soil','Grassland', 'Biofilm','Sediment','Freshwater',
         'Greenhouse','Wastewater','Game animal')


## Figure 6A
df = Kraken_Res[Kraken_Res$Type %in% envs,-c(1:3)]
df$Abundance = rowSums(df[,-1])
df =df[,c('Type','Abundance')]
m = ddply(df, .variables = 'Type', .fun = function(x){
  x = log10(x$Abundance)
  median(x[!is.infinite(x)])
})
m$V1[is.na(m$V1)] = -10
m = m[order(m$V1, decreasing = T),]
df$Type = factor(df$Type, levels = m$Type)

g = ggplot(data = df, mapping = aes(x=Type, y = Abundance, fill = Type)) +
  geom_boxplot() + scale_y_log10() + scale_fill_manual(values = env_cols) + 
  theme_bw() + xlab('') + ylab('Sum of relative abundance') + theme(legend.position = 'none',
                                                                    axis.text.x = element_text(angle = 45,hjust = 1))

# Figure 6B
df = Kraken_Res[Kraken_Res$Type %in% envs,-c(1:3)]
df = ddply(df, .variables = 'Type', .fun = function(x){
  colMeans(x[,-1])
})
library(pheatmap)
library(ggplotify)
mat = t(scale(df[,-1], center = T))
colnames(mat) = df$Type
rownames(mat) = gsub('^s__','', rownames(mat))
transmission = human_viruses[,c('Kraken_Name','Transmission')]
transmission$Kraken_Name = gsub('^s__','',transmission$Kraken_Name)
transmission = transmission[transmission$Kraken_Name %in% rownames(mat),]
transmission = transmission[order(transmission$Transmission),]
transmission = data.frame(Transmission = transmission$Transmission,
                          row.names = transmission$Kraken_Name, stringsAsFactors = F)
mat = mat[rownames(transmission),]
g = as.ggplot(pheatmap(mat,cluster_rows = F,annotation_row = transmission,
                       gaps_row = c(10, 11,12,18)))


# Figure 6C
total_pre =  apply(Kraken_Res[,-c(1:4)], 2, function(x){sum(x>0.0001)})
total_pre = data.frame('Kraken_Name'= names(total_pre), 'Freq' = as.numeric(total_pre),
                       stringsAsFactors = F)
#total_pre = merge(total_pre, human_viruses[,1:3])
total_pre = total_pre[order(total_pre$Freq, decreasing = T),]
total_pre = rbind(total_pre[1:15,], 
                  data.frame('Kraken_Name'='Others',Freq = sum(total_pre$Freq[16:nrow(total_pre)])))
total_pre$Kraken_Name = gsub('^s__','', total_pre$Kraken_Name)
total_pre$Kraken_Name = factor(total_pre$Kraken_Name, levels = rev(total_pre$Kraken_Name))
g1 = ggplot(total_pre, mapping = aes(x = Kraken_Name, y = Freq, fill = Kraken_Name)) +
  geom_bar(stat = 'identity')+ geom_text(aes(label=Freq),size=4,hjust=-.2)+
  scale_fill_tableau('Tableau 20',direction  = -1) + coord_flip()+
  theme_bw() + theme(legend.position = 'none')

## SARS-COV2
ind =  Kraken_Res$`s__Severe_acute_respiratory_syndrome-related_coronavirus`>0.0001
df = meta_info[meta_info$Sample_id %in% Kraken_Res$Sample_id[ind],]
df = data.frame(table(df$Type))
df$Var1 = as.character(df$Var1)
df = df[order(df$Freq, decreasing = T),]
df = rbind(df[df$Freq >=3,], data.frame(Var1='Others',Freq=sum(df$Freq[df$Freq <3])))
df$Var1[df$Var1 == ''] = 'Missing'
df$Var1 = factor(df$Var1, levels = df$Var1)
g2 = ggplot(df, mapping = aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = 'identity')+ geom_text(aes(label=Freq),size=4,vjust=-.2)+
  scale_fill_tableau('Tableau 20') + xlab('') + ylab('Number of Samples') +
  theme_bw() + theme(legend.position = 'none',
                     axis.text.x = element_text(angle = 45,hjust = 1))
g = g1 + theme(legend.position = "none")+
  annotation_custom(grob = ggplotGrob(g2), xmin = 1, xmax = 10, ymin = 100, ymax = 650)



# Figure 6D
Prevalence = data.frame()
m = apply(Kraken_Res[,-c(1:4)], 2, function(x){
  sum(x > 0.0001)
})
n = nrow(Kraken_Res)
for(env in envs){
  df = Kraken_Res[Kraken_Res$Type == env,-c(1:4)]
  if(nrow(df) < 30){next}
  k = nrow(df)
  q = apply(df, 2, function(x){sum(x > 0.0001)})
  pvalue = rep(NA, length(k))
  for(i in 1:length(q)){
    v = names(q)[i]
    pvalue[i] = phyper(q[v]-1, m[v], n, k, lower.tail=F)
  }
  df  = data.frame('Kraken_Name'=gsub('^s__','',names(q)), "sample_num" = nrow(df),
                   'observed'= as.numeric(q), pvalue = pvalue, stringsAsFactors = F)
  df$Type = env
  Prevalence = rbind(Prevalence, df)
}
Prevalence = Prevalence[Prevalence$observed > 0,]
Prevalence$Enrichment = 'No'
Prevalence$Enrichment[Prevalence$pvalue < 0.01] = 'Yes'
n = ddply(Prevalence, .variables = 'Kraken_Name', .fun = function(x){
  sum(x$Enrichment == 'Yes')
})
n = n[order(n$V1),]
Prevalence$Kraken_Name = factor(Prevalence$Kraken_Name, levels = n$Kraken_Name)
Prevalence$Enrichment = factor(Prevalence$Enrichment, levels = c('Yes','No'))
Prevalence$Prevalence = Prevalence$observed/Prevalence$sample_num
Prevalence$sample_num = as.character(Prevalence$sample_num)
Prevalence$Type = apply(Prevalence[,c('Type','sample_num')], 1, function(x){
  paste0(c(x[1],' (', x[2],')'),collapse = '')
})
n = ddply(Prevalence, .variables = 'Type', .fun = function(x){
  sum(x$Enrichment == 'Yes')
})
n = n[order(n$V1, decreasing = T),]

Prevalence$Type = factor(Prevalence$Type, levels = n$Type[c(1,5,4,2,3,6:10)])
Prevalence$Kraken_Name = factor(Prevalence$Kraken_Name, 
                                levels = levels(Prevalence$Kraken_Name)[c(1:16,20,19,18,17,22,21,23)])

g = ggplot(data = Prevalence, mapping = aes(x = Type, y = Kraken_Name)) +
  geom_point(aes(size = Prevalence, color = Enrichment)) +
  scale_color_manual(values = c('#DC0000FF','#C7C7C7FF')) + xlab('') + ylab('') +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust = 1))
















