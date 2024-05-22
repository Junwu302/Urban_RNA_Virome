library(plyr)
phyla_cols = c(Kitrinoviricota='#DC0000FF',Pisuviricota='#4DBBD5FF', Duplornaviricota='#F28E2B',
               Lenarviricota='#3C5488FF',Negarnaviricota='#67bf5c',urv.p.001='#B07AA1',
               urv.p.002='#a78e44', Unclassified='#C7C7C7FF')
host_cols = c(Bacteria = "#DC0000FF",Fungi="#177cb0",Plants="#00A087FF",Vertebrates="#a78e44",
              Invertebrates="#F39B7FFF",Algae="#B07AA1",Others='#C7C7C7FF')
load('../Urban_RNA_Virus_Data/RCR95.RData')
load('../Urban_RNA_Virus_Data/RCR95_Phylum.RData')
load('../Urban_RNA_Virus_Data/RdRP_Host.RData')
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


library(ggplot2)
library(Hmisc)
library(ggalluvial)
library(eoffice)
g = ggplot(df, aes(y= Freq, axis1 = Phylum, axis2 = Host)) +
  geom_alluvium(width= 1/12) +
  geom_stratum(width = 1/12, fill = c(phyla_cols[rev(levels(df$Phylum))],
                                      host_cols[rev(levels(df$Host))])) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size = 3) +
  scale_x_discrete(limits= c("Phylum", "Host"), expand = c(.05, .05)) +
  ylab("Number viruses in RCR95 set") + theme_bw()
topptx(g, file='../Urban_RNA_Virus_Figs/ICTV_Host.pptx', width = 6, height = 6)




# VIRON
df = RdRP_Host[RdRP_Host$Host == 'Vertebrates',c('Phylum','Virion')]
#sum(!is.na(df$virion_host))
df$Virion[df$Virion==''] = 'Others'
df$Virion[df$Virion=='non_human_mammalia'] = 'mammalia'
df$Virion[df$Virion=='non_mammalia_vertebrates'] = 'Others'
df$Virion = capitalize(df$Virion)
df$Virion = factor(df$Virion, levels = c('Others','Human','Mammalia'))
n = unique(df$Phylum)
n = n[order(n)]
g = ggplot(data = df, mapping = aes(x = Virion, fill = Phylum)) + 
  geom_bar() + scale_fill_manual(values = phyla_cols[n]) + scale_y_continuous (position = "right") +
  ylab("Number viruses in RCR95 set")+ xlab('') + theme_bw()
topptx(g, filename = '../Urban_RNA_Virus_Figs/Vertebrates_Host.pptx', width = 8, height = 6)



## Climate
library(plyr)
library(reshape2)
library(ggcorrplot)
load('../Urban_RNA_Virus_Data/meta_info.RData')
load('../Urban_RNA_Virus_Data/climate_info.RData')
load('../Urban_RNA_Virus_Data/RCR95_Phylum.RData')
load('../Urban_RNA_Virus_Data/RCR95_Abd.RData')

host = strsplit(RdRP_Host$Host, split = ';')
n = unlist(lapply(host, length))
df = data.frame(RdRP = rep(RdRP_Host$RdRP, n), Host = unlist(host),stringsAsFactors = F)

df = merge(df, RCR95_Abd, by.x = 'RdRP', by.y = 'RCR95')
df = ddply(df, .variables = 'Host', .fun = function(df){
  colSums(df[,-c(1:3)])
})
host = df$Host
df = cbind(data.frame('Sample_id' = colnames(df)[-1]), t(df[,-1]))
colnames(df)[-1] = host
df$Sample_id = gsub('^X','', gsub('\\.','-',df$Sample_id))


df = merge(meta_info[,c('Sample_id','City')], df)
df = merge(climate_info, df, by = 'City')
df = ddply(df[,-6], .variables = 'City', function(df){
  apply(df[,-1], 2, median)
})
climates = colnames(df)[2:5]
hosts = colnames(df)[6:ncol(df)]
hosts = hosts[hosts != 'Others']
cor_mat = matrix(NA, nrow = length(climates), ncol = length(hosts))
rownames(cor_mat) = climates
colnames(cor_mat) = hosts
pval_mat = cor_mat

for(i in climates){
  for (j in hosts) {
    r = cor.test(df[,i], df[,j])
    cor_mat[i,j] = r$estimate
    pval_mat[i,j] = r$p.value
  }
}
pval_mat = pval_mat[,-1]
cor_mat = cor_mat[,-1]

library(corrplot)
pdf('../Urban_RNA_Virus_Figs/Climate_Host_Cor.pdf', width = 8, height = 6)
corrplot(cor_mat, method = 'ellipse', p.mat = pval_mat,col = rev(RColorBrewer::brewer.pal(10,"RdYlBu")),
             insig = 'label_sig',sig.level=0.1)
dev.off()
