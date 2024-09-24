library(plyr)
load('../datap/hyla_cols')

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

# Figure A
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



## Figure 5B
library(plyr)
library(reshape2)
library(ggcorrplot)
load('../data/meta_info.RData')
load('../data/climate_info.RData')
load('../data/RCR95_Phylum.RData')
load('../data/RCR95_abd.RData')

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
corrplot(cor_mat, method = 'ellipse', p.mat = pval_mat,col = rev(RColorBrewer::brewer.pal(10,"RdYlBu")),
         insig = 'label_sig',sig.level=0.1)


# Figure 5C
# WDI vs viruse association
library(plyr)
library(reshape2)
library(pheatmap)
library(ggplotify)
library(ggplot2)
library(eoffice)
library(RColorBrewer)

load('..data/WDI_Host_Cor.RData')
load('../data/WDI_Meta.RData')
for(host in names(WDI_Cor)){
  df = WDI_Cor[[host]]
  df = df[df$Sig !='NoSig',]
  fl = paste0(c('../Supplementary_Files/',host,'_WDI.tsv'),collapse = '')
  write.table(df, file = fl, sep ='\t', row.names = F, quote = F)
}

df = data.frame()
for(host in names(WDI_Cor)){
  tmp = WDI_Cor[[host]]
  tmp = tmp[tmp$Sig != 'NoSig',]
  tmp$Host = host
  df = rbind(df, tmp)
}
df = merge(df, WDI_Meta[,c(1,5)], by.x = 'Series.Code', by.y = 'Code')
df$Topic = gsub(':.*$','',df$Topic) 
m = table(df$Topic)
df$Topic[df$Topic%in% names(m)[m <30]] = 'Others'

mat = matrix(NA, nrow = length(unique(df$Series.Code)), ncol = length(unique(df$Host)))
colnames(mat) = unique(df$Host)
rownames(mat) = unique(df$Series.Code)
for(host in names(WDI_Cor)){
  x = WDI_Cor[[host]]$Cor
  names(x) = WDI_Cor[[host]]$Series.Code
  mat[rownames(mat), host] = x[rownames(mat)]
}
row_ind = hclust(dist(mat))$order
col_ind = hclust(dist(t(mat)))$order

mat = matrix(NA, nrow = length(unique(df$Series.Code)), ncol = length(unique(df$Host)))
colnames(mat) = unique(df$Host)
rownames(mat) = unique(df$Series.Code)
for(i in 1:nrow(df)){
  host = df$Host[i]
  wdi = df$Series.Code[i]
  mat[wdi, host] = df$Cor[i]
}
annotation_row = data.frame(Topic = df$Topic[!duplicated(df$Series.Code)],
                            stringsAsFactors = F)
annotation_row$Topic = factor(annotation_row$Topic)
rownames(annotation_row) = df$Series.Code[!duplicated(df$Series.Code)]
ann_colors = c(brewer.pal(9,"Set1")[c(1:4)],'lightgray')
names(ann_colors) = c('Economic Policy & Debt','Environment','Health',
                      'Social Protection & Labor','Others')
ann_colors = list(Topic=ann_colors)

g= as.ggplot(pheatmap(mat[row_ind, col_ind], show_rownames = F, annotation_row = annotation_row,
                      annotation_colors = ann_colors,cluster_rows = F, 
                      cluster_cols = F, na_col = '#F8F8F8'))

# Figure 5D
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



#### Figure 5E
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


# Figure 5F
library(plyr)
load('../data/RCR95_Phylum.RData')
load('../data/UHGGspacer_match.RData')
load('../data/UHGGspacers.RData')
RCR95_Phylum$RCR95 = gsub('\\.[0-9]*$','',gsub('_[0-9]*$','',RCR95_Phylum$RCR95))
RCR95_Phylum = RCR95_Phylum[!duplicated(RCR95_Phylum$RCR95),]

uhgg_spacer_num = table(UHGGspacers$genome)
uhgg_spacer_num = uhgg_spacer_num[uhgg_spacer_num > 2]
UHGGspacers = UHGGspacers[UHGGspacers$genome %in% names(uhgg_spacer_num),]
UHGGspacer_match = UHGGspacer_match[UHGGspacer_match$mismatch <= 1 & 
                                      UHGGspacer_match$length > 20 &
                                      UHGGspacer_match$gapopen <= 1,]
UHGGspacer_match$sacc = gsub('\\.[0-9]$','',UHGGspacer_match$sacc)

UHGGspacer_match = merge(RCR95_Phylum[!duplicated(RCR95_Phylum$Contig),c('RCR95','Phylum')], UHGGspacer_match[,1:2], 
                         by.x = 'RCR95', by.y = 'sacc')
UHGGspacer_match = merge(UHGGspacer_match,UHGGspacers[,c(2,3,4)], by.x = 'qacc',by.y = 'spacer_id')

colnames(UHGGspacer_match) = c('Spacer_id','RCR95','Virus_Phylum','Species_rep','Lineage')
UHGGspacer_match$Virus_Phylum[is.na(UHGGspacer_match$Virus_Phylum)] = 'Unknown'

UHGGspacer_match$Spacer_Kingdom = unlist(lapply(strsplit(UHGGspacer_match$Lineage, split = ';'), function(x){
  x = gsub('_[A-Z]$','',gsub('^d__','',grep('^d__', x, value = T)))
  if(length(x) == 0){x = ''}
  return(x)
}))
UHGGspacer_match$Spacer_Phylum = unlist(lapply(strsplit(UHGGspacer_match$Lineage, split = ';'), function(x){
  x = gsub('_[A-Z]$','',gsub('^p__','',grep('^p__', x, value = T)))
  if(length(x) == 0){x = ''}
  return(x)
}))
UHGGspacer_match$Spacer_Class = unlist(lapply(strsplit(UHGGspacer_match$Lineage, split = ';'), function(x){
  x = gsub('_[A-Z]$','',gsub('^c__','',grep('^c__', x, value = T)))
  if(length(x) == 0){x = ''}
  return(x)
}))
UHGGspacer_match$Spacer_Genus = unlist(lapply(strsplit(UHGGspacer_match$Lineage, split = ';'), function(x){
  x = gsub('_[A-Z]$','',gsub('^g__','',grep('^g__', x, value = T)))
  if(length(x) == 0){x = ''}
  return(x)
}))
UHGGspacer_match$Spacer_Species = unlist(lapply(strsplit(UHGGspacer_match$Lineage, split = ';'), function(x){
  x = gsub('_[A-Z]$','',gsub('^s__','',grep('^s__', x, value = T)))
  if(length(x) == 0){x = ''}
  return(x)
}))

df = df[,c(2,10,3,8)]
colnames(df)[1:2] = c('from','to')
df = data.frame(table(UHGGspacer_match[!duplicated(UHGGspacer_match[,c('RCR95','Spacer_Species')]),
                                       c('Virus_Phylum','Spacer_Class')]))
df$Virus_Phylum = as.character(df$Virus_Phylum)
df$Spacer_Class = as.character(df$Spacer_Class)
df = df[df$Freq > 0,]
n1 = ddply(df, .variables = 'Spacer_Class', .fun = function(x){
  sum(x$Freq)
})
n1 = n1[order(n1$V1, decreasing = T),]
n1 = rbind(n1[n1$V1 >= 30,], data.frame(Spacer_Class='Others', V1 = sum(n1$V1[n1$V1 < 30])))
df$Spacer_Class[!(df$Spacer_Class %in% n1$Spacer_Class)] = 'Others'

n2 = ddply(df, .variables = 'Virus_Phylum', .fun = function(x){
  sum(x$Freq)
})
n2 = n2[order(n2$V1, decreasing = T),]
n2 = n2[c(1:4,6:8,5),]
df = ddply(df, .variables = c('Virus_Phylum','Spacer_Class'), .fun = function(x){sum(x$Freq)})
df$Spacer_Class = factor(df$Spacer_Class, levels = n1$Spacer_Class)
df$Virus_Phylum = factor(df$Virus_Phylum, levels = n2$Virus_Phylum)
colnames(df)[2] = 'Class'

g = ggplot(data = df, mapping = aes(x = Virus_Phylum, y = V1, fill = Class)) +
  geom_bar(stat ='identity') + scale_fill_tableau(palette = 'Tableau 10') +
  theme_bw() + theme(axis.text.x = element_text(angle = 90,hjust = 1)) +
  xlab('') + ylab('Number of species')









