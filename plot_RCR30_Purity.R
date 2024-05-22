library(plyr)
library(reshape2)
library(ggplot2)
library(eoffice)
load('../Urban_RNA_Virus_Data/RCR30.RData')
load('../Urban_RNA_Virus_Data/RCR95.RData')
RCR30 = RCR30[RCR30$RCR95 %in% RCR95$RCR95,]
num = table(RCR30$RCR30)
df = ddply(RCR30[RCR30$RCR30 %in% names(num)[num >= 20],], .variables = 'RCR30', .fun = function(x){
  x = data.frame(t(apply(x[,6:9],2,function(a){
    n_total = length(a)
    n_unknown = sum(is.na(a))
    a = a[!is.na(a)]
    n_max = table(a)
    n_max = as.numeric(n_max[order(n_max, decreasing = T)][1])
    purity = n_max/(n_total - n_unknown)
    return(c(n_total, n_unknown/n_total, n_max, purity))
  })))
  colnames(x) = c("n_total", "n_unknown", "n_max", "purity")
  x$tax = rownames(x)
  return(x)
})
df$tax = factor(df$tax, levels = c('Phylum','Class','Order','Family'))

g = ggplot(df, aes(x=purity, y = 100*(1-n_unknown), size = n_total)) +
  geom_point(alpha = .5, color='#1685a9', shape = 16)+ ylab('Annotation percentage %') + 
  facet_grid(.~tax) + theme_bw()
topptx(g, file = '../Urban_RNA_Virus_Figs/RCR30_Purity.pptx',
       width = 8, height = 4)
