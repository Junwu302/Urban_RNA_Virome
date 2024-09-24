library(ggplot2)
library(ggsci)
library(eoffice)

# Figure 2A
load('../data/RdRP_blast.RData')
load('../data/RCR95.RData')
th = c(0,seq(50, 100, by=10))
Res = data.frame()
for(i in 2:length(th)){
  df_blast = RdRP_blast[RdRP_blast$identity <th[i] & 
                          RdRP_blast$identity>=th[i-1],c(1,5:8)]
  df_geNomad = RCR95[RCR95$RdRP %in% df_blast$id,c(3,6:9)]
  df = merge(df_blast, df_geNomad, by.x = 'id', by.y = 'RdRP')
  df[is.na(df)] = ''
  frac = data.frame(Phylum = NA, Class= NA, Order=NA, Family=NA)
  for(j in 2:5){
    x = df[,j]
    y = df[,j+4]
    ind = !(is.na(y))
    frac[1,j-1] = sum(x[ind] == y[ind])/sum(ind)
  }
  Res = rbind(Res, frac)
}
Res$th = th[-1]
Res = melt(Res, id.vars = 'th')
g = ggplot(Res, mapping = aes(x = th, y = value, fill = variable)) +
  geom_bar(stat = 'identity', position = 'dodge2') +
  xlab('Identity') + ylab('Consistency ')+
  scale_fill_npg() + theme_bw()

