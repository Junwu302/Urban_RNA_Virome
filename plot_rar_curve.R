# Accumulation curve
# data used for plot is generated using the specaccum function in the vegan package
library(ggplot2)
library(ggthemes)
library(eoffice)

load('../Urban_RNA_Virus_Data/RvOTU_Accum.RData')
load('../Urban_RNA_Virus_Data/env_cols.RData')
df_all = data.frame(x = RvOTU_Accum$All$sites,
                    y = RvOTU_Accum$All$richness,
                    sd = RvOTU_Accum$All$sd, stringsAsFactors = F)
df_all$lower = df_all$y - df_all$sd
df_all$upper = df_all$y + df_all$sd
df_all$lower[df_all$lower < 0] = 0
g1 = ggplot(df_all, aes(x, y)) + geom_ribbon(aes(ymin = lower, ymax=upper, x = x), 
                                             fill = 'gray', alpha = 0.8) +
  geom_line(color ='black') + xlab('Number of samples')+ 
  ylab('Number of RvOTUs') + theme_bw()
topptx(g1, file = '../Urban_RNA_Virus_Figs/All_Accum.pptx', width = 8, height = 6)

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

g2 = ggplot(df_envs, aes(x, y)) + geom_ribbon(aes(ymin = lower, ymax=upper, x = x,fill = Env),alpha = 0.3) +
  geom_line(aes(color=Env),size = 1) + xlab('Number of samples')+ 
  ylab('Number of RvOTUs') + scale_color_manual(values = env_cols) + 
  scale_fill_manual(values = env_cols) +
  theme_bw()
topptx(g2, file = '../Urban_RNA_Virus_Figs/Env_Accum.pptx', width = 8, height = 6)

g2 = ggplot(df_envs, aes(x, y)) + geom_ribbon(aes(ymin = lower, ymax=upper, x = x,fill = Env),alpha = 0.3) +
  geom_line(aes(color=Env),size = 1) + facet_wrap( ~ Env, ncol=3,scales="free") +
  xlab('Number of samples')+ 
  ylab('Number of RvOTUs') + scale_color_manual(values = env_cols) + 
  scale_fill_manual(values = env_cols) +
  theme_bw()
topptx(g2, file = '../Urban_RNA_Virus_Figs/Env_Accum_grid.pptx', width = 8, height = 6)

df_envs = data.frame()
for(env in names(RvOTU_Accum)[-c(1,4,6,7)]){
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
g3 = ggplot(df_envs, aes(x, y)) + geom_ribbon(aes(ymin = lower, ymax=upper, x = x,fill = Env),alpha = 0.3) +
  geom_line(aes(color=Env),size = 1) + xlab('Number of samples')+ 
  ylab('Number of RvOTUs') + scale_color_manual(values = env_cols) + 
  scale_fill_manual(values = env_cols) +
  theme_bw()
topptx(g3, file = '../Urban_RNA_Virus_Figs/Env_Accum_sub.pptx', width = 8, height = 6)
