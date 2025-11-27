library(ggpubr)
df = readRDS(glue('100r=0.3.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p1 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.3'))

df = readRDS(glue('100r=0.5.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p2 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.5'))

df = readRDS(glue('100r=0.7.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p3 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.7'))

df = readRDS(glue('100r=0.9.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p4 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.9'))

df = readRDS(glue('r=0.5.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(-sample, -constant, -diagonal, -lw, -rblw, -test, -oas_c, -oracle_c, -corr_oas, -corr_oracle) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p2 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.5'))

df = readRDS(glue('r=0.7.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(-sample, -constant, -diagonal, -lw, -rblw, -test, -oas_c, -oracle_c, -corr_oas, -corr_oracle) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p3 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.7'))

df = readRDS(glue('r=0.9.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(-sample, -constant, -diagonal, -lw, -rblw, -test, -oas_c, -oracle_c, -corr_oas, -corr_oracle) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p4 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.9'))

ggarrange(p1, p2, p3, p4, nrow=2, ncol = 2, common.legend = TRUE, legend="bottom")

# comparison 1
df = readRDS(glue('r=0.3.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(lw, rblw, oas_c, oas_v, oracle_c, oracle_v,sample_size) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p1 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.3'))

df = readRDS(glue('r=0.5.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(lw, rblw, oas_c, oas_v, oracle_c, oracle_v,sample_size) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p2 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.5'))

df = readRDS(glue('r=0.7.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(lw, rblw, oas_c, oas_v, oracle_c, oracle_v,sample_size) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p3 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.7'))

df = readRDS(glue('r=0.9.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(lw, rblw, oas_c, oas_v, oracle_c, oracle_v,sample_size) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p4 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.9'))

ggarrange(p1, p2, p3, p4, nrow=2, ncol = 2, common.legend = TRUE, legend="bottom")

# comparison 2
df = readRDS(glue('r=0.3.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(oas_v, corr_oracle_c, oracle_v,corr_oas_c, sample_size) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p1 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.3'))

df = readRDS(glue('r=0.5.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(oas_v, corr_oracle_c, oracle_v,corr_oas_c, sample_size) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p2 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.5'))

df = readRDS(glue('r=0.7.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(oas_v, corr_oracle_c, oracle_v,corr_oas_c, sample_size) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p3 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.7'))

df = readRDS(glue('r=0.9.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(oas_v, corr_oracle_c, oracle_v,corr_oas_c, sample_size) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p4 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.9'))

ggarrange(p1, p2, p3, p4, nrow=2, ncol = 2, common.legend = TRUE, legend="bottom")

# comparison 3
df = readRDS(glue('r=0.3.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(oas_v, shafer, sample_size) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p1 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.3'))

df = readRDS(glue('r=0.5.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(oas_v, shafer, sample_size) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p2 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.5'))

df = readRDS(glue('r=0.7.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(oas_v, shafer, sample_size) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p3 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.7'))

df = readRDS(glue('r=0.9.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(oas_v, shafer, sample_size) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p4 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.9'))

ggarrange(p1, p2, p3, p4, nrow=2, ncol = 2, common.legend = TRUE, legend="bottom")

# comparison 4
df = readRDS(glue('r=0.3.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(oas_v_unknown_mean, shafer, sample_size) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p1 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.3'))

df = readRDS(glue('r=0.5.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(oas_v_unknown_mean, shafer, sample_size) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p2 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.5'))

df = readRDS(glue('r=0.7.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(oas_v_unknown_mean, shafer, sample_size) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p3 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.7'))

df = readRDS(glue('r=0.9.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  select(oas_v_unknown_mean, shafer, sample_size) %>% 
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p4 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.9'))

plot = ggarrange(p1, p2, p3, p4, nrow=2, ncol = 2, common.legend = TRUE, legend="bottom")
annotate_figure(plot, top = text_grob("sdlog=0.1",face = "bold", size = 14))
