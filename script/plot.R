plot_function = function(scen, destination="../results/"){
  scen_r = str_replace(scen, "_corr", "")
  scen_r = str_replace(scen_r, "_inv", "")
  
  res_sample_size = readRDS(paste0(destination,'coef_',scen, '.rds'))
  res_sample_size[setdiff(c("sd","r","n"), scen_r)] = NULL
  result = res_sample_size %>% group_by(.data[[scen_r]]) %>% summarise(across(everything(), ~mean(.)))
  result = pivot_longer(result,cols=!all_of(scen_r), names_to='method',values_to ='rho') %>% mutate(ours = ifelse(method %in% c('OASD','OASB'),TRUE,FALSE))
  p1 = ggplot(result %>% filter(!method %in% c('s')), aes(x=.data[[scen_r]],y=rho,group=method,color=method,shape=method, linetype=ours))+
    geom_line(data = result %>% filter(!method %in% c('s','Oracle'))) + 
    geom_line(data=result %>% filter(method %in% c('Oracle','Oracle_g')),color='black',show.legend=FALSE)+
    geom_point(data = result %>% filter(!method %in% c('s','Oracle')), size=1) +
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"))+
    scale_shape_manual(values = c("OD"=10,"OAS"=3,"LW"=14,"SS"=20,"OASD"=25, "OASB"=2),
                       labels=c("LW"=TeX(r'($italic(LW)$)'), "OAS"=TeX(r'($italic(OAS)$)'), "OASD"=TeX(r'($italic(OASD)$)'), "OD"=TeX(r'($italic(OD)$)'), "SS"=TeX(r'($italic(SS)$)'), "OASB"=TeX(r'($italic(OASB)$)')))+
    scale_color_manual(values = c("Oracle" = "black","OAS" = "blue","OD"="orange","LW"="brown","SS"="purple","OASD"="red", "OASB"="green"),
                       labels=c("LW"=TeX(r'($italic(LW)$)'), "OAS"=TeX(r'($italic(OAS)$)'), "OASD"=TeX(r'($italic(OASD)$)'), "OD"=TeX(r'($italic(OD)$)'), "SS"=TeX(r'($italic(SS)$)'), "OASB"=TeX(r'($italic(OASB)$)')))+ 
    xlab(scen) + ylab(TeX(r'($italic(Average)$ $\rho$)'))+ guides(linetype = "none") + 
    theme(
    text = element_text(size = 12),              # Global text size
    axis.title = element_text(size = 12,face = "italic"),         # Axis titles
    axis.text = element_text(size = 12),          # Axis text
    plot.title = element_text(size = 12),         # Plot title
    legend.title = element_text(size = 10, face="italic"),       # Legend title
    legend.text = element_text(size = 10),
    legend.position = "bottom", 
    legend.box = "horizontal", 
    legend.direction = "horizontal"
  )
  guide_plot <- cowplot::get_plot_component(p1, 'guide-box-bottom', return_all = TRUE)
  p1 = p1 + theme(legend.position = 'none')
  res_sample_size = readRDS(paste0(destination,scen,'.rds'))
  res_sample_size[setdiff(c("sd","r","n"), scen_r)] = NULL
  result = res_sample_size %>% group_by(.data[[scen_r]]) %>% summarise(across(everything(), ~mean(.))) %>% ungroup() %>% mutate(across(!all_of(c("s",scen_r)), ~(.data$s-.)/.data$s)*100)
  result = pivot_longer(result,cols=!all_of(scen_r),names_to='method',values_to ='PRIAL') %>% mutate(ours = ifelse(method %in% c('OASD','OASB'),TRUE,FALSE))
  p2 = ggplot(result %>% filter(!method %in% c('s')), aes(x=.data[[scen_r]],y=PRIAL,group=method,color=method, shape=method, linetype=ours))+
    geom_line(data = result %>% filter(!method %in% c('s','Oracle'))) + 
    geom_line(data=result %>% filter(method %in% c('Oracle','Oracle_g')),color='black',show.legend=FALSE)+
    geom_point(data = result %>% filter(!method %in% c('s','Oracle')), size=1)+
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"))+
    scale_shape_manual(values = c("OD"=10,"OAS"=3,"LW"=14,"SS"=20,"OASD"=25, "OASB"=2),
                       labels=c("LW"=TeX(r'($italic(LW)$)'), "OAS"=TeX(r'($italic(OAS)$)'), "OASD"=TeX(r'($italic(OASD)$)'), "OD"=TeX(r'($italic(OD)$)'), "SS"=TeX(r'($italic(SS)$)'), "OASB"=TeX(r'($italic(OASB)$)')))+
    scale_color_manual(values = c("Oracle" = "black","OAS" = "blue","OD"="orange","LW"="brown","SS"="purple","OASD"="red", "OASB"="green"),
                       labels=c("LW"=TeX(r'($italic(LW)$)'), "OAS"=TeX(r'($italic(OAS)$)'), "OASD"=TeX(r'($italic(OASD)$)'), "OD"=TeX(r'($italic(OD)$)'), "SS"=TeX(r'($italic(SS)$)'), "OASB"=TeX(r'($italic(OASB)$)')))+ 
    xlab(scen) + ylab(TeX(r'($italic(PRIAL)$)')) + guides(linetype = "none") + 
    theme(
    text = element_text(size = 12),              # Global text size
    axis.title = element_text(size = 12,face = "italic"),         # Axis titles
    axis.text = element_text(size = 12),          # Axis text
    plot.title = element_text(size = 12),         # Plot title
    legend.title = element_text(size = 12),       # Legend title
    legend.text = element_text(size = 12),         # Legend text
    legend.position = "none"
  )
  
  combined_p = (p2 | p1)/guide_plot + 
    plot_layout(heights = c(1, 0.1))
  combined_p
  ggsave(paste0(destination,scen,'_mse+coef.jpg'), width = 6, height = 4)
}

plot_corr_single_function=function(scen){
  scen_r = str_replace(scen, "_corr", "")
  scen_r = str_replace(scen_r, "_inv", "")
  
  res_sample_size = readRDS(paste0("results/",scen,'.rds'))
  res_sample_size[setdiff(c("sd","r","n"), scen_r)] = NULL
  result = res_sample_size %>% group_by(.data[[scen_r]]) %>% summarise(across(everything(), ~mean(.))) %>% ungroup() %>% mutate(across(!all_of(c("s",scen_r)), ~(.data$s-.)/.data$s)*100)
  result = pivot_longer(result,cols=!all_of(scen_r),names_to='method',values_to ='PRIAL')%>% mutate(ours = ifelse(method %in% c('OASD'),TRUE,FALSE))
  p = ggplot(result %>% filter(!method %in% c('s','OD')), aes(x=.data[[scen_r]],y=PRIAL,group=method,color=method,shape=method,linetype=ours))+
    geom_line(data = result %>% filter(!method %in% c('s','Oracle','OD'))) + 
    geom_line(data=result %>% filter(method %in% c('Oracle','Oracle_g')),color='black',show.legend=FALSE)+
    geom_point(data = result %>% filter(!method %in% c('s','Oracle','OD')), size=1)+
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"))+
    scale_shape_manual(values = c("OD"=10,"OAS"=3,"LW"=14,"SS"=20,"OASD"=25, "OASB"=2),
                       labels=c("LW"=TeX(r'($italic(LW)_{italic(corr)}$)'), "OAS"=TeX(r'($italic(OAS)_{italic(corr)}$)'), "OASD"=TeX(r'($italic(OASD)$)'), "OD"=TeX(r'($italic(OD)$)'), "SS"=TeX(r'($italic(SS)$)'), "OASB"=TeX(r'($italic(OASB)$)')))+
    scale_color_manual(values = c("Oracle" = "black","OAS" = "blue","OD"="orange","LW"="brown","SS"="purple","OASD"="red", "OASB"="green"),
                       labels=c("LW"=TeX(r'($italic(LW)_{italic(corr)}$)'), "OAS"=TeX(r'($italic(OAS)_{italic(corr)}$)'), "OASD"=TeX(r'($italic(OASD)$)'), "OD"=TeX(r'($italic(OD)$)'), "SS"=TeX(r'($italic(SS)$)'), "OASB"=TeX(r'($italic(OASB)$)')))+ 
    xlab(scen_r) + 
    ylab(TeX(r'($italic(PRIAL)$)')) + guides(linetype = "none") + 
    theme(
    text = element_text(size = 12),              # Global text size
    axis.title = element_text(size = 12,face = "italic"),         # Axis titles
    axis.text = element_text(size = 12),          # Axis text
    plot.title = element_text(size = 12),         # Plot title
    legend.title = element_text(size = 12, face='italic'),       # Legend title
    legend.text = element_text(size = 12, face='italic'),         # Legend text
    legend.position = "bottom"
  )
  return(p)
}

plot_inv_single_function=function(scen){
  scen_r = str_replace(scen, "_corr", "")
  scen_r = str_replace(scen_r, "_inv", "")
  
  res_sample_size = readRDS(paste0("results/",scen,'.rds'))
  res_sample_size[setdiff(c("sd","r","n"), scen_r)] = NULL
  result = res_sample_size %>% group_by(.data[[scen_r]]) %>% summarise(across(everything(), ~mean(.))) %>% ungroup() %>% mutate(across(!all_of(c("s",scen_r,"MP")), ~(.data$MP-.)/.data$MP)*100)
  result = pivot_longer(result,cols=!all_of(scen_r),names_to='method',values_to ='MSE') %>% mutate(ours = ifelse(method %in% c('OASD'),TRUE,FALSE))
  p = ggplot(result %>% filter(!method %in% c('s', 'MP')), aes(x=.data[[scen_r]],y=MSE,group=method,color=method,shape=method, linetype=ours))+
    geom_line(data = result %>% filter(!method %in% c('s','RBLW','OD','Oracle','MP'))) + 
    geom_line(data=result %>% filter(method %in% c('Oracle','Oracle_g')),color='black',show.legend=FALSE)+
    geom_point(data = result %>% filter(!method %in% c('s','RBLW','OD','Oracle','MP')), size=1) +
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"))+
    scale_shape_manual(values = c("OD"=10,"OAS"=3,"LW"=14,"SS"=20,"OASD"=25, "OASB"=2),
                       labels=c("LW"=TeX(r'($italic(LW)$)'), "OAS"=TeX(r'($italic(OAS)$)'), "OASD"=TeX(r'($italic(OASD)$)'), "OD"=TeX(r'($italic(OD)$)'), "SS"=TeX(r'($italic(SS)$)'), "OASB"=TeX(r'($italic(OASB)$)')))+
    scale_color_manual(values = c("Oracle" = "black","OAS" = "blue","OD"="orange","LW"="brown","SS"="purple","OASD"="red", "OASB"="green"),
                       labels=c("LW"=TeX(r'($italic(LW)$)'), "OAS"=TeX(r'($italic(OAS)$)'), "OASD"=TeX(r'($italic(OASD)$)'), "OD"=TeX(r'($italic(OD)$)'), "SS"=TeX(r'($italic(SS)$)'), "OASB"=TeX(r'($italic(OASB)$)')))+ 
    xlab(scen_r) + ylab(TeX(r'($italic(PRIAL)_{italic(INV)}$)')) + guides(linetype = "none") + 
    theme(
    text = element_text(size = 12),              # Global text size
    axis.title = element_text(size = 12,face = "italic"),         # Axis titles
    axis.text = element_text(size = 12),          # Axis text
    plot.title = element_text(size = 12),         # Plot title
    legend.title = element_text(size = 12, face='italic'),       # Legend title
    legend.text = element_text(size = 12, face='italic'),         # Legend text
    legend.position = "bottom"
  )
  return(p)
}

plot_agg_function=function(feature){
  plot_single_function = ifelse(feature=="corr", plot_corr_single_function, plot_inv_single_function)
  p1=plot_single_function(paste0("sd_",feature))
  guide_plot <- cowplot::get_plot_component(p1, 'guide-box-bottom', return_all = TRUE)
  p1 = p1 + theme(legend.position = 'none')
  p2=plot_single_function(paste0("n_",feature))
  p2 = p2 + theme(legend.position = 'none',axis.title.y = element_blank())
  p3=plot_single_function(paste0("r_",feature))
  p3 = p3 + theme(legend.position = 'none',axis.title.y = element_blank())
  
  combined_p = (p1 | p2 | p3)/guide_plot + 
    plot_layout(heights = c(1, 0.1))
  combined_p
  ggsave(paste0("results/",feature,'_tot.jpg'), width = 6, height = 4)
}

