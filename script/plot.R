plot_function <- function(scen_r, destination="results/") {
  
  # ------------------------------------------------------------
  # CONSTANTS: Define aesthetics here (no overrides later)
  # ------------------------------------------------------------
  
  shape_vals <- c(
    OD=NA, OAS=3, LW=14, SS=20, OASD=25
  )
  
  color_vals <- c(
    OD="black",
    OAS="blue",
    LW="brown",
    SS="orange",
    OASD="red"
  )
  
  linetype_vals <- c(
    OD ="dotted",
    "TRUE"="solid",      # ours == TRUE  (OASD, OASB)
    "FALSE"="dashed"     # all other estimators
  )
  
  tex_labs <- c(
    LW=TeX(r'($italic(LW)$)'),
    OAS=TeX(r'($italic(OAS)$)'),
    OASD=TeX(r'($italic(OASD)$)'),
    SS=TeX(r'($italic(SS)$)')
  )
  
  keep_cols = c("OAS","LW","SS","OASD")
  
  base_theme <- theme(
    text = element_text(size=12),
    axis.title = element_text(size=12, face="italic"),
    axis.text = element_text(size=12),
    legend.title = element_text(size=10, face="italic"),
    legend.text = element_text(size=10),
    legend.position = "bottom",
    legend.box = "horizontal"
  )
  
  # ------------------------------------------------------------
  # Helper for processing data (shared by rho/PRIAL)
  # ------------------------------------------------------------
  
  process_data <- function(df) {
    df[setdiff(c("sd","r","n","OASB","OB"), scen_r)] <- NULL
    df %>%
      group_by(.data[[scen_r]]) %>%
      summarise(across(everything(), mean)) %>%
      ungroup()
  }
  
  # ------------------------------------------------------------
  # Generic plotting helper
  # ------------------------------------------------------------
  
  make_plot <- function(result, scen_r, yvar, ylab_tex) {
    
    result_long <- result %>%
      pivot_longer(!all_of(scen_r), names_to="method", values_to=yvar) %>%
      mutate(
        ours = method %in% c("OASD"),
        linetype_key = ifelse(method %in% c("OD"),
                              method,  # Oracle types handled explicitly
                              as.character(ours))
      )
    
    ggplot(result_long %>% filter(method != "s"),
           aes(x=.data[[scen_r]],
               y=.data[[yvar]],
               group=method,
               shape=method,
               color=method,
               linetype=linetype_key)) +
      
      geom_line(size=0.5) +
      geom_point(size=0.5) +
      
      scale_shape_manual(values=shape_vals,  breaks = keep_cols, labels=tex_labs) +
      scale_color_manual(values=color_vals,  breaks = keep_cols,  labels=tex_labs) +
      scale_linetype_manual(values=linetype_vals,  breaks = keep_cols) +
      
      xlab(scen_r) +
      ylab(TeX(ylab_tex)) +
      guides(linetype="none") +
      base_theme
  }
  
  # ------------------------------------------------------------
  # Panel 1: Average correlation (rho)
  # ------------------------------------------------------------
  df1 <- readRDS(paste0(destination, "coef_", scen_r, ".rds"))
  res1 <- process_data(df1)
  p1   <- make_plot(res1, scen_r, "rho", r'($italic(Average)$ $\rho$)')
  
  guide_plot <- cowplot::get_plot_component(p1, 'guide-box-bottom', return_all=TRUE)
  p1 <- p1 + theme(legend.position="none")
  
  # ------------------------------------------------------------
  # Panel 2: PRIAL
  # ------------------------------------------------------------
  df2 <- readRDS(paste0(destination, scen_r, ".rds"))
  res2 <- process_data(df2) %>%
    mutate(across(!all_of(c("s", scen_r)), ~ (s - .x)/s * 100))
  
  p2 <- make_plot(res2, scen_r, "PRIAL", r'($italic(PRIAL)$)') +
    theme(legend.position="none")
  
  # ------------------------------------------------------------
  # Combine panels
  # ------------------------------------------------------------
  combined_p <- (p2 | p1) / guide_plot +
    plot_layout(heights=c(1, 0.1))
  
  print(combined_p)
  
  ggsave(
    paste0(destination, scen_r, "_mse+coef.jpg"),
    width=6, height=4
  )
}



plot_corr_single_function <- function(scen, destination="results/") {
  
  # ------------------------------------------------------------
  # Extract scenario name (remove _corr or _inv)
  # ------------------------------------------------------------
  scen_r <- scen %>%
    str_replace("_corr", "") %>%
    str_replace("_inv", "")
  
  # ------------------------------------------------------------
  # EXACT SAME AESTHETICS AS MAIN FUNCTION
  # ------------------------------------------------------------
  
  shape_vals <- c(
    OD=NA, OAS=3, LW=14, SS=20, OASD=25
  )
  
  color_vals <- c(
    OD="black",
    OAS="blue",
    LW="brown",
    SS="orange",
    OASD="red"
  )
  
  linetype_vals <- c(
    OD ="dotted",
    "TRUE"  = "solid",
    "FALSE" = "dashed"
  )
  
  tex_labs <- c(
    LW  = TeX(r'($italic(LW)$)'),
    OAS = TeX(r'($italic(OAS)$)'),
    OASD = TeX(r'($italic(OASD)$)'),
    SS  = TeX(r'($italic(SS)$)')
  )
  
  keep_cols <- names(tex_labs)  # only these appear in legend
  
  base_theme <- theme(
    text = element_text(size=12),
    axis.title = element_text(size=12, face="italic"),
    axis.text = element_text(size=12),
    legend.title = element_text(size=10, face="italic"),
    legend.text = element_text(size=10),
    legend.position = "bottom",
    legend.box = "horizontal"
  )
  
  # ------------------------------------------------------------
  # PROCESS DATA (same logic as main function)
  # ------------------------------------------------------------
  
  df <- readRDS(paste0(destination, "/", scen, ".rds"))
  df[setdiff(c("sd","r","n","OASB","OB"), scen_r)] <- NULL
  
  result <- df %>%
    group_by(.data[[scen_r]]) %>%
    summarise(across(everything(), mean)) %>%
    ungroup() %>%
    mutate(across(
      !all_of(c("s", scen_r)),
      ~ (s - .x) / s * 100
    ))
  
  # ------------------------------------------------------------
  # Prepare long format & linetype mapping
  # ------------------------------------------------------------
  
  result_long <- result %>%
    pivot_longer(
      !all_of(scen_r),
      names_to = "method",
      values_to = "PRIAL"
    ) %>%
    mutate(
      ours = method %in% c("OASD"),
      linetype_key = ifelse(method %in% c("OD"),
                            method,
                            as.character(ours))
    )
  
  # ------------------------------------------------------------
  # FINAL PLOT (exactly same structure as main function)
  # ------------------------------------------------------------
  
  p <- ggplot(
    result_long %>% filter(method != "s"),
    aes(x = .data[[scen_r]],
        y = PRIAL,
        group = method,
        shape = method,
        color = method,
        linetype = linetype_key)
  ) +
    geom_line(size=0.5) +
    geom_point(size=0.5) +
    
    scale_shape_manual(
      values = shape_vals,
      breaks = keep_cols,
      labels = tex_labs
    ) +
    scale_color_manual(
      values = color_vals,
      breaks = keep_cols,
      labels = tex_labs
    ) +
    scale_linetype_manual(
      values = linetype_vals,
      breaks = c("TRUE","FALSE")   # OD/OB excluded
    ) +
    
    xlab(scen_r) +
    ylab(TeX(r'($italic(PRIAL)$)')) +
    guides(linetype="none") +
    base_theme
  
  return(p)
}

plot_inv_single_function <- function(scen, destination="results/") {
  
  # ------------------------------------------------------------
  # Extract scenario name (strip _corr or _inv)
  # ------------------------------------------------------------
  scen_r <- scen %>%
    str_replace("_corr", "") %>%
    str_replace("_inv", "")
  
  # ------------------------------------------------------------
  # MATCH EXACTLY THE MAIN FUNCTION'S AESTHETICS
  # ------------------------------------------------------------
  
  shape_vals <- c(
    OD=NA, OAS=3, LW=14, SS=20, OASD=25
  )
  
  color_vals <- c(
    OD="black",
    OAS="blue",
    LW="brown",
    SS="orange",
    OASD="red"
  )
  
  linetype_vals <- c(
    OD ="dotted",
    "TRUE"  = "solid",   # OASD, OASB
    "FALSE" = "dashed"   # OAS, LW, SS
  )
  
  tex_labs <- c(
    LW  = TeX(r'($italic(LW)$)'),
    OAS = TeX(r'($italic(OAS)$)'),
    OASD= TeX(r'($italic(OASD)$)'),
    SS  = TeX(r'($italic(SS)$)')
  )
  
  keep_cols <- names(tex_labs)   # OD, OB excluded
  
  base_theme <- theme(
    text = element_text(size=12),
    axis.title = element_text(size=12, face="italic"),
    axis.text = element_text(size=12),
    legend.title = element_text(size=10, face="italic"),
    legend.text = element_text(size=10),
    legend.position = "bottom",
    legend.box = "horizontal"
  )
  
  # ------------------------------------------------------------
  # Load and process data
  # ------------------------------------------------------------
  
  df <- readRDS(paste0(destination, "/", scen, ".rds"))
  df[setdiff(c("sd","r","n","OASB","OB"), scen_r)] <- NULL
  
  result <- df %>%
    group_by(.data[[scen_r]]) %>%
    summarise(across(everything(), mean)) %>%
    ungroup() %>%
    mutate(across(
      !all_of(c(scen_r, "MP")),
      ~ (MP - .x) / MP * 100
    ))
  
  # ------------------------------------------------------------
  # Long format + linetype mapping
  # ------------------------------------------------------------
  
  result_long <- result %>%
    pivot_longer(
      !all_of(scen_r),
      names_to = "method",
      values_to = "MSE"
    ) %>%
    mutate(
      ours = method %in% c("OASD"),
      linetype_key = ifelse(
        method %in% c("OD"),
        method, 
        as.character(ours)
      )
    )
  
  # ------------------------------------------------------------
  # Final plot (exact style of main function)
  # ------------------------------------------------------------
  
  p <- ggplot(
    result_long %>% filter(!method %in% c("s","MP")),
    aes(
      x = .data[[scen_r]],
      y = MSE,
      group = method,
      shape = method,
      color = method,
      linetype = linetype_key
    )
  ) +
    geom_line(size=0.5) +
    geom_point(size=0.5) +
    
    scale_shape_manual(values=shape_vals,
                       breaks=keep_cols,
                       labels=tex_labs) +
    scale_color_manual(values=color_vals,
                       breaks=keep_cols,
                       labels=tex_labs) +
    scale_linetype_manual(values=linetype_vals,
                          breaks=c("TRUE","FALSE")) +
    
    xlab(scen_r) +
    ylab(TeX(r'($italic(PRIAL)_{italic(INV)}$)')) +
    guides(linetype = "none") +
    base_theme
  
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

