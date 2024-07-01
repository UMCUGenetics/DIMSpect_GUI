# functions for DIMSpect Shiny GUI
get_info_on_patient <- function(api, patient_id, sample_id = "", run_name = "", 
                                matrix_selected = "") {
  patient_table <- get_patient_info_db(api, patient_id, sample_id, run_name, matrix_selected)
  
  # if patient_table is empty, display an empty matrix
  if (is.null(patient_table)) {
    patient_table <- data.frame(matrix(ncol = 4, nrow = 0))
    patient_table <- as.data.frame(patient_table)
  }
  patient_table <- patient_table %>% rename("ID" = Sample.id, "Patient ID" = Sample.patient_id, 
                                            "Matrix" = Sample.type, "Run name" = name)
  patient_table <- patient_table %>% select(ID, `Patient ID`, `Run name`, Matrix)
  return(patient_table)
}

get_patient_info_db <- function(api, patient_id, sample_id, run, matrix) {
  patient_info <- api %>%
    req_url_path_append("patient/info") %>%
    req_url_query(pat_id = patient_id) %>%
    req_url_query(samp_id = sample_id) %>%
    req_url_query(run_name = run) %>%
    req_url_query(matrix_type = matrix) %>%
    req_perform() %>%
    resp_body_json(simplifyVector = T)
  
  patient_info <- as.data.frame.list(patient_info)
  return(patient_info)
}

get_run_zscore_data_db <- function(api, run_names, iden, zscore_min, zscore_max, samples) {
  run_zscore_data <- api %>%
    req_url_path_append("results/hmdb") %>%
    req_url_path_append(iden) %>%
    req_url_path_append("zscores") %>%
    req_url_path_append(zscore_min) %>%
    req_url_path_append(zscore_max) %>%
    req_url_query(run_names = run_names, .multi = "explode") %>%
    req_url_query(samples = samples, .multi = "explode") %>%
    req_perform() %>%
    resp_body_json(simplifyVector = T)
  
  run_zscore_data <- as.data.frame.list(run_zscore_data)
  run_zscore_data <- run_zscore_data %>% select(-c(DIMSResults.uuid, HMDB.uuid, DIMSResults.intensity))
  run_zscore_data$DIMSResults.polarity <- ifelse(run_zscore_data$DIMSResults.polarity == TRUE, "negative", "positive")
  run_zscore_data <- run_zscore_data %>% rename(Zscore = DIMSResults.z_score, 
                                                    scanmode = DIMSResults.polarity, m_z = DIMSResults.m_z, run_name = DIMSResults.run_name,
                                                    sample_id = DIMSResults.sample_id, HMDB_id = HMDB.hmdb_id, HMDB_name = HMDB.name,
                                                    HMDB_description = HMDB.description, HMDB_theor_mz = HMDB.theor_MZ)
  # TODO: verwijder distinct() na update DIMSdb
  run_zscore_data <- run_zscore_data %>% distinct()
  run_zscore_data <- run_zscore_data  %>% pivot_wider(names_from = sample_id, values_from = c(Zscore), names_glue = "{sample_id}_{.value}")
  run_zscore_data$HMDB_name <- gsub('"', '', run_zscore_data$HMDB_name)
  
  return(run_zscore_data)
}

get_patient_data <- function(api, patient_query_input, zscore_high, zscore_low, identified_only) {
  # get data for unique run names
  # print(patient_query_input)
  # TODO: check ID kolom
  all_run_names <- unique(patient_query_input$`Run name`)
  all_samples <- unique(patient_query_input$ID)
  
  run_data <- get_run_zscore_data_db(api, all_run_names, identified_only, zscore_low, zscore_high, all_samples)

  # sort data_filt_zscore on Z-score for query patient
  sample_cols <- paste0(all_samples, "_Zscore")
  run_data <- run_data %>% arrange(across(sample_cols, desc))
  return(run_data)
}

select_plot_data <- function(patient_data, selected_metabolites) {
  # patient_data <<- patient_data
  # info_columns <<- grep("HMDB_name|m_z", colnames(patient_data))
  # zscore_columns <<- grep("_Zscore", colnames(patient_data))
  # # selection of metabolites from table (https://yihui.shinyapps.io/DT-rows)
  # plot_data_selected <- patient_data[selected_metabolites, c(info_columns[2], zscore_columns)]
  
  plot_data_selected <- patient_data %>% select(HMDB_name, run_name, c(contains("_Zscore"))) %>% 
    slice(selected_metabolites)
  
  return(plot_data_selected)
}

create_violin_plots <- function(metab_data, sample_id) {
  colors_4plot <- c("#22E4AC", "#00B0F0", "#504FFF","#A704FD","#F36265","#DA0641")
  
  data_viool <- metab_data %>% select(-c(m_z, scanmode, HMDB_theor_mz, HMDB_description))
  
  data_viool <- data_viool %>% pivot_longer(cols = ends_with("_Zscore"), names_to = c("Sample"), values_to = "Z_score", values_drop_na = TRUE)
  data_viool$Sample <- gsub("_Zscore", "", data_viool$Sample)
  data_viool_pt_orig <- data_viool[which(data_viool$Sample %in% sample_id), ]

  data_viool$Z_score[data_viool$Z_score >  20] <-  20
  data_viool$Z_score[data_viool$Z_score <  -5] <-  -5
  
  data_viool_pt <- data_viool[which(data_viool$Sample %in% sample_id), ]
  data_viool_pt$Zscore_orig <- data_viool_pt_orig$Z_score
  
  data_viool <- data_viool[-which(data_viool$Sample %in% sample_id), ]
  
  violin_plot <- ggplot(data_viool, aes(x = Z_score, y = run_name)) + 
    xlim(-5, 20) + geom_violin(scale = "width") +
    theme(axis.text.y=element_text(size=rel(1.5)), plot.caption = element_text(size=rel(1)), 
          legend.position = "none", strip.text.y.left = element_text(angle = 0),
          strip.background = element_rect(colour = "black", fill = "white", linewidth = 1, linetype = "solid"),
          strip.text = element_text(size = 12)) +
    facet_grid(HMDB_name ~ ., scales="free", switch = "y", space = "free", labeller = as_labeller(split_label)) +
    geom_point(data = data_viool_pt, aes(fill = Zscore_orig), size = 5, shape=22) +
    geom_vline(xintercept = 2, col = "grey", lwd = 0.5, lty=2) +
    geom_vline(xintercept = -2, col = "grey", lwd = 0.5, lty=2) +
    scale_fill_gradientn(colors = colors_4plot, values = NULL, space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "colour") +
    geom_text(data = data_viool_pt, aes(16, label = paste0("Z = ", round(Zscore_orig, 2))), hjust = 0, vjust = +0.2, size = 5) +
    labs(x = "Z-scores", y = "Metabolites", color = "z-score") + 
    scale_y_discrete(position = "right")
  
  return(violin_plot)
}

split_label <- function(label) {
  str_replace_all(label, "(.{15})", "\\1\n")
}
