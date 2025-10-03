# functions for DIMSpect Shiny GUI
get_info_on_patient <- function(api, patient_id, sample_id = "", run_name = "",
                                matrix_selected = "") {
  patient_table <- get_patient_info_db(api, patient_id, sample_id, run_name, matrix_selected)

  # if patient_table is empty, display an empty matrix
  if (is.null(patient_table)) {
    patient_table <- data.frame(matrix(ncol = 4, nrow = 0))
    patient_table <- as.data.frame(patient_table)
  }
  patient_table <- patient_table %>% rename(
    "ID" = id, "Patient ID" = patient_id,
    "Matrix" = matrix, "Run name" = name
  )
  patient_table <- patient_table %>% select(ID, `Patient ID`, `Run name`, Matrix)
  patient_table$Matrix <- str_to_title(patient_table$Matrix)
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
  run_zscore_data$DIMSResults.polarity <- ifelse(run_zscore_data$DIMSResults.polarity == TRUE, "positive", "negative")
  run_zscore_data <- run_zscore_data %>% rename(
    Zscore = DIMSResults.z_score,
    scanmode = DIMSResults.polarity, m_z = DIMSResults.m_z, run_name = DIMSResults.run_name,
    sample_id = DIMSResults.sample_id, HMDB_id = HMDB.hmdb_id, HMDB_name = HMDB.name,
    HMDB_description = HMDB.description, HMDB_theor_mz = HMDB.theor_mz,
    ppm = DIMSResults.ppm_dev, chem_formula = HMDB.chem_formula, sec_hmdb_id = HMDB.sec_hmdb_id
  )

  # TODO: verwijder distinct() na update DIMSdb
  run_zscore_data <- run_zscore_data %>% distinct()
  run_zscore_data <- run_zscore_data %>% pivot_wider(names_from = sample_id, values_from = c(Zscore), names_glue = "{sample_id}_{.value}")
  run_zscore_data$HMDB_name <- gsub('"', "", run_zscore_data$HMDB_name)

  run_zscore_data$HMDB_id <- hmdb_id_to_hmdb_html(run_zscore_data$HMDB_id)

  return(run_zscore_data)
}

get_patient_data <- function(api, patient_query_input, zscore_high, zscore_low, identified_only) {
  # get data for unique run names
  # TODO: check ID kolom
  all_run_names <- unique(patient_query_input$`Run name`)
  all_samples <- unique(patient_query_input$ID)

  run_data <- get_run_zscore_data_db(api, all_run_names, identified_only, zscore_low, zscore_high, all_samples)

  adduct_pos <- list(
    "0" = "[M+H]+", "1" = "[M+Na]+", "2" = "[M+K]+", "3" = "[M+NaCl]+", "4" = "[M+NH4]+", "5" = "[M+2Na-H]+",
    "6" = "[M+CH3OH]+", "7" = "[M+KCl]+", "8" = "[M+NaK-H]+"
  )
  adduct_neg <- list(
    "0" = "[M-H]-", "1" = "[M+Cl]-", "2" = "[M+For]-", "3" = "[M+NaCl]-", "4" = "[M+KCl-", "5" = "[M+H2PO4]-",
    "6" = "[M+HSO4]-", "7" = "[M+Na-H]-", "8" = "[M+K-H]-", "9" = "[M-H2O]-", "10" = "[M-2H]-",
    "11" = "[M+I]-", "12" = "[M+Ac]-"
  )

  run_data <- run_data %>%
    select(-c(DIMSResults.row_hash, HMDB.hmdb_key, HMDB_description)) %>%
    relocate(c(HMDB_name, HMDB_id, chem_formula, adduct, HMDB_theor_mz, m_z, ppm), .before = everything()) %>%
    relocate(c(run_name, sec_hmdb_id, scanmode), .after = last_col()) %>%
    mutate(adduct = if_else(
      scanmode == "positive",
      unname(adduct_pos[as.character(adduct)]),
      unname(adduct_neg[as.character(adduct)])
    ))

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

  plot_data_selected <- patient_data %>%
    select(HMDB_name, run_name, c(contains("_Zscore"))) %>%
    slice(selected_metabolites)

  return(plot_data_selected)
}

format_data_violin_plot <- function(metab_data, sample_id) {
  patients <- sample_id
  data_viool <- metab_data %>% select(-c(m_z, scanmode, HMDB_theor_mz, HMDB_id, chem_formula, ppm, sec_hmdb_id, adduct))

  data_viool_pivot <- reshape2::melt(data_viool, id.vars = c("HMDB_name", "run_name"))
  data_viool_pivot <- data_viool_pivot %>% drop_na()
  colnames(data_viool_pivot) <- c("HMDB_name", "run_name", "Sample", "Zscore")
  data_viool_pivot <- data_viool_pivot %>%
    mutate(
      Sample = gsub("_Zscore", "", Sample),
      Type = ifelse(substr(Sample, 1, 1) == "P", "Patient", "Control"),
      Highlight = ifelse(Sample %in% patients, "Target Patient", "Other")
    )

  violin_data <- data_viool_pivot %>% filter(!Sample %in% patients)
  patient_data <- data_viool_pivot %>% filter(Sample %in% patients)
  control_data <- data_viool_pivot %>% filter(Type == "Control")

  data_violin_plot <<- list(
    violin_data = violin_data,
    control_data = control_data,
    patient_data = patient_data
  )

  return(data_violin_plot)
}

create_violin_plots <- function(list_data_violin) {
  violin_data <- list_data_violin$violin_data
  control_data <- list_data_violin$control_data
  patient_data <- list_data_violin$patient_data

  # Plot is flipped, by coord_flip(), so x = y
  ggplot(violin_data, aes(x = run_name, y = Zscore)) +
    geom_violin(aes(fill = run_name), color = "black", alpha = 0.8, trim = TRUE) +
    geom_jitter(
      aes(text = paste0(
        "<b>Patient:</b> ", Sample,
        "<br>",
        "<b>Metabolite:</b> ", HMDB_name,
        "<br>",
        "<b>Z-score:</b> ", round(Zscore, 3),
        "<br>",
        "<b>Run:</b> ", run_name
      )),
      height = 0, width = 0.1, color = "black"
    ) +
    geom_point(
      data = patient_data,
      aes(
        y = Zscore,
        color = Sample,
        text = paste0(
          "<b>Patient:</b> ", Sample,
          "<br>",
          "<b>Metabolite:</b> ", HMDB_name,
          "<br>",
          "<b>Z-score:</b> ", round(Zscore, 3),
          "<br>",
          "<b>Run:</b> ", run_name
        )
      ),
      position = position_jitter(width = 0.1, height = 0.1, seed = 1),
      size = 2.5, shape = 21, fill = "white", stroke = 1, show.legend = FALSE
    ) +
    coord_cartesian(xlim = c(min(violin_data$Zscore), 21)) +
    theme_bw() +
    labs(x = "Z-scores", color = "z-score") +
    facet_grid(HMDB_name ~ ., scales = "free", space = "free", labeller = as_labeller(split_label)) +
    theme(
      axis.text.y = element_text(size = rel(1.5)), plot.caption = element_text(size = rel(1)),
      legend.position = "none", strip.text.y.left = element_text(angle = 0),
      strip.background = element_rect(colour = "black", fill = "white", linewidth = 1, linetype = "solid"),
      strip.text = element_text(size = 12), axis.title.y = element_blank()
    ) +
    geom_hline(yintercept = 2, col = "grey", lwd = 0.5, lty = 2) +
    geom_hline(yintercept = -1.5, col = "grey", lwd = 0.5, lty = 2) +
    scale_y_discrete(position = "right") +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Set3") +
    coord_flip() # flip x and y axis
}

create_box_plots <- function(metab_data, sample_ids) {
  ggplot(metab_data, aes(
    x = Sample, y = Zscore,
    text = paste0(
      "<b>Patient:</b> ", Sample,
      "<br>",
      "<b>Metabolite:</b> ", HMDB_name,
      "<br>",
      "<b>Z-score:</b> ", round(Zscore, 3),
      "<br>",
      "<b>Run:</b> ", run_name
    )
  )) +
    geom_bar(aes(fill = Sample %in% sample_ids), stat = "identity") +
    scale_fill_manual(guide = "none", breaks = c(FALSE, TRUE), values = c("royalblue4", "red4")) +
    labs(x = "", y = "Z-score") +
    facet_grid(HMDB_name ~ run_name, scales = "free", switch = "y", labeller = labeller(HMDB_name = as_labeller(split_label))) +
    theme_bw() +
    geom_hline(yintercept = 2, col = "grey", lwd = 0.5, lty = 2) +
    geom_hline(yintercept = -1.5, col = "grey", lwd = 0.5, lty = 2) +
    theme(
      legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
      strip.background = element_rect(colour = "black", fill = "white", linewidth = 1, linetype = "solid"),
      strip.text = element_text(size = 12), axis.text = element_text(size = 12)
    ) +
    scale_y_continuous(position = "right")
}

format_data_single_pat_plot <- function(metab_data, metab_rows, sample_ids) {
  metab_data_samp <- metab_data %>% select(-c(HMDB_theor_mz, m_z, ppm, scanmode, HMDB_id, sec_hmdb_id))
  colnames(metab_data_samp) <- str_remove(colnames(metab_data_samp), "_Zscore")

  metabs <- metab_data_samp %>% slice(metab_rows)
  metab_data_samp <- metab_data_samp %>% filter(HMDB_name %in% metabs$HMDB_name & adduct %in% metabs$adduct)

  metab_data_long <- reshape2::melt(metab_data_samp, id.vars = c("HMDB_name", "run_name", "chem_formula", "adduct"))
  metab_data_long <- metab_data_long %>% drop_na()
  metab_data_long <- metab_data_long %>%
    rename(Sample = variable, Zscore = value) %>%
    mutate_at(c("adduct", "Sample"), as.character) %>%
    mutate(
      type = ifelse(grepl("^C", Sample), "Control", "Patient"),
      Sample = ifelse(type == "Control", "Controls", Sample),
      Zscore = as.numeric(Zscore)
    ) %>%
    drop_na() %>%
    group_by(Sample, run_name) %>%
    mutate(group_size = n()) %>%
    arrange(Sample)
  return(metab_data_long)
}

create_single_pat_plot <- function(metab_data) {
  ggplot(metab_data, aes(
    x = Sample, y = Zscore,
    text = paste0(
      "<b>Patient:</b> ", Sample,
      "<br>",
      "<b>Metabolite:</b> ", HMDB_name,
      "<br>",
      "<b>Z-score:</b> ", round(Zscore, 3),
      "<br>",
      "<b>Run:</b> ", run_name
    )
  )) +
    geom_boxplot(data = subset(metab_data, group_size > 2), aes(fill = type, colour = type)) +
    geom_point(data = subset(metab_data, group_size <= 2), aes(colour = type, fill = type)) +
    scale_fill_manual(values = c("Control" = "green", "Patient" = "#b20000")) +
    scale_color_manual(values = c("Control" = "black", "Patient" = "#b20000")) +
    labs(x = "", y = "Z-score") +
    facet_grid(HMDB_name ~ run_name, scales = "free", switch = "y", labeller = labeller(HMDB_name = as_labeller(split_label))) +
    theme_bw() +
    theme(
      legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
      strip.background = element_rect(colour = "black", fill = "white", linewidth = 1, linetype = "solid"),
      strip.text = element_text(size = 12), axis.text = element_text(size = 12)
    ) +
    geom_hline(yintercept = 2, col = "grey", lwd = 0.5, lty = 2) +
    geom_hline(yintercept = -1.5, col = "grey", lwd = 0.5, lty = 2) +
    scale_y_continuous(position = "right")
}

format_data_single_pat_table <- function(metab_data, sample_ids) {
  sample_zscores <- paste0(sample_ids, "_Zscore")
  metab_data_pat <- metab_data %>% select(HMDB_id, HMDB_name, chem_formula, adduct, HMDB_theor_mz, m_z, ppm, all_of(sample_zscores), sec_hmdb_id, scanmode, run_name)
  
  metab_data_pat_wide <- metab_data_pat %>%
    pivot_wider(
      id_cols = c(HMDB_id, HMDB_name, chem_formula, adduct, HMDB_theor_mz, m_z, ppm, sec_hmdb_id, scanmode),
      names_from = run_name,
      values_from = contains("Zscore"),
      names_glue = "{.value}-{run_name}"
    )

  metab_data_pat_wide <- metab_data_pat_wide[, !sapply(metab_data_pat_wide, function(x) all(is.na(x)))]
  
  zscore_cols <- grep("_Zscore-", colnames(metab_data_pat_wide), value = TRUE)
  metab_data_pat_wide$avg_value <- rowMeans(metab_data_pat_wide[, zscore_cols], na.rm = TRUE)
  
  ungrouped_cols <- c("HMDB_name", "HMDB_id", "chem_formula", "adduct", "HMDB_theor_mz", "m_z", "ppm", "sec_hmdb_id", "scanmode", "avg_value")
  col_order <- c(ungrouped_cols, colnames(metab_data_pat_wide)[grep("_Zscore-", colnames(metab_data_pat_wide))])
  metab_data_pat_wide <- metab_data_pat_wide[, col_order]
  return(metab_data_pat_wide)
}

make_layout_single_pat_tab <- function(cols_sample_patient_table) {
  ungrouped_cols <- c("Name", "HMDB ID", "Molecular formula", "Adduct", "Theoretical mass", "M/Z value",
                    "ppm", "Secondary HMDB IDs", "Scanmode", "Average Z-score")
  
  # Groepeer de sample-run kolommen
  cols_sample_run <- cols_sample_patient_table[!cols_sample_patient_table %in% ungrouped_cols]
  cols_sample_run_split <- strsplit(cols_sample_run, "-")
  
  sample_runs_list <- list()
  for (sample_runs in cols_sample_run_split) {
    sample <- sub("_Zscore$", "", sample_runs[[1]])
    sample_runs_list[[sample]] <- c(sample_runs_list[[sample]], sample_runs[[2]])
  }
  
  # Eerste rij: ongroepeerde kolommen (met rowspan = 2) + sample namen (met colspan)
  header_row_1 <- c(
    lapply(ungrouped_cols, function(col) tags$th(rowspan = 2, col)),
    lapply(names(sample_runs_list), function(sample) {
      tags$th(colspan = length(sample_runs_list[[sample]]), sample)
    })
  )
  
  # Tweede rij: run-namen
  header_row_2 <- lapply(unlist(sample_runs_list, use.names = FALSE), tags$th)
  
  # Bouw de sketch
  sketch <- withTags(
    table(
      class = "display",
      thead(
        tr(header_row_1),
        tr(header_row_2)
      )
    )
  )
  
  return(sketch)
}


split_label <- function(label) {
  str_replace_all(label, "(.{15})", "\\1\n")
}

get_metab_info <- function(api, hmdb_id, hmdb_name) {
  if (nchar(hmdb_id) > 1) {
    metab_info <- api %>%
      req_url_path_append("/hmdb/info/hmdb_id") %>%
      req_url_path_append(hmdb_id) %>%
      req_perform() %>%
      resp_body_json(simplifyVector = T)
  } else {
    metab_info <- api %>%
      req_url_path_append("/hmdb/info/hmdb_name") %>%
      req_url_path_append(hmdb_name) %>%
      req_perform() %>%
      resp_body_json(simplifyVector = T)
  }
  metab_info <- as.data.frame.list(metab_info)

  metab_info$hmdb_id <- hmdb_id_to_hmdb_html(metab_info$hmdb_id)
  metab_info <- metab_info %>% select(uuid, hmdb_id, name, chem_formula, theor_mz, sec_hmdb_id, description)
  return(metab_info)
}

is_valid_hmdb_id <- function(input_hmdb_id) {
  grepl("^HMDB\\d{5}(\\d{2})?$", input_hmdb_id)
}

hmdb_id_to_hmdb_html <- function(hmdb_id) {
  paste0('<a href=\"https://hmdb.ca/metabolites/', hmdb_id, '\" target=_\"blank\">', hmdb_id, "</a>")
}

change_colnames_patient_table <- function(colnames_old) {
  sample_names <- colnames_old[grepl("_Zscore", colnames_old)]
  sample_names <- gsub("_Zscore", " Zscore", sample_names)
  colnames_new <- c(
    "Name", " HMDB ID", "Molecular formula", "Adduct", "Theoretical mass", "M/Z value",
    "ppm", sample_names, "Run name", "Secondary HMDB IDs", "Scanmode"
  )
  return(colnames_new)
}

change_colnames_single_pat_table <- function(colnames_old) {
  sample_names <- colnames_old[grepl("_Zscore", colnames_old)]
  colnames_new <- c("Name", "HMDB ID", "Molecular formula", "Adduct", "Theoretical mass", "M/Z value",
                    "ppm", "Secondary HMDB IDs", "Scanmode", "Average Z-score", sample_names)
  return(colnames_new)
}

get_metab_zscores <- function(api, metab_info_df, zscore_high, zscore_low) {
  metab_zscores <- api %>%
    req_url_path_append("results/metab") %>%
    req_url_path_append(zscore_low) %>%
    req_url_path_append(zscore_high) %>%
    req_url_query(hmdb_uuids = metab_info_df$uuid, .multi = "explode") %>%
    req_perform() %>%
    resp_body_json(simplifyVector = T)

  metab_zscores <- as.data.frame.list(metab_zscores)

  metab_zscores_info <- left_join(metab_zscores, metab_info_df, by = join_by(DIMSResultsHMDBLink.hmdb_id == uuid))

  metab_zscores_info <- metab_zscores_info %>%
    select(
      hmdb_id, name, DIMSResultsHMDBLink.adduct, DIMSResults.sample_id, DIMSResults.z_score, DIMSResults.run_name,
      DIMSResults.polarity, DIMSResults.m_z
    ) %>%
    rename(
      zscore = DIMSResults.z_score,
      scanmode = DIMSResults.polarity, m_z = DIMSResults.m_z, run_name = DIMSResults.run_name,
      sample_id = DIMSResults.sample_id, adduct = DIMSResultsHMDBLink.adduct
    ) %>%
    mutate_if(is.numeric, round, 5)

  metab_zscores_info$scanmode <- ifelse(metab_zscores_info$scanmode == TRUE, "positive", "negative")
  metab_zscores_info <- change_num_to_adduct(metab_zscores_info)

  return(metab_zscores_info)
}

change_num_to_adduct <- function(df) {
  adduct_pos <- list(
    "0" = "[M+H]+", "1" = "[M+Na]+", "2" = "[M+K]+", "3" = "[M+NaCl]+", "4" = "[M+NH4]+", "5" = "[M+2Na-H]+",
    "6" = "[M+CH3OH]+", "7" = "[M+KCl]+", "8" = "[M+NaK-H]+"
  )
  adduct_neg <- list(
    "0" = "[M-H]-", "1" = "[M+Cl]-", "2" = "[M+For]-", "3" = "[M+NaCl]-", "4" = "[M+KCl-", "5" = "[M+H2PO4]-",
    "6" = "[M+HSO4]-", "7" = "[M+Na-H]-", "8" = "[M+K-H]-", "9" = "[M-H2O]-", "10" = "[M-2H]-",
    "11" = "[M+I]-", "12" = "[M+Ac]-"
  )

  df <- df %>%
    mutate(adduct = if_else(
      scanmode == "positive",
      unname(adduct_pos[as.character(adduct)]),
      unname(adduct_neg[as.character(adduct)])
    ))
  return(df)
}

get_metab_plot <- function(metab_data) {
  ggplot(metab_data, aes(
    x = sample_id, y = zscore,
    text = paste0(
      "<b>Patient:</b> ", sample_id,
      "<br>",
      "<b>Metabolite:</b> ", name,
      "<br>",
      "<b>Adduct:</b> ", adduct,
      "<br>",
      "<b>Z-score:</b> ", round(zscore, 3),
      "<br>",
      "<b>Run:</b> ", run_name
    )
  )) +
    geom_point() +
    theme_bw() +
    theme(
      legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
      strip.background = element_rect(colour = "black", fill = "white", linewidth = 1, linetype = "solid"),
      strip.text = element_text(size = 12)
    )
}
