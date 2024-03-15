# functions for DIMSpect Shiny GUI
get_info_on_patient <- function(patient_id = "", sample_id = "", run_name = "", 
                                matrix_selected = "") {
  # patient <- api %>%
  #     req_url_path_append('patients') %>%
  #     req_url_path_append(patient_id) %>%
  #     req_perform() %>%
  #     resp_body_json()
  # start selection either from patient_id or from sample_id
  if (patient_id != "" | sample_id !="") {
    patient_table <- sample_info[which(grepl(patient_id, sample_info[ , "sample_name"]) & 
                                       grepl(sample_id,  sample_info[ , "sample_name"])), ]
  }
  # filter on run name
  if (run_name != "" & !is.null(patient_table)) {
    patient_table <- patient_table[grep(run_name, patient_table[ , "run_name"]), ]
  }
  # filter for matrix
  if (!is.null(patient_table)) {
    patient_table <- patient_table[grep(matrix_selected, patient_table[ , "matrix"]), ]
  }
  # if patient_table is empty, display an empty matrix
  if (is.null(patient_table)) {
    patient_table <- sample_info[0, ]
  }
  patient_df <- as.data.frame(patient_table)
  return(patient_df)
}

# query all patients
get_all_patients <- function(patient_ids) {
  all_patients <- api %>%
    req_url_path_append('patients') %>%
    req_perform() %>%
    resp_body_json(simplifyVector=TRUE)
  return(all_patients)
}

get_patient_data <- function(patient_query_input, zscore_high = 2, zscore_low = -1.5, 
                             identified_only = "Identified only") {
  # get data for unique run names
  all_run_names <- unique(patient_query_input$run_name)
  all_samples <- unique(patient_query_input$sample_name)
  # get data from (pseudo)database
  for (run_number in 1:length(all_run_names)) {
    query_run_name <- all_run_names[run_number]
    data_selected_run <- pseudo_database[[query_run_name]]
    # TODO: multiple runs into a list data_all_runs. Get this from real database.
  }
  # TODO: check if data_selected_run is not empty
  
  # filter data for identified or unidentified peak groups
  # TODO: this filter doesn't work yet.
  empty_rows <- which(data_selected_run$assi_HMDB == "")
  if (identified_only == "Identified only" && length(empty_rows) > 0) {
    data_selected_run <- data_selected_run[-empty_rows, ]
  } else if (identified_only == "Unidentified only" && length(empty_rows) > 0) {
    data_selected_run <- data_selected_run[empty_rows, ]
  }
  
  # reduce data to metabolites with aberrant Z-scores for query sample_ids.
  for (query_sample_id in all_samples) {
    zscore_col <- paste0(query_sample_id, "_Zscore")
    zscore_colnr <- which(colnames(data_selected_run) == zscore_col)
    select_rows <- data_selected_run[ , zscore_colnr] > zscore_high | data_selected_run[ , zscore_colnr] < zscore_low
    data_filt_zscore <- data_selected_run[select_rows, ]
    # TODO: if there's more than one sample, store more than one data_filt_zscore in a list
  }
  select_columns <- c(which(colnames(data_filt_zscore) == "HMDB_code"), 
                      which(colnames(data_filt_zscore) == "mzmed.pgrp"),
                      which(colnames(data_filt_zscore) == "assi_HMDB"),
                      which(colnames(data_filt_zscore) == "scanmode"),
                      grep("Zscore", colnames(data_filt_zscore)))
  data_filt_zscore <- data_filt_zscore[ , select_columns]
  # change column name for "assi_HMDB" and "mzmed.pgrp"
  colnames(data_filt_zscore)[which(colnames(data_filt_zscore) == "assi_HMDB")] <- "HMDB_name"
  colnames(data_filt_zscore)[which(colnames(data_filt_zscore) == "mzmed.pgrp")] <- "m_z"
  
  # sort data_filt_zscore on Z-score for query patient
  zscore_col <- paste0(query_sample_id, "_Zscore")
  zscore_colnr <- which(colnames(data_filt_zscore) == zscore_col)
  sort_order <- sort(data_filt_zscore[ , zscore_colnr], index.return=TRUE)
  data_sorted <- data_filt_zscore[sort_order$ix, ]
  
  # Many entries contain several HMDB IDs. Use only first one in first column
  data_sorted$full_HMDB_code <- data_sorted$HMDB_code
  data_sorted$full_HMDB_name <- data_sorted$HMDB_name
  for (row_index in 1:nrow(data_sorted)) {
    data_sorted$HMDB_code[row_index] <- strsplit(data_sorted$HMDB_code[row_index], ";")[[1]][1]
    data_sorted$HMDB_name[row_index] <- strsplit(data_sorted$HMDB_name[row_index], ";")[[1]][1]
  }
  
  return(data_sorted)
}

# Create violin plots (adapted from code from DIMS pipeline)
create_violin_plots <- function(sample_id, data_perrun) {
  
  # set parameters for plots
  plot_height <- 9.6 
  plot_width <- 6
  fontsize <- 1 
  circlesize <- 0.8 
  colors_4plot <- c("#22E4AC", "#00B0F0", "#504FFF","#A704FD","#F36265","#DA0641")
  #                   green     blue      blue/purple purple    orange    red
  
  # page headers:
  # page_headers <- names(metab_perpage)
  
  # create a violin plot of all metabolites in data_perrun
  zscore_col <- paste0(sample_id, "_Zscore")
  zscore_colnr <- which(colnames(data_perrun) == zscore_col)
  pt_list_2plot <- data_perrun[ , zscore_colnr]
  data_sorted_min1column <- data_sorted[ , -zscore_colnr]
  # shorten entries in column HMDB_code
  for (row_index in 1:nrow(data_sorted_min1column)) {
    data_sorted_min1column$HMDB_code[row_index] <- strsplit(data_sorted_min1column$HMDB_code[row_index][[1]], ";")[[1]][1]
  }
  # put data in long format. Something goes wrong here.
  metab_list_2plot <- reshape2::melt(data_sorted_min1column, id.vars = "HMDB_code")
  srt <- sort(metab_list_2plot$HMDB_code, index.return=TRUE)
  metab_list_2plot <- metab_list_2plot[srt$ix, ]
  metab_list_2plot$value[metab_list_2plot$value >  20] <-  20
  metab_list_2plot$value[metab_list_2plot$value <  -5] <-  -5
  # for (row_number in 1:nrow(data_sorted)) {
  # extract original data for patient of interest (pt_name) before cut-offs
  # pt_list_2plot_orig <- data_sorted[ , zscore_colnr]
  # cut off Z-scores higher than 20 or lower than -5 (for nicer plots)
  #metab_list_2plot$value[metab_list_2plot$value >  20] <-  20
  #metab_list_2plot$value[metab_list_2plot$value <  -5] <-  -5
  # extract data for patient of interest (pt_name)
  #pt_list_2plot <- data_sorted[ , zscore_colnr]
  # restore original Z-score before cut-off, for showing Z-scores in PDF
  # pt_list_2plot$value_orig <- pt_list_2plot_orig$value
  # remove patient of interest (pt_name) from list; violins will be made up of controls and other patients
  # data_sorted_1column <- data_sorted[row_number, -zscore_colnr]
  # put in long format for ggplot
  # metab_list_2plot <- reshape2::melt(data_sorted_1column, id.vars = "HMDB_code")
  
  # draw violin plot. This is code for 20 violin plots; modify.
  ggplot_object <- ggplot(metab_list_2plot, aes(x=value, y=HMDB_code)) +
    theme(axis.text.y=element_text(size=rel(fontsize)), plot.caption = element_text(size=rel(fontsize))) +
    # xlim(-5, 20) +
    geom_violin(scale="width") +
    geom_point(data = pt_list_2plot, aes(color=value), size = 3.5*circlesize, shape=22, fill="white") +
    scale_fill_gradientn(colors = colors_4plot, values = NULL, space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "colour") +
    # add Z-score value for patient of interest at x=16
    geom_text(data = pt_list_2plot, aes(16, label = paste0("Z=", round(value_orig, 2))), hjust = "left", vjust = +0.2, size = Z_size) +
    # add labels. Use font Courier to get all the plots in the same location.
    labs(x = "Z-scores", y = "Metabolites", subtitle = sub_perpage, color = "z-score") + 
    theme(axis.text.y = element_text(family = "Courier", size=6)) +
    # do not show legend
    theme(legend.position="none") +
    # add title 
    ggtitle(label = paste0("Results for patient ", pt_name)) + 
    # labs(x = "Z-scores", y = "Metabolites", title = paste0("Results for patient ", pt_name), subtitle = sub_perpage, color = "z-score") + 
    # add vertical lines
    geom_vline(xintercept = 2, col = "grey", lwd = 0.5, lty=2) +
    geom_vline(xintercept = -2, col = "grey", lwd = 0.5, lty=2)
  
  suppressWarnings(print(ggplot_object))
  
} # end create_violin_plots

