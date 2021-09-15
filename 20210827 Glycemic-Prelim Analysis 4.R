library(tidyverse)

### ==== Phase 0. Reading in the Data (Always run this FIRST) ====
    qt <- read.delim('20210817_quant_table_glycemic.txt', as.is = T, check.names = F)
    qtc <- cov_dratio_filter(qt, num_info_col = 7, cov_row = 1, dratio_row = 3)[[1]]
    bad_biomarkers <- cov_dratio_filter(qt, num_info_col = 7, cov_row = 1, dratio_row = 3)[[2]]
    
### Data-Specific Cleaning    
    qtc <- qtc %>% remove_pool()
    qtc$`Time in Minutes` <- paste0(qtc$`Time in Minutes`, ' Min')
    qtc$Glucose <- as.numeric(qtc$Glucose)
    qtc <- qtc[-which(qtc$Subjects == 'Subject_14'),] ## Incomplete Subject Data

### ==== Phase 1. QC Metric Scripts: Column Column CoV and D-Ratio Improvements ====
    ### Pre & Post effects of column correction by MRMkit
    # CoV & D-Ratio Comparison
    
    plot_cov_dratio_comparisons('quant_table_pre_correction.txt', '20210817_quant_table_glycemic.txt',
                                cov_row = 1, dratio_row = 3,
                                pre_num_info_col = 4, post_num_info_col = 7)
    file.rename('CoV & D-Ratio Comparison Graphs', '1. CoV & D-Ratio Comparison Graphs')
    
    # MWUT between different columns
    p1 <- plot_columnDiff_mwutVal('quant_table_pre_correction.txt', 4, study_codename = 'Glycemic', sort_colname = 'Column')
    p2 <- plot_columnDiff_mwutVal('20210817_quant_table_glycemic.txt', 7, study_codename = 'Glycemic', sort_colname = 'Column')
    
    output_plot(p1, foldername = '1. Column-Correction Plots', filename = 'Pre-Corrected C1 vs C2', filetype = 'png')
    output_plot(p2, foldername = '1. Column-Correction Plots', filename = 'Post-Corrected C1 vs C2', filetype = 'jpeg')

    rm(p1, p2)
### ==== Phase 2.1. RM ANOVA for all 4 RCT Groups (& Comparison Contrasts) ====
    qtc_log <- qtc %>% log_df(7)
    qtc2 <- qtc_log %>% outliers_begone(7)
    qtc2_berry <- qtc2 %>% group_by(Classification) %>% group_split()
    
    total_RM_master <- data.frame(mtbl = colnames(qtc2)[8:ncol(qtc2)])
    for (i in 1:length(qtc2_berry)){
      my_df <- qtc2_berry[[i]]
      type <- unique(my_df$Classification)
      print(type)
      
      num_info_col <- 7
      
      # Factor Casting
      my_df <- my_df %>% within({
        `Subjects` <- as.factor(`Subjects`)
        Classification <- as.factor(Classification)
        `Time in Minutes` <- as.factor(`Time in Minutes`)
      })
      
      RManova_master <- data.frame(matrix(nrow = 0, ncol = 2)) ### <<< Colnum changes according to no. of factors & interactions
      names(RManova_master) <- c('Metabolites', paste0(type,'_Effect_Time'))
      
      for (j in (num_info_col+1):ncol(my_df)){
        metabolite <- paste0('`', colnames(my_df)[j], '`')
        
        input_data <- cbind(my_df[,c(1:num_info_col)], my_df[,j])
        if(any(is.na(input_data[,8]))){input_data <- input_data[-which(is.na(input_data[,8])),]}
      
      
      test_anova <- summary(aov(eval(parse(text = metabolite)) ~ `Time in Minutes`+ ### Multiply this by other factors to see interaction and effects of other factors alone 
                                  Error(`Subjects`/(`Time in Minutes`)), ### Which means divide this by the factors as well
                                data = input_data))
      
      pval_list <- str_remove_all(metabolite, '`')
      for (k in 2){ ## Change the range accordingly <<<<
        pval <- test_anova[[k]][[1]]$`Pr(>F)`[1]  ### We only look at the effect caused by the interaction between PID & Time (Where Time is nested in PID - Time changes with the same PID)
        pval_list <- append(pval_list, pval)
        }
      RManova_master[nrow(RManova_master)+1,] <- pval_list
      
      rm(k, j, pval, pval_list, test_anova, metabolite)
      }
      
      total_RM_master <- total_RM_master %>% cbind(RManova_master[,2])
      names(total_RM_master)[ncol(total_RM_master)] <- paste0(type,'_Effect_Time')
      
      rm(type, RManova_master, my_df, i, input_data)
    }
    
    write.csv(total_RM_master, '2.1. Glycemic-RM_ANOVA-Pvalues.csv', row.names = F)
    
    rm(qtc2, qtc2_berry)
    
### ==== Phase 2.1.1. Filtering for the significant values for contrasts ====
    ### We only want those with significant changes in time across 1 group of contrast, and not other.
    
    for (i in 2:5){
      total_RM_master[,i] <- as.numeric(unlist(total_RM_master[,i]))
      rm(i)
    }
    RM_contrasts <- list(c('P-E', 'Q-E'), c('P-F', 'Q-F'), c('P-E', 'P-F'), c('Q-E', 'Q-F'))
    for(contrast in RM_contrasts){
      first_colnum <- which(grepl(contrast[1], colnames(total_RM_master)))
      second_colnum <- which(grepl(contrast[2], colnames(total_RM_master)))
      my_df <- total_RM_master[,c(1, first_colnum, second_colnum)]
      my_df[,c(2,3)] <- my_df[,c(2,3)] %>% apply(MARGIN = 2, FUN = as.numeric)
      
      my_df <- my_df[which(my_df[,2] < 0.05 | my_df[,3] < 0.05),] # Selecting only the significant values
      my_df <- my_df[-which(my_df[,2] < 0.05 & my_df[,3] < 0.05),] # Selecting significant values for only EITHER category (Change in 1, no change in another)
      my_df <- my_df[order(my_df[,2]),]
      write.csv(my_df, paste0('2.1.1. ', contrast[1], ' vs ', contrast[2], ' - Pvalues.csv'), row.names = F)
      
      rm(first_colnum, second_colnum, my_df, contrast)
    }
    
    rm(total_RM_master, RM_contrasts)
    
### ==== Phase 2.2. Plotting the graphs for all 4 RCT Groups ====
    dual_boxplot_contrasts <- list(c('P', 'Q'), c('E', 'F'), c('P-E', 'P-F'), c('Q-E', 'Q-F'),
                                   c('P-E', 'Q-E'), c('P-F', 'Q-F'))
    
    for (contrast in dual_boxplot_contrasts){
      label <- paste0(contrast[1], ' vs ', contrast[2])
      
      qtc1 <- qtc %>% filter(grepl(contrast[1], qtc$Classification)) %>% add_column(Classification2 = contrast[1], .before = 'Classification')
      qtc2 <- qtc %>% filter(grepl(contrast[2], qtc$Classification)) %>% add_column(Classification2 = contrast[2], .before = 'Classification')
      
      qtc_plotting <- qtc1 %>% rbind(qtc2)
      batch_png_plot(qtc_plotting, 'Time in Minutes', 8, toLog10 = T, remove_outliers = 'quantile',
                 plot_type = 'multi_boxplot', colName_toFill = 'Classification2', legend_position = 'right')
      
      file.rename('plotsOutputDir', paste0('2.2.1. Plots_', label))

      rm(label, qtc1, qtc2, qtc_plotting, contrast)
    }
    
    ### For whole of qtc with 4 boxplots each at each timepoint:
    plot_graph(qtc, 'Time in Minutes', 123, toLog10 = T, remove_outliers = 'quantile',
               plot_type = 'multi_boxplot', colName_toFill = 'Classification_ID', legend_position = 'bottom')
    batch_png_plot(qtc, 'Time in Minutes', 6, toLog10 = T, remove_outliers = 'quantile',
                   plot_type = 'multi_boxplot', colName_toFill = 'Classification_ID', legend_position = 'right')
    file.rename('plotsOutputDir', '2.2.2. Multi Boxplots')   
    
    
### ==== Phase 2.3. Generation of P-values between each timepoint for P vs Q ====
    num_info_col <- 6
    qtc2 <- qtc3 <-  qtc_log %>% outliers_begone(7)
    # qtc2 <- qtc2 %>% add_column(PairSorting = paste(qtc2$Classification_ID, qtc2$`Subjects`, sep = '_'), .after = 'Classification_ID')
    qtc2_berry <- qtc2 %>% group_by(`Time in Minutes`) %>% group_split()
    
    mtbl <- colnames(qtc2)[(num_info_col+1):ncol(qtc2)]
    mycontrasts <- list(c('P', 'Q'), c('E', 'F'), c('P-E', 'P-F'), c('Q-E', 'Q-F'),
                        c('P-E', 'Q-E'), c('P-F', 'Q-F'))
    
    container <- data.frame(mtbls = mtbl)
    for (i in 1:length(qtc2_berry)){
      my_df <- qtc2_berry[[i]]
      timepoint <- unique(my_df$`Time in Minutes`)
      
      my_df <- my_df %>% add_column(SortPairs = paste0(my_df$`Subjects`, '_', my_df$Classification_ID), .after = 'Classification_ID')
      pvals_df <- a_stat_test(my_df, 'Classification', 7, mycontrasts = mycontrasts, test_type = 't.test',
                              paired = T, colName_SortPairs = 'SortPairs')
      container <- container %>% cbind(pvals_df[2:(length(mycontrasts) + 1)])
      names(container)[(ncol(container) - 5):ncol(container)] <- paste0(timepoint, '_',
                                                                        names(container)[(ncol(container) - 5):ncol(container)])
      rm(timepoint, i)
    }
    names(container)[1] <- 'Metabolites'
    
### Sorting the metabolites
    ordering <- ordering2 <- 1
    for (i in 2:7){ordering <- ordering %>% append(seq(i, ncol(container), 6))}
    container <- container[,ordering]
    for (i in seq(2, ncol(container), 5)){ordering2 <- ordering2 %>% append(c(i, i+3, i+4, i+1, i+2))}
    container <- container[,ordering2]
    rm(ordering, ordering2)
    
    write.csv(container, '2.3. GIN_T-test_Pvalues_AllContrasts(At Each Timepoint).csv', row.names = F)

### ==== Phase 2.4. Calculating iAUC and AUC (IMPT) ====
    rm(qt, bad_biomarkers) # Who needs these?
    
    auc_berry <- qtc %>% group_by(`Subjects`, Classification_ID) %>% group_split()
    
    auc_df <- qtc[1,-c(2, 4:6)] ## Just a container.
    auc_df <- auc_df %>% add_column('AUC_Type' = NA, .after = 'Classification')
    
    master_df_ref_list <- list()
    
    ## Calculate iAUC, AUC and net iAUC (reason being that biomarkers can have lower levels than BL)
    for (i in 1:length(auc_berry)){ ## For every Patient, every RCT...
      this_patient_rct <- auc_berry[[i]]
      unique_id <- paste(unique(this_patient_rct$`Subjects`), unique(this_patient_rct$Classification_ID), sep = '_') 
      print(unique_id)
      
      to_add_auc <- to_add_iauc <- to_add_dauc <- to_add_nauc <- c(unique(this_patient_rct$`Subjects`),
                                                                   unique(this_patient_rct$Classification))
      to_add_auc <- to_add_auc %>% append('AUC')
      to_add_iauc <- to_add_iauc %>% append('iAUC')
      to_add_nauc <- to_add_nauc %>% append('Net iAUC')
      to_add_dauc <- to_add_dauc %>% append('dAUC')
      
      time_series <- this_patient_rct$`Time in Minutes` %>% str_remove_all(' Min') %>% as.numeric()
      
      instance_df_ref_list <- list()
      
      for (j in 7:ncol(this_patient_rct)){ ## For every metabolite...
        mtbl_name <- colnames(this_patient_rct[j])
        mtbl_counts <- unlist(this_patient_rct[,j])
        
        if(any(is.na(mtbl_counts))){ ## NA/Inf Catcher
          mtbl_counts[is.na(mtbl_counts)] <- mean(mtbl_counts, na.rm = T)
        }
        if(any(is.infinite(mtbl_counts))){
          mtbl_counts[which(is.infinite(mtbl_counts))] <- mean(mtbl_counts[is.finite(mtbl_counts)], na.rm = T)
        }
        
        auc <- iauc <- dauc <- nauc <- vector()
        for (k in 1:(length(mtbl_counts)-1)){ ## For each timepoint...
          bl_count <- mtbl_counts[1]
          initial_count <- mtbl_counts[k]
          next_count <- mtbl_counts[k+1]
          time_diff <- time_series[k+1] - time_series[k]
          
          ## AUC
          auc <- auc %>% append((initial_count + next_count)/2 * time_diff) ## Whole area under the 2 points
          
          ## iAUC/dAUC/nAUC (incremental and 'decremental' and net AUC) === THINK triangles and trapeziums!
          if((initial_count >= bl_count & next_count >= bl_count)|
             (initial_count <= bl_count & next_count <= bl_count)){ ## if both points in question are above/below BL values
            initial_diff <- initial_count - bl_count
            next_diff <- next_count - bl_count
            my_area <- (initial_diff + next_diff)/2 * time_diff
            
            nauc <- nauc %>% append(my_area)
            if (my_area < 0){
              iauc <- iauc %>% append(0)
              dauc <- dauc %>% append(my_area)
            }else{
              iauc <- iauc %>% append(my_area)
              dauc <- dauc %>% append(0)
            }
            
            rm(initial_diff, next_diff, my_area)
            
          }else if((initial_count >= bl_count & next_count < bl_count)|
                   (initial_count < bl_count & next_count >= bl_count)){ ## if both points aren't all above/below BL values (one above one below)
            first_time_segment <- (initial_count - bl_count)/(initial_count - next_count) * time_diff 
            second_time_segment <- (next_count - bl_count)/(next_count - initial_count) * time_diff ## these time segments should both be +ve
            my_first_area <- first_time_segment * (initial_count - bl_count) / 2
            my_second_area <- second_time_segment * (next_count - bl_count) / 2
            
            nauc <- nauc %>% append(my_first_area + my_second_area)
            
            if (my_first_area > 0){ # and 2nd area is < 0
              iauc <- iauc %>% append(my_first_area)
              dauc <- dauc %>% append(my_second_area)
            }else if (my_second_area > 0){ # and 1st area is < 0
              iauc <- iauc %>% append(my_second_area)
              dauc <- dauc %>% append(my_first_area)
            }
            
            rm(first_time_segment, second_time_segment, my_first_area, my_second_area)
          }
          rm(initial_count, bl_count, next_count, time_diff, k)
        }
        
        df_area_reference <- data.frame(metabolite = rep(mtbl_name, 4), auc, iauc, dauc, nauc)
        instance_df_ref_list <- instance_df_ref_list %>% append(list(df_area_reference))
        
        to_add_auc[mtbl_name] <- sum(auc)
        to_add_iauc[mtbl_name] <- sum(iauc)
        to_add_dauc[mtbl_name] <- sum(dauc)
        to_add_nauc[mtbl_name] <- sum(nauc)
        
        rm(mtbl_name, mtbl_counts, j, auc, iauc, dauc, nauc, df_area_reference)
      }
      
      auc_df <- auc_df %>% rbind(to_add_auc, to_add_dauc, to_add_iauc, to_add_nauc)
      master_df_ref_list <- master_df_ref_list %>% append(list(unique_id, instance_df_ref_list)) # Should probably add a label...
      
      rm(instance_df_ref_list, to_add_auc, to_add_dauc, to_add_iauc, to_add_nauc,
         time_series, this_patient_rct, unique_id, i)
    }
    
    auc_df <- auc_df[-1,] ## Get rid of that filler first row
    
    rm(auc_berry)
    
#### Special Outlier Treatments (since we have negative numbers). followed by Mean + SE calculations
    auc_df[,c(5:ncol(auc_df))] <- lapply(auc_df[,c(5:ncol(auc_df))], as.numeric)
    auc_df_berry <- auc_df %>% group_by(AUC_Type) %>% group_split()
    
    master_mean_se_ci <- data.frame(Metabolites = colnames(auc_df)[5:ncol(auc_df)])
    auc_pvalues_fc <- data.frame(Metabolites = colnames(auc_df)[5:ncol(auc_df)])
    
    for (i in 1:length(auc_df_berry)){
      my_df <- auc_df_berry[[i]]
      auc_name <- unique(my_df$AUC_Type)
      
      ## Damn Outliers
      if (auc_name == 'AUC'){
        
        ## Log --> outlier remove --> unlog
        my_df[,c(5:ncol(my_df))] <- my_df[,c(5:ncol(my_df))] %>% log(base = 10)
        my_df <- my_df %>% outliers_begone(4, method = 'quantile')
        my_df[,c(5:ncol(my_df))] <- my_df[,c(5:ncol(my_df))] %>% lapply(function(x){return(10**x)})
        
      }else if (auc_name == 'dAUC'){
        
        ## Kill the 0s --> absolute --> log --> outlier remove --> unlog --> negate
        my_df[my_df == 0] <- NA
        my_df[,c(5:ncol(my_df))] <- (my_df[,c(5:ncol(my_df))] * -1) %>% log(base = 10) %>% outliers_begone(4) %>% lapply(function(x){return(-1 * (10**x))})
        
      }else if (auc_name == 'iAUC'){
        
        ## Kill the 0s --> log --> outlier remove --> unlog
        my_df[my_df == 0] <- NA
        my_df[,c(5:ncol(my_df))] <- my_df[,c(5:ncol(my_df))] %>% log(base = 10) %>% outliers_begone(4) %>% lapply(function(x){return(10**x)})
      }else if (auc_name == 'Net iAUC'){

        ## Standardization --> outlier remove --> destandardization
        for (j in 5:ncol(my_df)){

          ## Using min max normalization *****
          mtbl_list <- unlist(my_df[,j])
          mymin <- min(mtbl_list)
          mymax <- max(mtbl_list)
          mtbl_list <- mtbl_list %>% lapply(function(y){return ((y-mymin)/(mymax-mymin))}) %>% unlist()

          ## Separate outlier removal function for individual vectors...
          q1 <- as.numeric(quantile(mtbl_list, 0.05, na.rm = T))
          q2 <- as.numeric(quantile(mtbl_list, 0.95, na.rm = T))
          mtbl_list[which(mtbl_list < q1 | mtbl_list > q2)] <- NA

          ## Reverse min max *****
          mtbl_list <- mtbl_list %>% lapply(function(y){return (y*(mymax-mymin) + mymin)}) %>% unlist()

          my_df[,j] <- mtbl_list

          rm(mtbl_list, mymin, mymax, q1, q2, j)
        }
      }
      
      ## Plotting, Means & SE are to be done now
      
      ### <><><> Plotting done here (Toggle which sequence you wish to plot)
      mycontrasts <- list(c('P-E', 'P-F'), c('Q-E', 'Q-F'), c('P-E', 'Q-E'), c('P-F', 'Q-F'))
      ## Contrasts must be arranged in a list. Each contrast as elements within a vector.

      batch_png_plot(my_df, 'Classification', 4, toLog10 = F, remove_outliers = 'none', contrasts = mycontrasts,
                     y_label = auc_name, folder_name = paste0('2.4. GIN_', auc_name, '_T-tests_Contrasts_Indiv'),
                     test_type = 'wilcox.test', height = 8, width = 8)

      ### New Contrasts
      mycontrasts2 <- list(c('P', 'Q'), c('E', 'F'))
      my_df2 <- rbind(my_df %>% filter(grepl('P', Classification)), my_df %>% filter(grepl('Q', Classification)),
                      my_df %>% filter(grepl('E', Classification)), my_df %>% filter(grepl('F', Classification)))
      my_df2$Classification[1:length(my_df2$Classification)/2] <- my_df2$Classification[1:length(my_df2$Classification)/2] %>% str_remove_all('-.*')
      my_df2$Classification[(length(my_df2$Classification)/2 + 1) : length(my_df2$Classification)] <- my_df2$Classification[(length(my_df2$Classification)/2 + 1) : length(my_df2$Classification)] %>% str_remove_all('^.*-')

      batch_png_plot(my_df2, 'Classification', 4, toLog10 = F, remove_outliers = 'none', contrasts = mycontrasts2,
                     y_label = auc_name, height = 8, width = 8, test_type = 'wilcox.test',
                     folder_name = paste0('2.4. Glycemic_', auc_name, '_T-tests_Contrasts'))

      rm(contrasts2, my_df2)
      
      ### <><><> Mean and SE must be in each RCT group (Indicate 95% CI too)
      my_df_berry <- my_df %>% group_by(Classification) %>% group_split()
      
      for(k in 1:length(my_df_berry)){
        my_df2 <- my_df_berry[[k]]
        rct <- unique(my_df2$Classification)
        
        print(paste(auc_name, rct, sep = '_'))
        
        means <- se <- confi_interval <- vector()
        for (j in 5:ncol(my_df2)){
          mtbl_counts <- unlist(my_df2[,j])
          mtbl_mean <- mean(mtbl_counts, na.rm = T)
          mtbl_se <- sd(mtbl_counts, na.rm = T) / sqrt(length(mtbl_counts))
          mtbl_ci <- paste0(round((mtbl_mean - 1.96*mtbl_se), 2), '-', round((mtbl_mean + 1.96*mtbl_se), 2)) ## 95% to 1.96*SE
          
          means <- means %>% append(mtbl_mean)
          se <- se %>% append(mtbl_se)
          confi_interval <- confi_interval %>% append(mtbl_ci)
          
          rm(j, mtbl_counts, mtbl_mean, mtbl_se, mtbl_ci)
        }
        master_mean_se_ci <- master_mean_se_ci %>% cbind(means, se, confi_interval)
        names(master_mean_se_ci)[c((ncol(master_mean_se_ci)-2):ncol(master_mean_se_ci))] <- c(paste0(auc_name,'_', rct, '_means'),
                                                                                              paste0(auc_name,'_', rct, '_se'),
                                                                                              paste0(auc_name, '_', rct, '_95%_CI'))
        rm(means, se, my_df2, rct, confi_interval, k)
      }
      
      rm(my_df_berry)
      
      ### <><><> Performing T-Testings (Mann-Whitney)
      mycontrasts <- list(c('P-E', 'P-F'), c('Q-E', 'Q-F'), c('P-E', 'Q-E'), c('P-F', 'Q-F'),
                          c('P', 'Q'), c('E', 'F')) # Include the general groups.
      
      for (j in 1:length(mycontrasts)){ ### Isolating by contrasts
        first_element <- mycontrasts[[j]][1]
        second_element <- mycontrasts[[j]][2]
        
        my_df_1 <- my_df[which(grepl(first_element, my_df$Classification)),]
        my_df_2 <- my_df[which(grepl(second_element, my_df$Classification)),]
        
        pvalue_container <- fc_container <- vector()
        for (k in 5:ncol(my_df)){ ### All biomarkers
          first_contrast_list <- unlist(my_df_1[,k])
          second_contrast_list <- unlist(my_df_2[,k])
          
          if(sum(is.na(first_contrast_list)) / length(first_contrast_list) > 0.75 | ## Mass NA catcher
             sum(is.na(second_contrast_list)) / length(second_contrast_list) > 0.75){
            pvalue_container <- pvalue_container %>% append(NA)
            fc_container <- fc_container %>% append(NA)
            next()
            rm(k)
          }
          
          if(auc_name == 'AUC' | auc_name == 'Net iAUC'){ ### Only this would make sense... Q/P, F/E + Fold Change of AUC + Adjusted P-value (Labelling)
            fc <- mean(second_contrast_list, na.rm = T) / mean(first_contrast_list, na.rm = T)
            fc_container <- fc_container %>% append(fc)
            rm(fc)
          }
          
          pval <- wilcox.test(first_contrast_list, second_contrast_list)$p.value
          pvalue_container <- pvalue_container %>% append(pval)
          
          rm(k, pval, first_contrast_list, second_contrast_list)
        }
        
        auc_pvalues_fc <- auc_pvalues_fc %>% cbind(pvalue_container)
        colnames(auc_pvalues_fc)[ncol(auc_pvalues_fc)] <- paste0('pval_', auc_name, '_', first_element, 'vs', second_element)
        
        if(auc_name == 'AUC' | auc_name == 'Net iAUC'){
          auc_pvalues_fc <- auc_pvalues_fc %>% cbind(fc_container)
          colnames(auc_pvalues_fc)[ncol(auc_pvalues_fc)] <- paste0('fc_', auc_name, '_', first_element, 'vs', second_element)
        }
        
        rm(my_df_1, my_df_2, first_element, second_element, j, pvalue_container, fc_container)
      }
      
      rm(i, auc_name, my_df)
    }
    
    auc_pvalues_fc <- auc_pvalues_fc[,-c(14:25)] ## Remove dAUC and iAUC - who needs them?
    
    write.csv(auc_pvalues_fc, '2.4. Glycemic_MWUT_PVals_FC_AUC.csv', row.names = F)
    write.csv(master_mean_se_ci, '2.4. Glycemic_Mean_SE_CI_AUC.csv', row.names = F)
    write.csv(auc_df, '2.4. Glycemic_AUCs.csv', row.names = F)
    
    rm(master_mean_se_ci, mycontrasts, auc_df_berry)

### Testing for Net iAUC
    netauc_df <- auc_df %>% filter(AUC_Type == 'Net iAUC')
    test_contrasts <- list(c('P', 'Q'), c('E', 'F'))
    
    to_plot <- data.frame(Metabolites = colnames(netauc_df[5:ncol(netauc_df)]))
    
    for (contrast in test_contrasts){
      df1 <- netauc_df[which(grepl(contrast[1], netauc_df$Classification)),]
      df2 <- netauc_df[which(grepl(contrast[2], netauc_df$Classification)),]
      my_fc <- my_pvals <- vector()
      for (i in 5:ncol(netauc_df)){
        list1 <- df1[,i] %>% unlist()
        list2 <- df2[,i] %>% unlist()
        my_pvals <- my_pvals %>% append(wilcox.test(list1, list2)$p.value)
        fc <- mean(list2, na.rm = T) / mean(list1, na.rm = T)
        my_fc <- my_fc %>% append(fc)
        
        rm(list1, list2, i)
      }
      
      to_plot <- to_plot %>% add_column(my_pvals) %>% add_column(my_fc)
      names(to_plot)[(ncol(to_plot)-1):ncol(to_plot)] <- c(paste0('NetiAUC_Pvals_', contrast[2], 'vs', contrast[1]),
                                                           paste0('NetiAUC_FC_', contrast[2], 'vs', contrast[1]))
      rm(df1, df2)
    }
    to_plot$NetiAUC_Pvals_QvsP <- to_plot$NetiAUC_Pvals_QvsP %>% p.adjust(method = 'BH')
    to_plot$NetiAUC_Pvals_FvsE <- to_plot$NetiAUC_Pvals_FvsE %>% p.adjust(method = 'BH')
    
    plot_volcano(to_plot[,c(1, 4, 5)], volc.title = 'F vs E', x.title = 'Log2(Fold Change of Net iAUC) [F/E]',
                 y.title = '-lg(P-Value) [BH-Adjusted]', toLog2FC = T, toNegLog10Pval = T, x_limit = c(-15, 15),
                 text_size = 3, repel_labels = F)
    
### ===== Phase 2.4.1. P-Adjusted (BH) Values =====
    pval_colnums <- which(grepl('pval', colnames(auc_pvalues_fc)))
    auc_pvalues_fc_adj <- auc_pvalues_fc
    
    for(i in pval_colnums){
      column_name <- colnames(auc_pvalues_fc)[i]
      bh_pvals <- unlist(auc_pvalues_fc[,i]) %>% p.adjust(method = 'BH')
      auc_pvalues_fc_adj <- auc_pvalues_fc_adj %>% add_column(bh_pvals, .after = column_name)
      names(auc_pvalues_fc_adj)[which(colnames(auc_pvalues_fc_adj) == 'bh_pvals')] <- paste0('P.Adjusted(BH)_', column_name)
      
      rm(column_name, bh_pvals, i)
    }        
    
    write.csv(auc_pvalues_fc_adj, '2.4.1. Glycemic_Pvalues_BH_FC.csv', row.names = F)
    
### ==== Phase 2.4.2 Calculating the Volcano Plots ====
    ## (P-value & FC across contrasts - P-Q, E-F)
    ## Following up from previous section... we need auc_pvalues_fc
    adj_order <- 1
    for (i in seq(2, ncol(auc_pvalues_fc_adj), 3)){
      adj_order <- adj_order %>% append(c(i+1, i+2))
    }
    auc_volc_adj <- auc_pvalues_fc_adj[,adj_order]
    rm(adj_order)
    
    for (i in seq(2, length(auc_volc_adj), 2)){
      j <- i+1
      
      my_df_adj <- auc_volc_adj[,c(1, i, j)]
      my_df_adj[,2] <- log(my_df_adj[,2], base = 10)*-1
      
      if(any(grepl('Net iAUC', colnames(my_df_adj)))){
        my_df_adj_nauc <- my_df_adj
      }
      
      my_df_adj[,3] <- log(my_df_adj[,3], base = 2) # Note: Net iAUC may also contain negative values which cannot be logged.
      
      df_names <- colnames(my_df_adj)[2] %>% str_remove_all('^.*_') %>% str_split('vs') %>% unlist()
      auc_type <- colnames(my_df_adj)[2] %>% str_remove_all('.*pval_') %>% str_remove_all('_.*')
      file_name <- paste0('(', auc_type, ') ', df_names[2], ' vs ', df_names[1])
      df_title <- paste0(df_names[2], '/', df_names[1])
      x_title <- paste0('log2(Fold Change of ', auc_type, ') [', df_title, ']')
      y_title <- '-lg(P-Value) [BH-Adjusted]'
      print(paste0(auc_type, ' ', df_title))
      
      read_plot_volcano(dataframe = my_df_adj, volc.title = file_name,
                        x.title = x_title, y.title = y_title, text_size = 3,
                        fc.threshold.left = -0.75, fc.threshold.right = 0.75, repel_labels = F) # Volcy for adjusted p-values
      
      if(exists('my_df_adj_nauc')){
        x_title_nauc <- x_title %>% str_remove_all('log2\\(') %>% str_remove('\\)')
        file_name <- file_name %>% paste0(' (Non-Logged FC)')
        read_plot_volcano(dataframe = my_df_adj_nauc, volc.title = file_name,
                          x.title = x_title_nauc, y.title = y_title, text_size = 3,
                          fc.threshold.left = -0.75, fc.threshold.right = 0.75, x_limit = c(-15, 15), repel_labels = F)
        
        rm(my_df_adj_nauc)
      }
      
      rm(i, j, my_df_adj, df_title, df_names, file_name, x_title, y_title)
    }
    
    file.rename('Volcano Plot Outputs', '2.4.2. AUC Volcano Plots')
    write.csv(auc_volc_adj, "2.4.2. (Net)AUC Volcano Plotting Data.csv", row.names = F)
    
### ==== Phase 2.5. PCA Analysis ====
    library(ggfortify)
### -=-=-=-=- EF -=-=-=-=-
    qtc_forpca <- qtc[,-c(1:2, 4:6)]
    qtc_forpca[is.na(qtc_forpca)] <- 0
    
    pvalmaster_EF <- data.frame(mtbl = auc_pvalues_fc_adj$Metabolites, pvals_EF = auc_pvalues_fc_adj$`P.Adjusted(BH)_pval_AUC_EvsF`)
    pvalmaster_EF <- pvalmaster_EF[order(pvalmaster_EF$pvals_EF),]
    candidates_EF <- pvalmaster_EF$mtbl[c(1:50)]
    
    # Filtering out for outlier samples/biomarkers
    # qtc_forpca <- qtc_forpca[,-which(colnames(qtc_forpca) %in% c('Metabolite 43'))]
    
    qtc_forpca <- qtc_forpca[-which(row.names(qtc_forpca) %in% c(410, 778, 790, 551, 94, 96, 93,
                                                                 676, 95)),]
    
    qtc.pca <- prcomp(qtc_forpca[,c(which(colnames(qtc_forpca) %in% candidates_EF))])
    summary(qtc.pca)
    p1 <- autoplot(qtc.pca, data = qtc_forpca, colour = 'Classification', label = T, label.size = 4,
             loadings = F, loadings.label = T, loadings.label.size = 4,
             scale = 0, frame = T, frame.type = 'norm')
    p1
    
    output_plot(p1, foldername = '2.5. PCA Plots', filename = 'PCA Plot - E vs F')

### -=-=-=-=- PQ -=-=-=-=-
    qtc_forpca <- qtc[,-c(1:2, 4:6)]
    qtc_forpca[is.na(qtc_forpca)] <- 0
    
    pvalmaster_PQ <- data.frame(mtbl = auc_pvalues_fc_adj$Metabolites, pvals_PQ = auc_pvalues_fc_adj$`P.Adjusted(BH)_pval_AUC_PvsQ`)
    pvalmaster_PQ <- pvalmaster_PQ[order(pvalmaster_PQ$pvals_PQ),]
    candidates_PQ <- pvalmaster_PQ$mtbl[c(1:50)]
    
    ## Filter the outliers
    candidates_PQ <- candidates_PQ[-which(candidates_PQ %in% c(
    'Metabolite 43',
    'Metabolite 44'
    ))]
    qtc_forpca <- qtc_forpca[-which(row.names(qtc_forpca) %in% c(93, 94, 95, 96, 98, 324, 551, 676)),]
    
    qtc.pca <- prcomp(qtc_forpca[,c(which(colnames(qtc_forpca) %in% candidates_PQ))])
    summary(qtc.pca)
    p2 <- autoplot(qtc.pca, data = qtc_forpca, colour = 'Classification', label = T, label.size = 4,
                   loadings = T, loadings.label = T, loadings.label.size = 4,
                   scale = 0, frame = T, frame.type = 'norm')
    p2
    output_plot(p2, foldername = '2.5. PCA Plots', filename = 'PCA Plot - P vs Q')

### ==== Phase 2.6. Correlation Analysis ====
    
    qtc2 <- qtc %>% outliers_begone(7, method = 'quantile')
    qtc2$`Time in Minutes` <- qtc2$`Time in Minutes` %>% str_remove_all(' Min') %>% as.numeric() # Correct ordering for Timing
    ### For individual RCT groups & by individual time
    myclassifications <- c('Overall', 'P', 'Q', 'P-E', 'P-F', 'Q-E', 'Q-F')
    
    master_corr_df <- master_timecorr_df <- data.frame(Metabolites = colnames(qtc2)[7:ncol(qtc2)])
    for(classification in myclassifications){
      
      if(classification == 'Overall'){
        my_df <- qtc2
      }else{
        my_df <- qtc2[which(grepl(classification, qtc2$Classification)),]
      }
      
      ### For individual RCT groups
      corr_df <- behold_correlations(my_df, 6, 'Glucose', method = 'pearson')
      names(corr_df)[c(2,3)] <- paste0(classification, '_', names(corr_df[c(2,3)]), '_withGlucose')
      master_corr_df <- master_corr_df %>% add_column(corr_df[,c(2,3)])
      
      ### For individual timepoints
      my_df_berry <- my_df %>% group_by(`Time in Minutes`) %>% group_split()
      for (i in 1:length(my_df_berry)){
        my_df_time <- my_df_berry[[i]]
        timepoint <- unique(my_df_time$`Time in Minutes`)
        corr_time_df <- behold_correlations(my_df_time, 6, 'Glucose', method = 'spearman')
        names(corr_time_df)[c(2,3)] <- paste(classification, timepoint, names(corr_time_df[c(2,3)]), 'withGlucose', sep = '_')
        master_timecorr_df <- master_timecorr_df %>% add_column(corr_time_df[,c(2,3)])
        
        print(paste(classification, timepoint, sep = '_'))
        rm(timepoint, i, my_df_time, corr_time_df)
      }
      
      rm(classification, my_df, corr_df, my_df_berry)
    }
    
    ## We can replace qtc2 with qtc_log to see the correlation between logged values. (With outliers removed)
    
    write.csv(master_corr_df, '2.6. Glycemic-Correlation of Metabolites with Glucose Level.csv', row.names = F)
    write.csv(master_timecorr_df, '2.6. Glycemic-Correlation of Metabolites with Glucose Level (Timepoints).csv', row.names = F)
    
    rm(qtc2, master_corr_df, master_timecorr_df, myclassifications)
    
### ==== Phase 2.6.1, T-Tests between 0 and 60 Min ====
    qtc2 <- qtc_log %>% outliers_begone(7) %>% filter(`Time in Minutes` == '0 Min' | `Time in Minutes` == '60 Min')
    
    myclassifications <- c('Overall', 'P', 'Q', 'P-E', 'P-F', 'Q-E', 'Q-F')
    mycontrasts <- list(c('\\b0 Min', '\\b60 Min'))
    
    pval_0to60_master <- data.frame(colnames(qtc2)[7:ncol(qtc2)])
    names(pval_0to60_master) <- 'Metabolites'
    
    for(classification in myclassifications){
      if(classification == 'Overall'){
        my_df <- qtc2
      }else{
        my_df <- qtc2[which(grepl(classification, qtc2$Classification)),]
      }
      print(classification)
      
      p_values <- a_stat_test(my_df, 'Time in Minutes', 6, mycontrasts = mycontrasts, grepl_fixed = F, paired = T, colName_SortPairs = 'Subjects')
      pval_0to60_master <- pval_0to60_master %>% cbind(p_values[,2])
      names(pval_0to60_master)[ncol(pval_0to60_master)] <- paste0(classification, '_0Minvs60Min_PVals')
      
      rm(classification, my_df, p_values)
    }
    
    write.csv(pval_0to60_master, '2.6.1. Glycemic - PValues(0vs60Min).csv', row.names = F)
    
    rm(pval_0to60_master, qtc2, myclassifications, mycontrasts)
    
### ==== Phase 2.7. Visualization: Line graphs over time for each biomarker ====
    
### Coming up with the organized plotting data post-median calculation
    
    qtc_log_plotting_berry <- qtc_log %>% group_by(Classification, `Time in Minutes`) %>% group_split()
    glucose_medians <- classi <- time_min <- 'temp'
    mtbl_plotter <- t(data.frame(rows = colnames(qtc_log_plotting_berry[[1]][8:ncol(qtc_log_plotting_berry[[1]])])))
    for(i in 1:length(qtc_log_plotting_berry)){
      df <- qtc_log_plotting_berry[[i]]
      classi <- classi %>% append(unique(df$Classification))
      time_min <- time_min %>% append(unique(df$`Time in Minutes`))
      glucose_medians <- glucose_medians %>% append(median(df$Glucose, na.rm = T))
      
      meddy_list <- vector()
      for (j in 8:ncol(df)){ ### Plugging into mtbl_plotter first, then combining the other classification columns
        mymedian <- median(unlist(df[,j]), na.rm = T)
        meddy_list <- meddy_list %>% append(mymedian)
        rm(mymedian, j)
      }
      
      mtbl_plotter <- mtbl_plotter %>% rbind(meddy_list)
      rm(df, i, meddy_list)
    }
    qtc_log_plotting <- data.frame(cbind(classi, time_min, glucose_medians, mtbl_plotter))
    rm(classi, time_min, glucose_medians, mtbl_plotter)
    
    names(qtc_log_plotting)[1:3] <- c('Classification', 'Time in Minutes', 'Glucose')
    names(qtc_log_plotting)[4:ncol(qtc_log_plotting)] <- qtc_log_plotting[1,][4:ncol(qtc_log_plotting)]
    qtc_log_plotting <- qtc_log_plotting[-1,]
    rownames(qtc_log_plotting) <- NULL
    qtc_log_plotting[3:ncol(qtc_log_plotting)] <- lapply(qtc_log_plotting[3:ncol(qtc_log_plotting)],
                                                         function(x){return(as.numeric(x))})
    
    qtc_log_plotting <- within(qtc_log_plotting, c(`Time in Minutes` <-  factor(`Time in Minutes`, levels = c("0 Min", "30 Min", "60 Min", '120 Min', '180 Min')),
                                                   Classification <-  as.factor(Classification)))
    names(qtc_log_plotting)[2] <- 'Time'
    qtc_log_plotting$Time <- qtc_log_plotting$Time %>% str_remove_all(' Min') %>% as.numeric()
    
### The Plotting

    ## Testplot
    plot_graph_line(qtc_log_plotting, 'Classification', 'Time', 'Metabolite 77', x_label = 'Time in Minutes')
    
    for (i in 3:ncol(qtc_log_plotting)){
      mtbl <- colnames(qtc_log_plotting)[i]
      print(mtbl)
      if(mtbl == 'Glucose'){mylabel <-  ''}else{mylabel <- 'Log10(Counts)'}
      
      ### Name cleaning for filename producing
      mtbl_name <- mtbl  ## We create a separate container for our filename because...
      mtbl_name <- mtbl_name %>% str_replace_all('/', ' Or ')  ## Only the mtbl names are affected by the slash or colon
      mtbl_name <- mtbl_name %>% str_replace_all(':', '.') ## Our graph labels/titles need not have them

      if(!file.exists('2.7. Lineplots_GIN')){dir.create('2.7. Lineplots_GIN')}
      png(filename = paste0('2.7. Lineplots_GIN/', mtbl_name, '.png'), width = 11, height = 8, units = 'in', res = 300)
      p1 <- plot_graph_line(qtc_log_plotting, 'Classification', 'Time', mtbl_name = mtbl, y_label = mylabel)
      print(p1)
      dev.off()
      
      rm(mtbl, mtbl_name, p1, i)
    }
    
    rm(qtc_log_plotting, qtc_log_plotting_berry)

### ==== Phase 2.8 Correlation and Linear Regression ====
    ### Trying to find correlations between glucose values at all timepoints and baseline biomarker levels.
    qtc_corr <- qtc_log %>% outliers_begone(7)
    qtc_corr$`Time in Minutes` <- qtc_corr$`Time in Minutes` %>% str_remove_all(' Min') %>% as.numeric()
    qtc_corr_berry <- qtc_corr %>% group_by(`Time in Minutes`) %>% group_split()
    
    myclassifications <- c('P', 'Q', 'E', 'F', 'P-E', 'P-F', 'Q-E', 'Q-F')
    bl_values <- qtc_corr_berry[[1]]
    master_pearson <- master_spearman <- master_lm <- data.frame(Metabolites = colnames(bl_values)[8:ncol(bl_values)], ' ' = '')
    for(classif in myclassifications){ ## For every classification/RCT group...
      bl_values_rct <- bl_values[which(grepl(classif, bl_values$Classification)),]
      for(my_df in qtc_corr_berry){ ## For every timepoint in qtc_corr_berry...
        timepoint <- unique(my_df$`Time in Minutes`) %>% paste0('-Min')
        my_glucose <- my_df[which(grepl(classif, my_df$Classification)),]$Glucose ## Timepoint, RCT specific for glucose.
        bl_values_rct$Glucose <- my_glucose
        
        print(paste(classif, timepoint, 'Glucose', sep = '_'))
        
        ### Performing correlations using behold_correlations function
        corr_df_pearson <- behold_correlations(bl_values_rct, 7, 'Glucose', method = 'pearson')
        corr_df_spearman <- behold_correlations(bl_values_rct, 7, 'Glucose', method = 'spearman')
        
        names(corr_df_pearson)[c(2,3)] <- names(corr_df_pearson)[c(2,3)] %>% paste('Glucose', timepoint, classif, sep = '_')
        names(corr_df_spearman)[c(2,3)] <- names(corr_df_spearman)[c(2,3)] %>% paste('Glucose', timepoint, classif, sep = '_')
        
        master_pearson <- master_pearson %>% cbind(corr_df_pearson[,c(2,3)]) %>% cbind(' ' = '') ## Creating an empty column for separation of categories.
        master_spearman <- master_spearman %>% cbind(corr_df_spearman[,c(2,3)]) %>% cbind(' ' = '')
        
        ### Performing LM
        coef <- se <- pval <- rsquare <- vector()
        for(i in 8:ncol(bl_values_rct)){
          mylm <- summary(lm(bl_values_rct$Glucose ~ unlist(bl_values_rct[,i])))
          coef <- coef %>% append(mylm$coefficients[2,1]) 
          se <- se %>% append(mylm$coefficients[2,2])
          pval <- pval %>% append(mylm$coefficient[2,4])
          rsquare <- rsquare %>% append(mylm$r.squared)
          rm(mylm, i)
        }
        
        master_lm <- master_lm %>% cbind(coef, se, pval, rsquare)
        names(master_lm)[(ncol(master_lm)-3):ncol(master_lm)] <- names(master_lm)[(ncol(master_lm)-3):ncol(master_lm)] %>% paste(
          'Glucose', timepoint, classif, sep = '_'
        )
        master_lm <- master_lm %>% cbind(' '= '')
        
        rm(corr_df_pearson, corr_df_spearman, timepoint, my_glucose, my_d, coef, se, pval, rsquare)
      }
      rm(bl_values_rct)
    }
    
    names(master_lm)[2] <- names(master_spearman)[2] <- names(master_pearson)[2] <- ''
    
    write.csv(master_pearson, '2.8. Glycemic - Pearson Correlation (BL Biomarkers vs Glucose at each Timepoint).csv', row.names = F)
    write.csv(master_spearman, '2.8. Glycemic - Spearman Correlation (BL Biomarkers vs Glucose at each Timepoint).csv', row.names = F)
    write.csv(master_lm, '2.8. Glycemic - LM (BL Biomarkers vs Glucose at each Timepoint).csv', row.names = F)

### ==== Phase X.1. Using mean comparison  ====
## Calculates the mean instead of AUC & variants for each patient's levels of each biomarker
    
    ### Container creation
    qtc_mean <- qtc[,c(1, 3, 4)]
    qtc_mean$ToSort <- paste(qtc_mean$`Subjects`, qtc_mean$Classification, qtc_mean$Classification_ID, sep = '-')    
    qtc_mean <- unique(qtc_mean$ToSort) ## Container for the means
    qtc_mean <- data.frame(t(data.frame(lapply(qtc_mean, function(x){return (strsplit(x, '-'))}))))
    row.names(qtc_mean) <- NULL
    names(qtc_mean) <- c('Subjects', 'Circadian', 'GI Level', 'Classification_ID')
    
    ### Filling it up
    for (i in 10:ncol(qtc)){
      data_list <- qtc[,i]
      mtbl_name <- colnames(qtc)[i]
      
      mymeans <- vector()
      for (j in seq(1, length(data_list), 5)){
        mymeans <- mymeans %>% append(mean(data_list[j:(j+4)], na.rm = T))
      }
      
      meaned_df <- data.frame(mymeans)
      names(meaned_df) <- mtbl_name
      qtc_mean <- qtc_mean %>% cbind(meaned_df)
      
      rm(i, j, mymeans, data_list, mtbl_name, meaned_df)
    }
    
    qtc_mean <- qtc_mean[order(qtc_mean$Classification_ID),]
    qtc_mean <- qtc_mean %>% add_column(Classification = paste(qtc_mean$Circadian, qtc_mean$`GI Level`, sep = '-'), .before = 'Classification_ID')
    
    ### BoxPlotting of the 4 RCT Groups
    mycontrasts <- list(c('P-E', 'P-F'), c('Q-E', 'Q-F'), c('P-E', 'Q-E'), c('P-F', 'Q-F'))
    plot_graph(qtc_mean, 'Classification', 6, toLog10 = F, remove_outliers = 'none',
               palette = c('deeppink', 'aquamarine3', 'purple', 'cyan3'),
               contrasts = mycontrasts)
    batch_png_plot(qtc_mean, 'Classification', 6, toLog10 = T, remove_outliers = 'quantile',
                   palette = c('deeppink', 'aquamarine3', 'purple', 'cyan3'),
                   contrasts = mycontrasts)
    file.rename('plotsOutputDir', 'X.1. Glycemic - Barplots of Means')
    rm(mycontrasts)
    
    ### Isolating P-Values
    mycontrasts <- list(c('E', 'F'), c('P', 'Q'), c('P-E', 'P-F'),
                        c('Q-E', 'Q-F'), c('P-E', 'Q-E'), c('P-F', 'Q-F'))
    qtc_mean_log <- cbind(qtc_mean[,c(1:6)], log(qtc_mean[,c(7:ncol(qtc_mean))], base = 10)) %>% outliers_begone(6)
    df_pvals <- a_stat_test(qtc_mean_log, 'Classification', 5, mycontrasts = mycontrasts, test_type = 't.test')

    write.csv(df_pvals, 'X.1. GIN-Pvals of RCT Groups (Mean).csv', row.names = F)  
    
    rm(df_pvals, qtc_mean_log, qtc_mean)
    
