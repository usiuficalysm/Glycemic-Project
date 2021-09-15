# importFrom(magrittr, "%>%") to be included in NAMESPACE
# May as well import the whole tidyverse for good measure (ggplot2, stringr, dplyr etc.)
# ggpubr is also a possibility

###### ============ MODIFIED: 030821 1745hrs =============
# ^^^^^^^^^

### ==== Outlier Removal ====
## Individual vector outlier removal

#' Removes Outliers for Individual Vectors
#'
#' @param myvector Numeric Vector. Vector of data to have outliers removed
#' @param method String. Method of outlier removal (eg. 0.05-0.95 Quantile, 2SD etc.) Default = 'quantile'.
#'
#' @return A vector with outliers removed.
#' @export
#'
#' @examples outliers_bye(data_vector, method = 'quantile')
outliers_bye <- function(myvector, method = 'quantile'){
  if(method == 'quantile'){
    q1 <- as.numeric(quantile(myvector, 0.05, na.rm = T))
    q2 <- as.numeric(quantile(myvector, 0.95, na.rm = T))
    myvector[which(myvector < q1 | myvector > q2)] <- NA
    
    return(myvector)
    
  }else if (method == 'sd'){
    sd_col <- sd(myvector, na.rm = T)
    mean_col <- mean(myvector, na.rm = T)
    outliers <- which(abs(myvector - mean_col) > 2*sd_col) ## Toggle your SD filter here!
    myvector[outliers] <- NA
    myvector[myvector == Inf] <- NA
    
    return(myvector)
    
  }else{
    warning("Outlier method is unknown. Please specify the method.")
  }
}

### Note: This function works for removing outliers from ALL columns (after info_cols) in the dataset

#' Remove Outliers for Dataframe
#'
#' @param dataset Dataframe. Dataframe containing data of multiple variables to be dealt a hand in outlier removal. Usually derived from the quant_table output from MRMkit.
#' @param num_info_col Integer. Number of columns containing metadata.
#' @param method String. Method of outlier removal (eg. 0.05-0.95 Quantile, 2SD etc.) Default = 'quantile'.
#'
#' @return The same dataframe with the vectors of all variables having outliers removed.
#' @export
#'
#' @examples outliers_begone(df, num_info_col = 3, method = 'quantile')
outliers_begone <- function(dataset, num_info_col, method = 'quantile'){
  if(method == 'quantile' | method == 'sd'){
    for(i in (num_info_col+1):ncol(dataset)){
      mtbl_vector <- unlist(dataset[,i])
      mtbl_vector <- mtbl_vector %>% outliers_bye(method = method)
      dataset[,i] <- mtbl_vector
    }
    return (dataset)
  }else{
    warning("Outlier method is unknown. Please specify the method. Returning the original dataframe...")
    return(dataset)
  }
}

### ==== Plotters ====
#' Plot a Single-Variable Graph from Multi-Variable Dataset.
#'
#' @param qtc Dataframe. Dataframe to be used for plotting.
#' @param colName_toSortBy String. Column name containing the different categories of data to apply contrasts on.
#' @param mtbl_col_num Integer. Column number of variable in question. 
#' @param toLog10 Boolean. Do you want to perform a Log10 transformation on the data? Default = False.
#' @param remove_outliers String. Indicate method for outlier_removal c('sd', 'quantile'). Default = 'sd'. Indicate gibberish for no outlier-removal.
#' @param palette String Vector. Main palette for coloring the different categories.
#' @param palette_fill String Vector. Secondary palette for special graph types. (eg. Multi-barplot)
#' @param plot_type String. c('std_boxplot', 'violin', 'multi_boxplot'). Default = 'std_boxplot'.
#' @param colName_toFill String. Secondary categorical variable for special graph types. (eg. Multi-barplot)
#' @param contrasts List. Input a list of contrasts to be compared in stat_compare_means. eg. (list(c('Type A', 'Type B'), c('Type A', 'Type C'))). Default = NULL.
#' @param test_type String. c('t.test', 'wilcox.test', 'anova', 'kruskal.test') The test type to be used in stat_comapare means. Default = 't.test'.
#' @param y_label String. y_axis label. Default = 'Counts'.
#' @param title_size Integer. Title size. Default = 36.
#' @param y_axis_text_size Integer. y-axis text size. Default = 11.
#' @param x_axis_text_size Integer. x-axis text size. Default = 18.
#' @param y_axis_title_size Integer. y-axis title size. Default = 14.
#' @param legend_position String. Where do you want the legend to be? Default is 'none'.
#'
#' @return ggplot2 object of the desired plot of a target variable.
#' @export
#'
#' @examples plot_graph(qtc, 'Classification', 132, toLog10 = T, remove_outliers = 'sd', plot_type = 'std_boxplot', contrasts = list(c('A', 'B'), c('A', 'C')), test_type = 't.test')
plot_graph <- function(qtc, colName_toSortBy, mtbl_col_num, toLog10 = F, remove_outliers = 'sd',
                       palette = c('#FF61F5','#FF3587','#FF3131', '#FFA56A', '#FFC101', '#81d100',
                                   '#00a827', '#00d692', '#00c9c9', '#1495ff', '#1434ff','#7900db'),
                       palette_fill = c('aquamarine', 'deeppink', 'cyan1', 'purple',
                                        'sienna1', 'seagreen2', 'cornflowerblue', 'red3'),
                       plot_type = 'std_boxplot', colName_toFill = NULL,
                       contrasts = NULL,
                       test_type = 't.test',
                       y_label = 'Counts',
                       title_size = 36, y_axis_text_size = 11, x_axis_text_size = 18, y_axis_title_size = 14,
                       legend_position = 'none'){
  require(tidyverse) # Mainly ggplot2, dplyr, stringr
  require(ggpubr)   # stat_compare_means
  
  ### -- Name cleaning --
  # names(qtc) <- str_replace_all(names(qtc), '/', ' Or ')
  names(qtc) <- str_replace_all(names(qtc), ':', '.')
  
  ### -- Label naming/Transformations --
  test_name <- colnames(qtc[mtbl_col_num])
  new_name <- str_remove_all(test_name, '.Ratio.')
  print(new_name)
  test_name <- paste0('`', test_name, '`')
  
  names(qtc)[names(qtc) == colName_toSortBy] <- 'Contrast'
  qtc$Contrast <- factor(qtc$Contrast, unlist(unique(qtc$Contrast)), ordered = T)
  
  if(!is.null(colName_toFill)){
    names(qtc)[names(qtc) == colName_toFill] <- 'Classifier'
    qtc$Contrast_Fill <- factor(qtc$Contrast, unlist(unique(qtc$Contrast_Full)), ordered = T)
  }
  
  ### -- Creating comparison lists for the contrasts to be made when plotting (for stat.tests) --
  if(class(contrasts) != 'list' & !is.null(contrasts)){
    if(contrasts == 'all'){
      factor_list <- levels(qtc$Contrast)
      comparison_list <- list()
      for(i in 1:(length(factor_list)-1)){
        myvector <- factor_list[i]
        for(j in (i+1):length(factor_list)){
          myvector2 <- myvector %>% append(factor_list[j])
          comparison_list <- comparison_list %>% c(list(myvector2))
        }
      }
    }else{break}
  }else if(!is.null(contrasts)){
    comparison_list <- contrasts
  }else{
    comparison_list <- NA
  }
  
  
  ### -- Log10 or Not --
  if(toLog10){
    qtc[mtbl_col_num] <- log(qtc[mtbl_col_num], base = 10)
    y_label <- paste0('Log10(', y_label, ')')
  }
  
  ### -- Outlier Removal --
  if(tolower(remove_outliers) == 'sd'){
  ### SD Chopping
      sd_col <- sd(unlist(qtc[mtbl_col_num]), na.rm = T)
      mean_col <- mean(unlist(qtc[mtbl_col_num]), na.rm = T)
      outliers <- which(abs(qtc[mtbl_col_num] - mean_col) > 2*sd_col) ## Toggle your SD filter here!
      argh <- unlist(qtc[mtbl_col_num])
      argh[outliers] <- NA
      argh[argh == Inf] <- NA
      qtc[mtbl_col_num] <- argh
  }else if (remove_outliers == 'quantile'){
  ### Quantile/Percentile Filtering
      leftbound <- quantile(unlist(qtc[mtbl_col_num]), 0.05, na.rm = T) ## Toggle your Quantile Threshold here!
      rightbound <- quantile(unlist(qtc[mtbl_col_num]), 0.95, na.rm = T)
      outliers <- which(unlist(qtc[mtbl_col_num]) < leftbound | unlist(qtc[mtbl_col_num]) > rightbound)
      argh <- unlist(qtc[mtbl_col_num])
      argh[outliers] <- NA
      argh[argh == Inf] <- NA
      qtc[mtbl_col_num] <- argh
  }else{
    warning('WARNING: Outlier removal method is unspecified/unclear. It will not be performed. Sucks to be you. :(')
  }
  
  ### -- Initializing plot types --
  ## Standard Plots
  std_boxplot <- geom_boxplot(fill = NA, colour = 'gray6', lwd = 1.2, width = 0.1, outlier.color = 'cyan', outlier.shape = NA)
  std_pointplot <- geom_point(size = 1.5, position = position_jitter(width = 0.075), aes(color = Contrast))
  std_violinplot <- geom_violin(trim = FALSE, aes(color = Contrast), color = NA, alpha = 0.3)
  
  ## Standard Colorations
  std_colors <- scale_color_manual(values = palette)
  std_fill_colors <- scale_fill_manual(values = palette_fill)
 
  ## Multi Boxplot
  multi_boxplot <- geom_boxplot(colour = 'gray6', lwd = 1.2, width = 0.4, outlier.color = 'cyan',
                                outlier.shape = NA, alpha = 0.3)
  multi_fill_colors <- scale_fill_manual(name = 'Classification', values = palette_fill)
  
  ## Unused
  std_dotplot <- geom_dotplot(binaxis = 'y', stackdir = 'centerwhole', dotsize = 0.5,
                              position = position_jitter(width = 0, height = 0.05))
  
  
  ### -- Choosing the Plot Contrast --
  switch(plot_type,
         std_boxplot = {
           element1 <- std_pointplot
           element2 <- std_boxplot
           color1 <- std_colors
           fill1 <- std_fill_colors
         },
         violin = {
           element1 <- std_violinplot
           element2 <- std_pointplot
           color1 <- std_colors
           fill1 <- std_fill_colors
         },
         multi_boxplot = {
           element1 <- multi_boxplot
           element2 <- multi_boxplot
           color1 <- std_colors
           fill1 <- multi_fill_colors 
         }
         )
  
  if(is.null(colName_toFill)){
    fill_type <- 'Contrast'
  }else{
    fill_type <- 'Classifier'
  }
  
  
  ### == Actual Plot ==
  p1 <- ggplot(qtc, aes(x = Contrast, y = eval(parse(text = test_name)), fill = get(fill_type))) +
    theme_minimal() +
    ggtitle(new_name) + theme(plot.title = element_text(size = title_size, hjust = 0.5, face = 'bold'),
                              axis.text.x = element_text(size = x_axis_text_size, face = 'bold'),
                              axis.title.y = element_text(size = y_axis_title_size, face = 'bold'),
                              axis.text.y = element_text(size = y_axis_text_size, face = 'bold'),
                              legend.position = legend_position) +
    xlab('') + ylab(y_label) +
    element1 + element2 +
    stat_compare_means(comparisons = comparison_list, method = test_type, ### Unfortunately, this had to be added directly.
                       bracket.size = 1.15, size = 4.5, na.rm = T,
                       vjust = -0.2) +
    stat_compare_means(comparison = comparison_list, method = test_type, ### For the Labels (eg. ***)
                       label = 'p.signif', size = 8, hide.ns = T,
                       vjust = 1.45) +
    color1 + fill1
  
  return(p1)
}

## Note that names have to fulfill the file-naming criterion in Windows (eg. no ':', '/' etc.).

#' Creates a Folder and Outputs Single-Variable Graphs for All Variables in a Dataset with Metadata (left-aligned).
#' 
#' @param qtc Dataframe. Dataframe to be used for plotting.
#' @param colName_toSortBy String. Column name containing the different categories of data to apply contrasts on.
#' @param num_info_col Integer. Number of columns containing metadata (non-variable columns)
#' @param toLog10 Boolean. Do you want to perform a Log10 transformation on the data? Default = False.
#' @param remove_outliers String. Indicate method for outlier_removal c('sd', 'quantile'). Default = 'sd'. Indicate gibberish for no outlier-removal.
#' @param palette String Vector. Main palette for coloring the different categories.
#' @param palette_fill String Vector. Secondary palette for special graph types. (eg. Multi-barplot)
#' @param contrasts List. Input a list of contrasts to be compared in stat_compare_means. eg. (list(c('Type A', 'Type B'), c('Type A', 'Type C'))). Default = NULL.
#' @param test_type String. c('t.test', 'wilcox.test', 'anova', 'kruskal.test') The test type to be used in stat_comapare means. Default = 't.test'.
#' @param y_label String. y_axis label. Default = 'Counts'.
#' @param height Integer. Height of plot picture output in specified units. Default is 8 (in).
#' @param width Integer. Width of plot picture output in specified units. Default is 11 (in).
#' @param units String. c('in', 'px', 'cm', 'mm') Specified units for plot picture dimensions. Default is 'in'.
#' @param res Integer. Resolution in ppi recorded in bitmap file. Default is 300.
#' @param title_size Integer. Title size. Default = 36.
#' @param y_axis_text_size Integer. y-axis text size. Default = 11.
#' @param x_axis_text_size Integer. x-axis text size. Default = 18.
#' @param y_axis_title_size Integer. y-axis title size. Default = 14.
#' @param legend_position String. Where do you want the legend to be? Default is 'none'.
#' @param folder_name String. Intended folder name for output plots in the same directory as script. Default is 'plotsOutputDir'.
#'
#' @return A folder created in the same script directory containing the Single-Variable Plots of all Variables within a Metadata-containing dataset.
#' @export
#'
#' @examples batch_png_plot(qtc, 'Classification', 5, toLog10 = T, remove_outliers = 'sd', plot_type = 'std_boxplot', contrasts = list(c('A', 'B'), c('A', 'C')), test_type = 't.test', folder_name = 'Dataset)
batch_png_plot <- function(qtc, colName_toSortBy, num_info_col, toLog10 = F, remove_outliers = 'sd',
                           palette = c('#FF61F5','#FF3587','#FF3131', '#FFA56A', '#FFC101', '#81d100',
                                       '#00a827', '#00d692', '#00c9c9', '#1495ff', '#1434ff','#7900db'),
                           palette_fill = c('aquamarine', 'deeppink', 'cyan1', 'purple',
                                            'sienna1', 'seagreen2', 'cornflowerblue', 'red3'),
                           plot_type = 'std_boxplot', colName_toFill = NULL,
                           contrasts = NULL, test_type = 't.test', y_label = 'Counts',
                           height = 8, width = 11, units = 'in', res = 300,
                           title_size = 36, y_axis_text_size = 11, x_axis_text_size = 18, y_axis_title_size = 14,
                           legend_position = 'none',
                           folder_name = 'plotsOutputDir'){
  
  require(tidyverse)
  
  # Name cleaning
  names(qtc) <- str_replace_all(names(qtc), '/', ' Or ') ## Mainly for semi-quantitative (perhaps qualitative as well) metabolomics where biomarkers are sometimes indicated with a '/'
  names(qtc) <- str_replace_all(names(qtc), ':', '.')
  
  dir.create(folder_name)
  
  for (i in seq((num_info_col+1), ncol(qtc))){
    
    ### ------------------ Regular x1 per JPEG
      png(filename = paste0('./', folder_name, '/', colnames(qtc)[i], '.png'), width = width, height = height, units = units, res = res)
      p1 <- plot_graph(qtc, colName_toSortBy, i, toLog10 = toLog10,
                       remove_outliers = remove_outliers, palette = palette, palette_fill = palette_fill,
                       plot_type = plot_type, colName_toFill = colName_toFill,
                       contrasts = contrasts, test_type = test_type,
                       y_label = y_label,
                       title_size = title_size, y_axis_text_size = y_axis_text_size,
                       x_axis_text_size = x_axis_text_size, y_axis_title_size = y_axis_title_size,
                       legend_position = legend_position)
      print(p1)
    
    # ### ------------------ Only if you want to print x4 on 1 JPEG? Maybe find a way to graph out the pathways as specified
    # if (you_want_x4){
    #  require(cowplot) 
    #  jpeg(file = paste(colnames(qtc)[i], (i-num_info_col)/4+1, '.jpg', sep = ''), width = 11, height = 8, units = 'in')
    #   plot_list <- vector('list')
    #   if (i+3 > ncol(qtc)){
    #     for (j in seq(i, ncol(qtc))){
    #       p1 <- plot_graph(qtc, col_name, j)
    #       plot_list[[j-i+1]] <- as_grob(p1)
    #     }
    #     print(plot_grid(plotlist = plot_list, ncol = 2, nrow = 2))
    #   }else{
    #     for (j in seq(i,i+3)){
    #       p1 <- plot_graph(qtc, col_name, j)
    #       plot_list[[j-i+1]] <- as_grob(p1)
    #     }
    #     print(plot_grid(plotlist = plot_list, ncol = 2, nrow = 2))
    #   }
    # }
    dev.off()
  }
}

# You MUST correctly indicate whether the data should be logged. The volcano output follows -lg(P-value) against log2(FC).
#' Plots a Volcano Plot of -lg(P-Value) against log2(Fold Change).
#'
#' @param df Dataframe. Dataframe to be used for plotting.
#' @param pval.threshold Integer. The threshold for the -lg(P-Value) value. Default is 1.3.
#' @param fc.threshold.left Integer. The negative threshold for the log2(Fold Change) value. Default is -1.
#' @param fc.threshold.right Integer. The positive threshold for the log2(Fold Change) value. Default is 1.
#' @param volc.title String. The title of the Volcano Plot. Default is blank ('').
#' @param text_size Integer. The text size of the labelled metabolites fulfilling the mentioned thresholds. Default is 2.5.
#' @param toLog2FC Boolean. Do you want to log2 the FC column?
#' @param toNegLog10Pval Boolean. Do you want to -log10 the Pval column?
#' @param x.title String. The x-axis title of the Volcano Plot. Default is 'log2(Fold Change)'.
#' @param y.title String. The y-axis title of the Volcano Plot. Default is '-lg(P-Value)'.
#' @param x_limit Numeric Vector. The x-limit window to set for the volcano plot. toLog2FC must be F. Default is NULL.
#' @param repel_labels Boolean. Do you want to repel the labels for better organization? geom_text_repel will be used. Default is F.
#' @param margins Numeric Vector of 4. To create margins around the graphs to include labels. c(y.max, x.max, y.min, x.min). Default is a vector of NAs.
#'
#' @return A ggplot2 object of a Volcano Plot of -lg(P-Value) against log2(Fold Change).
#' @export
#'
#' @examples plot_volcano(df, volc.title = 'Group A vs Group B')
plot_volcano <- function (df, pval.threshold = 1.3, fc.threshold.left = -1, fc.threshold.right = 1,
             volc.title = '', x.title = 'log2(Fold Change)', y.title = '-lg(P-Value)',
             text_size = 2.5, toLog2FC = F, toNegLog10Pval= F, x_limit = NULL, repel_labels = F,
             margins = NULL){
  require(tidyverse)
  require(ggrepel)
  ## Checks for foldchange & p-value configuration. Consider all cases. X must be FC. Y must be Pval.
  ## Untransformed Pval lies from 0-1. FC usually hovers around 1, can be more or less than 1 but always positive.
  ## Transformed Pval lies from 0 to infinity. FC can be negative now and hovers around 0.
  if(toNegLog10Pval){ # Entirely Raw / Untransformed Pval => can filter for 0 < Pval < 1; let's hope your FC data ain't < 1 only.
    if(all((unlist(df[,2]) > 0) & (unlist(df[,2] <= 1)))){df <- df[,c(1,3,2)]}
  }else if(!toLog2FC){ # All Transformed => log2(FC) has both +ve and -ve, while -log(p-value) must always be positive
    if(any(unlist(df[,3] < 0))){df <- df[,c(1,3,2)]}
  }else{ # Transformed Pval, non-transformed FC - we'll use distribution mean comparisons. FC should have mean closest to 1.
    absmean1 <- abs(mean(unlist(df[,2]), na.rm = T) - 1)
    absmean2 <- abs(mean(unlist(df[,3]), na.rm = T) - 1)
    if (absmean2 < absmean1){df <- df[,c(1,3,2)]} #This means the df[,3] contains the FC, which should not be the case.
  }
  
  ### -- Perform Transformations--- (The columns should be in the right positions)
  if(toLog2FC){df[,2] <- log(unlist(df[,2]), base = 2)}
  if(toNegLog10Pval){df[,3] <- log(unlist(df[,3]), base = 10) * -1}
    
  names(df) <- c('Name', 'X', 'Y') ### X is FC, Y is -log(p-value)... well it should be!
  
  ### -- Initializing color scheme --
  df$color<- 'black'
  df$color[which(df$Y > pval.threshold)] <- 'green4'
  df$color[which(df$X > fc.threshold.right & df$Y > pval.threshold)] <- 'red'
  df$color[which(df$X < fc.threshold.left & df$Y > pval.threshold)] <- 'blue'
  
  df <- df[order(df$color),]
  
  ## Only provide labels to those which are colored to avoid overcluttering
  df$labels <- df$Name
  df$labels[which(df$color == 'black')] <- ''
  
  ## Setting xlim and ylim if they exceed certain thresholds (which skews the graphs badly) - only if non-Logged FC
  if(!is.null(x_limit)){
    if(!any(is.na(x_limit))){x_limit <- sort(x_limit)}
    if(!toLog2FC){
      if(!is.na(x_limit[1])){
        if(min(df$X) > x_limit[1]){x_limit[1] <- NA} # To not set an absurdly large window
      }
      if(!is.na(x_limit[2])){
        if(max(df$X) < x_limit[2]){x_limit[2] <- NA}
      }
    }
    if(all(is.na(x_limit))){rm(x_limit)}
  }
  
  p1 <- ggplot2::ggplot(df, aes(x = X, y = Y, color = color)) + geom_point()
  
  if(repel_labels){ ### Repelling or not repelling
    p1 <- (p1 + coord_cartesian(clip = 'off')
          + geom_text_repel(aes(fontface = 'bold'),
                            label = unlist(df$labels),
                            size = text_size, alpha = 0.7, vjust = 1.3,
                            force = 1.2, hjust = 0.5,
                            segment.size = 0.6, segment.alpha = 0.5, max.overlaps = Inf)) ## unleash the labels
  }else{
    p1 <- p1 + geom_text(aes(fontface = 'bold'),
                         label = unlist(df$labels),
                         size = text_size, alpha = 0.7, vjust = 1.3,)
  }
  
  p1 <- (p1 + geom_hline(yintercept = pval.threshold, linetype = 'dotted', color = 'gray1', size = 1) ## pval limit
         + geom_vline(xintercept = fc.threshold.right, linetype = 'dotted', color = 'gray1', size = 1) ## fc limit right
         + geom_vline(xintercept = fc.threshold.left, linetype = 'dotted', color = 'gray1', size = 1) ## fc limit left
         + geom_vline(xintercept = 0, linetype = 'solid', color = 'gray3', size = 1) ## 0 Point 
         + xlab(x.title) + ylab(y.title) + ggtitle(volc.title) ## Labelling
         + theme_minimal() + scale_color_manual(values = c(unique(df$color)))
         + theme(plot.title = element_text(face = 'bold', size = 16, hjust = 0.5),
                 axis.title.x = element_text(face = 'bold', size = 12),
                 axis.title.y = element_text(face = 'bold', size = 12),
                 legend.position = 'none')) ## Font size/type configuration
  
  if(!is.null(margins)){p1 <- p1 + theme(plot.margin = unit(margins, 'in'))}
  if(!is.null(x_limit) & exists('x_limit')){p1 <- p1 + xlim(x_limit)}
  
  return(p1)
}

## File must contain metabolite names, -lg(P-Values), log2(Fold Change) values.
## If the plot looks weird please look at your data again. Does it need to be transformed? Is it already Transformed?

#' Reads in a file of .txt of .csv and outputs the Corresponding Volcano Plot in a New Directory. (Builds on plot_volcano)
#'
#' @param filename String. Filename of .txt/.csv file to read in. Either fill this or dataframe in. Default is NULL.
#' @param dataframe Dataframe. Dataframe to use should we choose not to read in a file. Default is NULL.
#' @param pval.threshold Integer. The threshold for the -lg(P-Value) value. Default is 1.3. 
#' @param fc.threshold.left Integer. The negative threshold for the log2(Fold Change) value. Default is -1.
#' @param fc.threshold.right Integer. The positive threshold for the log2(Fold Change) value. Default is 1.
#' @param volc.title String. The title of the Volcano Plot. Default is blank ('').
#' @param text_size Integer. The text size of the labelled metabolites fulfilling the mentioned thresholds. Default is 2.5.
#' @param height Integer. Height of plot picture output in specified units. Default is 8 (in).
#' @param width Integer. Width of plot picture output in specified units. Default is 11 (in).
#' @param units String. c('in', 'px', 'cm', 'mm') Specified units for plot picture dimensions. Default is 'in'.
#' @param res Integer. Resolution in ppi recorded in bitmap file. Default is 300.
#' @param toLog2FC Boolean. Do you want to log2 the FC column?
#' @param toNegLog10Pval Boolean. Do you want to -log10 the Pval column?
#' @param x.title String. The x-axis title of the Volcano Plot. Default is 'log2(Fold Change)'.
#' @param y.title String. The y-axis title of the Volcano Plot. Default is '-lg(P-Value)'.
#' @param x_limit Numeric Vector. The x-limit window to set for the volcano plot. toLog2FC must be F. Default is NULL.
#' @param repel_labels Boolean. Do you want to repel the labels for better organization? geom_text_repel will be used. Default is F.
#' @param margins Numeric Vector of 4. To create margins around the graphs to include labels. c(y.max, x.max, y.min, x.min). Default is a vector of NAs.
#'
#' @return Outputs a ggplot2 Volcano Plot in a new directory ('Volcano Plot Outputs')
#' @export
#'
#' @examples read_plot_volcano(df_volcano.txt, volc.title = 'Group A vs Group B')
read_plot_volcano <- function(filename = NULL, dataframe = NULL, pval.threshold = 1.3, fc.threshold.left = -1, fc.threshold.right = 1,
                              toLog2FC = F, toNegLog10Pval= F, x_limit = NULL, repel_labels = F, margins = NULL,
                              volc.title = '', x.title = 'log2(Fold Change)', y.title = '-lg(P-Value)',
                              text_size = 2.5, height = 8, width = 11, units = 'in', res = 300){
  if(!is.null(filename)){
    extension <- substr(filename, (nchar(filename)-3), nchar(filename))
    if(extension == '.txt'){
      df <- read.delim(filename, as.is = T, check.names = F) ### Note that the data should have 3 columns: Name, X and Y values
    }else if(extension == '.csv'){
      df <- read.csv(filename, as.is = T, check.names = F)
    }else{
      stop('This function does not support this filetype :(, please convert it to .txt or .csv filetype instead.')
    }
  }else if(!is.null(dataframe)){
    df <- dataframe
  }else{
    stop('You indecisive cupcake! D: Do you want to read in a file? Or use a dataframe? Only one of the can be filled and the other must be NULL.')
  }
  
  ### -- Check for filetype --

  if(!file.exists('Volcano Plot Outputs')){dir.create('Volcano Plot Outputs')}
  
  png(filename = paste0('Volcano Plot Outputs/VolcanoPlot_', volc.title, '.png'), width = width, height = height, units = units, res = res)
  p1 <- plot_volcano(df, pval.threshold = pval.threshold, fc.threshold.left = fc.threshold.left, fc.threshold.right = fc.threshold.right,
                     toLog2FC = toLog2FC, toNegLog10Pval = toNegLog10Pval, x_limit = x_limit,
                     repel_labels = repel_labels, margins = margins,
                     volc.title = volc.title, x.title = x.title, y.title = y.title, text_size = text_size)
  print(p1)
  dev.off()
}


#' Plots Single-Variable Line Graphs Across Time (of Different Categories)
#'
#' @param qtc Dataframe. Dataframe containing the values of the Variables (should be single/median values), their categorical classifications and time (metadata).
#' @param colName_toSortBy String. Column name containing the different categories of data to apply contrasts on.
#' @param xAxisName String. Column name for the X-axis. Generally should be related to Time analysis - indicating the different time points.
#' @param mtbl_name String. Column name of variable to be plotted in that graph.
#' @param y_label String. Y-label for graph. Default is 'Log10(Counts)'.
#' @param x_label String. X-label for graph. Default is 'Time'.
#'
#' @return A ggplot2 plot of the line graph across time.
#' @export
#'
#' @examples plot_graph_line(qtc, 'Classification', 'Time in Minutes', 'Biomarker 33')
plot_graph_line <- function(qtc, colName_toSortBy, xAxisName, mtbl_name, y_label = 'Log10(Counts)', x_label = 'Time'){
  require(tidyverse)
  
  mtbl_name <- paste0('`', mtbl_name, '`')  ## We cannot eval-parse spaces or slashes!
  
  p1 <- ggplot(data = qtc, aes(x = get(xAxisName), y = eval(parse(text = mtbl_name)),
                               group = get(colName_toSortBy))) +
    theme_minimal() +
    ggtitle(str_remove_all(mtbl_name, '`')) + theme(plot.title = element_text(size = 30, hjust = 0.5, face = 'bold'),
                                               # axis.title.x = element_text(size = 14, face = 'bold', vjust = -1),
                                               axis.text.x = element_text(size = 11, face = 'bold'),
                                               axis.title.x = element_text(size = 18, face = 'bold'),
                                               axis.title.y = element_text(size = 18, face = 'bold'),
                                               axis.text.y = element_text(size = 11, face = 'bold'),
                                               legend.position = 'right') +
    ylab(y_label) + xlab(x_label) + 
    geom_line(aes(linetype = Classification, color = Classification), lwd = 1.5) + 
    geom_point(aes(color = Classification, shape = Classification), size = 9) + 
    scale_color_manual(values = c('aquamarine3', 'deeppink', 'cyan3', 'purple')) + ## Note that manual inputs MUST always be in order of aes specification
    scale_shape_manual(values = c(15, 16, 17, 18))
  
  return(p1)
}



### ==== Utility Scripts ====

#' Integrates Relevant Ratios into the MRMKit quant_table Output
#'
#' @param my_df Dataframe. quant_table by MRMkit
#' @param ratio_filename String. .csv file containing the ratios in a particular format (follow the format)
#'
#' @return quant_table with the relevant ratios calculated and appended to the end of the quant_table
#' @export
#'
#' @examples add_ratios(quant_table, 'our_ratios.csv')
add_ratios <- function(my_df, ratio_filename){ ## Ensure you have a ratio_file of suitable columns (to specify)
  ratio_list <- read.csv(ratio_filename, as.is = T, check.names = F)
  for (i in 1:nrow(ratio_list)){
    biomarker1 <- ratio_list$biomarker1[i]
    biomarker2 <- ratio_list$biomarker2[i]
    ratio_name <- ratio_list$ratio_name[i]
    
    if(!(biomarker1 %in% colnames(my_df)) | !(biomarker2 %in% colnames(my_df))){next()}
    
    if(ratio_list$Comments[i] == '' | ratio_list$Comments[i] == 'ET'){  ## The ET is a problem, I will have to move it out the comments column soon.
      my_df[,ratio_name] <- my_df[,biomarker1] / my_df[,biomarker2]
    }else if (ratio_list$Comments[i] == 'Special'){
      biomarker3 <- ratio_list$biomarker3[i]
      my_df[,ratio_name] <- my_df[,biomarker1] / (my_df[,biomarker2] + my_df[,biomarker3]) ## This only applies to GABR
    }
  }
  return(my_df)
}


#' Remove Pooled Samples from the quant_table MRMkit output.
#'
#' @param qt Dataframe. quant_table from MRMkit output to remove pooled samples from.
#' @param sorting_column String. The column to filter through and select the pooled samples.
#' @param pool_id String. The ID/name which indicates the presence of a pooled sample.
#'
#' @return A dataframe identical to qt with all pooled samples removed.
#' @export
#'
#' @examples remove_pool(qt, sorting_column = 'Classification', pool_id = 'pool')
remove_pool <- function(qt, sorting_column = 'type', pool_id = 'BQC'){
  return(qt[-which(unlist(qt[sorting_column]) == pool_id),])
}

#' Performs Log Transformation on Dataframe Variable Data
#'
#' @param qtc Dataframe. The dataframe to be log-transformed.
#' @param num_info_col Integer. The number of metadata columns in qtc. Non-variables not to be logged.
#' @param base Integer. Log-Base to perform. Default is 10.
#'
#' @return A dataframe that has undergone log-transformation according to the specified base.
#' @export
#'
#' @examples log_df(qtc, num_info_col = 10, base = 10)
log_df <- function(qtc, num_info_col, base = 10){
  qtc[qtc == 0] <- NA
  qtc_log <- cbind(qtc[,c(1:num_info_col)], log(qtc[,c((num_info_col+1):ncol(qtc))], base = base))
  return(qtc_log)
}

### ==== Statistics ====
#' Generates Dataframe of P-values from Specified Contrasts in a Dataset
#'
#' @param qtc Dataframe. The dataframe for reference to perform statistical tests on.
#' @param colName_toSortBy String. Column name containing the different categories of data to apply contrasts on.
#' @param num_info_col Integer. Number of metadata (non-variable) columns to ignore when performing statistical tests.
#' @param mycontrasts List. List of contrasts to perform statistical tests on (eg. list(c('A', 'B'), c('B', 'C'))). Default is 'all', where all contrasts possible in the column of colName_toSortBy will be considered.
#' @param test_type String. Statistical test to perform. c('t.test', 'wilcox.test') currently. Default is 't.test'
#' @param paired Boolean. Do you want paired statistical tests? Default is F.
#' @param colName_SortPairs String. Column name to order data for paired statistical tests. Ensure that this is also within the front info columns. Default is NULL.
#' @param grepl_fixed Boolean. Contrast filtering uses grepl, do you want fixed = T? (Matches phrase exactly in grepl).
#'
#' @return A dataframe containing the p-values of the statistical tests between all specified contrasts and of all variables within the dataframe.
#' @export
#'
#' @examples a_stat_test(qtc, 'Classification', 6, mycontrasts = list(c('A', 'B'), c('A', 'C')), test_type = 'wilcox.test')
a_stat_test <- function(qtc, colName_toSortBy, num_info_col, mycontrasts = 'all', grepl_fixed = F,
                        test_type = 't.test', paired = F, colName_SortPairs = NULL){
  final_table <- data.frame(marker_name = colnames(qtc)[(num_info_col+1):ncol(qtc)])
  
  if(class(mycontrasts) != 'list'){
    if(mycontrasts == 'all'){
      categories <- unique(unlist(qtc[colName_toSortBy]))
      mycontrasts <- list()
      for(i in 1:(length(categories)-1)){
        j <- i+1
        while(j <= length(categories)){
          mycontrasts <- append(mycontrasts, list(c(categories[i], categories[j])))
          j <- j+1
        }
      }
    }else{
      warning('That is a ridiculous hat you wear! Please input a list of contrasts for mycontrasts or leave as default for all possible contrasts.')
      break
    }
  }
  
  for(contrast in mycontrasts){
    stop = F
    cat1 <- contrast[1]
    cat2 <- contrast[2]
    contrast_name <- paste0(cat1, ' vs ', cat2, ' Pvals')
    print(str_remove_all(contrast_name, '\\\\b'))
    
    qtc1 <- qtc[which(grepl(cat1, unlist(qtc[colName_toSortBy]), fixed = grepl_fixed)),]
    qtc2 <- qtc[which(grepl(cat2, unlist(qtc[colName_toSortBy]), fixed = grepl_fixed)),]
    
    if(!is.null(colName_SortPairs)){
      qtc1 <- qtc1[order(unlist(qtc1[colName_SortPairs])),]
      qtc2 <- qtc2[order(unlist(qtc2[colName_SortPairs])),]
      if(any(is.na(qtc1)) | any(is.na(qtc2))){
        warning('Your data has NA values. Statistical Pairings may not work as intended.. D:')
      }
    }
    
    
    pval_container <- vector()
    for(i in (num_info_col + 1):ncol(qtc)){
      if(test_type == 't.test'){
        pval <- t.test(as.numeric(unlist(qtc1[,i])), as.numeric(unlist(qtc2[,i])), paired = paired)$p.value
        pval_container <- pval_container %>% append(pval)
      }else if(test_type == 'wilcox.test'){
        pval <- wilcox.test(as.numeric(unlist(qtc1[,i])), as.numeric(unlist(qtc2[,i])), paired = paired)$p.value
        pval_container <- pval_container %>% append(pval)
      }else{
        warning('Unknown test type.')
        stop = T
        break
      }
    }
    
    if(stop){break}
    
    final_table <- final_table %>% cbind(pval_container)
    names(final_table)[which(names(final_table) == 'pval_container')] <- contrast_name
  }
  
  return(final_table)
}

#' Calculates Correlation between Target Metabolites and All Other Metabolites
#'
#' @param qtc Dataframe. Your quant_table of interest.
#' @param num_info_col Integer. Number of information columns/metadata in quant_table.
#' @param corr_element String. Name of Metabolites to test correlations against.
#' @param method String. c('pearson', 'spearman', 'kendall'). Method of correlation. Default is 'pearson'.
#'
#' @return Dataframe containing the correlation coefficients and p-values of all indicated metabolites.
#' @export
#'
#' @examples behold_correlations(qtc, 8, 'Metabolites A', method = 'kendall')
behold_correlations <- function(qtc, num_info_col, corr_element, method = 'pearson'){
  require(stringr)
  master_cor <- data.frame(mtbls = colnames(qtc)[(num_info_col+1):ncol(qtc)])
  coef <- pvals <- vector()
  for(i in (num_info_col+1):ncol(qtc)){
    other_element <- unlist(qtc[,i])
    cor.object <- cor.test(unlist(qtc[corr_element]), other_element, method = method)
    coef <- coef %>% append(as.numeric(cor.object$estimate))
    pvals <- pvals %>% append(as.numeric(cor.object$p.value))
  }
  master_cor <- master_cor %>% cbind(coef, pvals)
  names(master_cor) <- c('Metabolites', paste0(stringr::str_to_title(method), '_Correlation'), paste0('P_Values'))
  return(master_cor)
}

#' Calculates Fold Change of Given Contrasts from quant_table Dataframe
#'
#' @param qtc Dataframe. The quant_table to be used.
#' @param colName_toSortBy String. Column name containing the different categories of data to apply contrasts on.
#' @param num_info_col Integer. Number of information columns/metadata in quant_table.
#' @param mycontrasts List. List of contrasts to perform statistical tests on (eg. list(c('A', 'B'), c('B', 'C'))). Default is 'all', where all contrasts possible in the column of colName_toSortBy will be considered.
#' @param toLog2 Boolean. Do you want to Log2 the FC calculated?
#'
#' @return A dataframe containing the fold changes between all specified contrasts and of all variables within the dataframe.
#' @export
#'
#' @examples
calculate_fc <- function(qtc, colName_toSortBy, num_info_col, mycontrasts = 'all', toLog2 = F){
  fc_list <- data.frame(Metabolites = colnames(qtc)[(num_info_col+1):ncol(qtc)])
  
  if(class(mycontrasts) != 'list'){
    if(mycontrasts == 'all'){
      categories <- unique(unlist(qtc[colName_toSortBy]))
      mycontrasts <- list()
      for(i in 1:(length(categories)-1)){
        j <- i+1
        while(j <= length(categories)){
          mycontrasts <- append(mycontrasts, list(c(categories[i], categories[j])))
          j <- j+1
        }
      }
    }else{
      warning('That is a ridiculous hat you wear! Please input a list of contrasts for mycontrasts or leave as default for all possible contrasts.')
      break
    }
  }
  
  for(contrast in mycontrasts){
    df1 <- qtc[which(grepl(contrast[1], unlist(qtc[colName_toSortBy]))),]
    df2 <- qtc[which(grepl(contrast[2], unlist(qtc[colName_toSortBy]))),]
    comp_name <- paste0(contrast[2], '/', contrast[1])
    print(comp_name)
    
    fc_list_baby <- vector()
    for(i in (num_info_col+1):ncol(qtc)){
      fc <- mean(unlist(df2[,i]), na.rm = T) / mean(unlist(df1[,i]), na.rm = T)
      fc_list_baby <- fc_list_baby %>% append(fc)
    }
    
    if(toLog2){
      fc_list_baby <- log(fc_list_baby, base = 2)
    }
    
    fc_list <- fc_list %>% cbind(fc_list_baby)
    names(fc_list)[ncol(fc_list)] <- comp_name
  }
  
  return(fc_list)
}


### ==== QC Scripts ====

#' Filters out Biomarkers with Bad CoV and D-Ratio Values.
#'
#' @param qt Dataframe. The raw quant_table output from MRMkit to filter out biomarkers from.
#' @param num_info_col Integer. The number of metadata columns (non-variable).
#' @param cov_row Integer. Which row contains the CoV values to filter? Default is 7.
#' @param dratio_row Integer. Which row contains the D-Ratio values to filter? Default is 9.
#' @param cov_threshold Integer. The threshold to filter CoV from, accepting biomarkers with CoV values less than or equals to this value. Default is 30.
#' @param dratio_threshold Integer. Identical to cov_threshold but for the D-Ratio instead. Default is 50.
#'
#' @return A dataframe with biomarkers having BOTH undesirable CoV and D-Ratio values filtered out.
#' @export
#'
#' @examples cov_dratio_filter(qt, 3)
cov_dratio_filter <- function(qt, num_info_col, cov_row = 7, dratio_row = 9,
                              cov_threshold = 30, dratio_threshold = 50){
  cov_list <- qt[cov_row,-c(1:num_info_col)]
  cov_less30 <- colnames(cov_list)[which(as.numeric(unlist(cov_list)) <= cov_threshold)]
  dratio_list <- qt[dratio_row, -c(1:num_info_col)]
  dratio_less50 <- colnames(dratio_list)[which(as.numeric(unlist(dratio_list)) <= dratio_threshold)]
  
  survived_biomarkers <- union(cov_less30, dratio_less50)
  bad_biomarkers <- colnames(qt)[-c(1:num_info_col, which(colnames(qt) %in% survived_biomarkers))]
  bad_biomarkers <- bad_biomarkers[order(bad_biomarkers)] # Generates a list of rejected biomarkers for reference.
  
  qtc <- qt[-c(1:10),] # Okay generally quant_table output should have 10 rows of extra info (CoV, D-Ratio etc.)
  qtc <- qtc[,c(1:num_info_col, which(colnames(qtc) %in% survived_biomarkers))]
  row.names(qtc) <- NULL
  
  cat(paste0('Unique Analyte List: ', length(cov_list), '\nSurvived CoV < 30: ', length(cov_less30), '\nSurvived D-Ratio < 50: ',
               length(dratio_less50), '\nTotal Survived: ', length(survived_biomarkers),
               '\nTotal Eliminated: ', length(bad_biomarkers)))
  
  return(list(qtc, bad_biomarkers))
}

### == COV/D-Ratio Comparison (pre-post correction) ==
## Requires 2x quant_tables (pre and post correction)

#' Initializes the Dataframe required for QC Comparison Plotting from MRMKit quant_table output.
#'
#' @param pre_qtc_filename String. Filename for the pre-corrected MRMKit-generated quant_table.
#' @param post_qtc_filename String. Filename for the post-corrected MRMKit-generated quant_table.
#' @param pre_num_info_col Integer. Number of info-columns/metadata in the pre-corrected MRMKit-generated quant_table.
#' @param post_num_info_col Integer. See above. Applies to post-corrected MRMKit-generated quant_table.
#'
#' @return A dataframe containing the list of CoVs and D-Ratios for each metabolite, for either pre-corrected or post-corrected.
#' @export
#'
#' @examples initialize_cov_dratio_lists('quant_table(pre).txt', 'quant_table(post).txt')
initialize_cov_dratio_lists <- function(pre_qtc_filename, post_qtc_filename, pre_num_info_col = 3, post_num_info_col = 3){
  
  ## Pre and Post qtcs should ideally be the RAW output from MRMKit
  bad_qtc <- read.delim(pre_qtc_filename, as.is = T, check.names = F)
  qt <- read.delim(post_qtc_filename, as.is = T, check.names = F)
  
  bad_cov_list <- data.frame(t(bad_qtc[7, -c(1:pre_num_info_col)]))
  bad_dratio_list <- data.frame(t(bad_qtc[9, -c(1:pre_num_info_col)]))
  bad_list <- cbind(bad_cov_list, bad_dratio_list)
  bad_list$Classification <- 'Before Correction'
  
  qt_cov_list <- data.frame(t(qt[7,-c(1:post_num_info_col)]))
  qt_dratio_list <- data.frame(t(qt[9, -c(1:post_num_info_col)]))
  good_list <- cbind(qt_cov_list, qt_dratio_list)
  good_list$Classification <- 'After Correction'
  
  names(bad_list) <- names(good_list) <- c('CoV', 'D-Ratio', 'Classification')
  
  cd_df <- rbind(bad_list, good_list)
  return(cd_df)
}

## Requires initialize_cov_dratio_lists to use

#' Creates a Directory and Outputs Histograms comparing the Distribution of CoV and D-Ratio Pre and Post Correction by MRMkit
#'
#' @param pre_qtc_filename String. Filename for the pre-corrected MRMKit-generated quant_table.
#' @param post_qtc_filename String. Filename for the post-corrected MRMKit-generated quant_table.
#' @param pre_num_info_col Integer. Number of info-columns/metadata in the pre-corrected MRMKit-generated quant_table.
#' @param post_num_info_col Integer. See above. Applies to post-corrected MRMKit-generated quant_table.
#' @param height Integer. Height of plot picture output in specified units. Default is 8 (in).
#' @param width Integer. Width of plot picture output in specified units. Default is 11 (in).
#' @param units String. c('in', 'px', 'cm', 'mm') Specified units for plot picture dimensions. Default is 'in'.
#' @param res Integer. Resolution in ppi recorded in bitmap file. Default is 300.
#'
#' @return Creates Histograms Comparing CoV and D-Ratio Distributions Pre vs Post Correction.
#' @export
#'
#' @examples plot_cov_dratio_comparisons('quant_table(pre).txt', 'quant_table(post).txt')
plot_cov_dratio_comparisons <- function(pre_qtc_filename, post_qtc_filename, pre_num_info_col = 3, post_num_info_col = 3,
                                        width = 11, height = 8, units = 'in', res = 300){
  
  cd_df <- initialize_cov_dratio_lists(pre_qtc_filename = pre_qtc_filename,
                                       post_qtc_filename = post_qtc_filename,
                                       pre_num_info_col = pre_num_info_col,
                                       post_num_info_col = post_num_info_col)
  
  bad_cov_mean <- round(mean(cd_df$CoV[which(cd_df$Classification == 'Before Correction')]), 1)
  good_cov_mean <- round(mean(cd_df$CoV[which(cd_df$Classification == 'After Correction')], na.rm = T), 1)
  
  ### Please change your labels accordingly...
  
  cov_hist <- ggplot(cd_df, aes(x = CoV, fill = Classification)) + theme_minimal() +
    geom_histogram(binwidth = 5, color = 'black', alpha = 0.4, position = 'identity') +
    scale_fill_manual(values = c('blue', 'red')) +
    geom_vline(xintercept = 30, linetype = 'dashed', color = 'black', size = 1) + # The Cutoff of 30%
    geom_text(mapping = aes(x = 30, y = 60, label = '30%', hjust = -0.3, fontface = 'italic'), size = 6) +
    geom_vline(xintercept = good_cov_mean, linetype = 'dotted', color = 'blue', size = 1) +
    geom_vline(xintercept = bad_cov_mean, linetype = 'dotted', color = 'red', size = 1) +
    
    ylab('Cumulative Frequency') + xlab('CoV Value') +
    ggtitle(paste0('Histogram of CoV Comparison;\nMean of CoV: ', bad_cov_mean, ' --> ', good_cov_mean)) +
    theme(plot.title = element_text(size = 20, hjust = 0.5, face = 'italic'),
          axis.text.x = element_text(size = 18),
          axis.title.x = element_text(size = 20, face = 'bold'),
          axis.title.y = element_text(size = 20, face = 'bold'),
          axis.text.y = element_text(size = 18),
          legend.position = 'right')
  
  bad_dratio_mean <- round(mean(cd_df$`D-Ratio`[which(cd_df$Classification == 'Before Correction')]), 1)
  good_dratio_mean <- round(mean(cd_df$`D-Ratio`[which(cd_df$Classification == 'After Correction')], na.rm = T), 1)
  
  dratio_hist <- ggplot(cd_df, aes(x = `D-Ratio`, fill = Classification)) + theme_minimal() +
    geom_histogram(binwidth = 5, color = 'black', alpha = 0.4, position = 'identity') +
    scale_fill_manual(values = c('blue', 'red')) +
    geom_vline(xintercept = 50, linetype = 'dashed', color = 'black', size = 1) + # The Cutoff of 50%
    geom_text(mapping = aes(x = 50, y = 60, label = '50%', hjust = -0.3, fontface = 'italic'), size = 6) +
    geom_vline(xintercept = good_dratio_mean, linetype = 'dotted', color = 'blue', size = 1) +
    geom_vline(xintercept = bad_dratio_mean, linetype = 'dotted', color = 'red', size = 1) +
    
    ylab('Cumulative Frequency') + xlab('D-Ratio Value') +
    ggtitle(paste0('Histogram of D-Ratio Comparison;\nMean of D-Ratio: ', bad_dratio_mean, ' --> ', good_dratio_mean)) +
    theme(plot.title = element_text(size = 20, hjust = 0.5, face = 'italic'),
          axis.text.x = element_text(size = 18),
          axis.title.x = element_text(size = 20, face = 'bold'),
          axis.title.y = element_text(size = 20, face = 'bold'),
          axis.text.y = element_text(size = 18),
          legend.position = 'right') + 
    xlim(0, 200) # D-Ratio appears to have astronomical outliers
  
  if(!file.exists('CoV & D-Ratio Comparison Graphs')){dir.create('CoV & D-Ratio Comparison Graphs')}
  
  png(filename = 'CoV & D-Ratio Comparison Graphs/CoV_Graph.png', width = width, height = height, units = units, res = res)
  print(cov_hist)
  print('CoV!')
  dev.off()
  
  png(filename = 'CoV & D-Ratio Comparison Graphs/D-Ratio_Graph.png', width = width, height = height, units = units, res = res)
  print(dratio_hist)
  print('D-Ratio!')
  dev.off()
}


#' Plots Histogram of Mann Whitney U-Test P-Values (Test between Samples from Different Columns in the quant_table)
#'
#' @param filename String. The filename of the quant_table to be used.
#' @param num_info_col Integer. The number of metadata (non-variable) information columns to be ignored.
#' @param study_codename String. Reference to the naming of .mzML files - what is the codeword for the study?
#' @param num_info_rows Integer. Number of metric rows that MRMkit outputs in quant_table. Default is 10.
#' @param c1_odd Boolean. Is Column 1 on an odd number? Default is T.
#' @param sort_colname String. The column name to sort the MS columns of samples by. Default is 'filename'.
#'
#' @return A ggplot2 object of the histogram of MWUT P-Values between samples of different columns.
#' @export
#'
#' @examples plot_columnDiff_mwutVal('quant_table.txt', 3, study_codename = 'ABC', sort_colname = 'filename')
plot_columnDiff_mwutVal <- function(filename, num_info_col, study_codename, num_info_rows = 10, c1_odd = T, sort_colname = 'filename'){
  qt <- read.delim(filename, as.is = T, check.names = F)
  qtc <- qt[-c(1:num_info_rows),]
  mtbl <- colnames(qtc)[c((num_info_col + 1):ncol(qtc))]
  
  if(sort_colname == 'filename'){
    qtc$Sorting <- qtc[sort_colname] %>% unlist() %>% str_remove(study_codename) %>% str_remove('[.]mzML')
  }else{ # It should be just 'batch' or something with just numbers indicating the batch-column split of 1,2,1,2
    qtc$Sorting <- qtc[sort_colname] %>% unlist()  
  }
  
  mwut_pvals <- vector()
  print('Calculating MWUT P-Vals...')
  for (i in (num_info_col + 1):(ncol(qtc)-1)){ # We do not want the Sorting Column
    to_analyze <- cbind(qtc$Sorting, qtc[,i]) %>% apply(MARGIN = c(1,2), as.numeric) %>% as.data.frame()
    names(to_analyze) <- c('batch', 'values')
    if(c1_odd){
      c1 <- to_analyze[which(to_analyze$batch %% 2 == 1),]
      c2 <- to_analyze[which(to_analyze$batch %% 2 == 0),]
    }else{ # Why would this even be the case
      c1 <- to_analyze[which(to_analyze$batch %% 2 == 0),]
      c2 <- to_analyze[which(to_analyze$batch %% 2 == 1),]
    }
    mwut_pvals <- mwut_pvals %>% append(wilcox.test(unlist(c1$values), unlist(c2$values))$p.value)
  }
  
  print('Done!')
  mwut_df <- data.frame(mtbl, mwut_pvals)
  
  affected_mtbls <- length(which(mwut_df$mwut_pvals < 0.05))
  mwut_histogram <- ggplot(mwut_df, aes(x = mwut_pvals)) + theme_minimal() +
    geom_histogram(binwidth = 0.01, color = 'black', fill = 'pink') +
    geom_vline(xintercept = 0.05, linetype = 'dashed', color = 'red') +
    geom_text(mapping = aes(x = 0.05, y = affected_mtbls*2/3, label = '0.05', hjust = -0.3, color ='red', fontface = 'italic'), size = 6) +
    ylab('Cumulative Frequency') + xlab('Mann-Whitney Test P-Value') +
    ggtitle(paste0('Histogram for Mann-Whitney U-Test;\nNo. of metabolites affected: ', affected_mtbls)) +
    theme(plot.title = element_text(size = 20, hjust = 0.5, face = 'italic'),
          axis.text.x = element_text(size = 18),
          axis.title.x = element_text(size = 20, face = 'bold'),
          axis.title.y = element_text(size = 20, face = 'bold'),
          axis.text.y = element_text(size = 18),
          legend.position = 'none')
  
  return(mwut_histogram)
}

#' Outputs a Plot into a New Directory
#'
#' @param p1 ggplot2 plot. The plot object to be output in the specified directory.
#' @param foldername String. Directory name to output the plots in. If the directory is present, the plot will be output directly into it. Default is 'plotsdir'.
#' @param filename String. The filename of the plot to be output. Default is 'myplot'.
#' @param filetype String. c('png', 'jpg', 'pdf'). What is the output file type? Default is 'png'.
#' @param height Integer. Height of plot picture output in specified units. Default is 8 (in).
#' @param width Integer. Width of plot picture output in specified units. Default is 11 (in).
#' @param units String. c('in', 'px', 'cm', 'mm') Specified units for plot picture dimensions. Default is 'in'.
#' @param res Integer. Resolution in ppi recorded in bitmap file. Default is 300.
#'
#' @return Plot will be printed as a file (of specified filetype) in the specified directory.
#' @export
#'
#' @examples output_plot(p1, foldername = 'myfolder', filename = 'myplot', filetype = 'png')
output_plot <- function(p1, foldername = 'plotsDir', filename = 'myplot', filetype = 'png', 
                        height = 8, width = 11, units = 'in', res = 300){
  
  if(!file.exists(foldername)){dir.create(foldername)}
  
  if(filetype == 'png'){
    png(filename = paste0(foldername, '/', filename, '.png'), width = width, height = height, units = units, res = res)
    print(p1)
    dev.off()
  }else if(grepl('jpg|jpeg', filetype)){
    jpeg(filename = paste0(foldername, '/', filename, '.jpg'), width = width, height = height, units = units, res = res)
    print(p1)
    dev.off()
  }else if (filetype == 'pdf'){
    pdf(file = paste0(foldername, '/', filename, '.pdf'), width = width, height = height)
    print(p1)
    dev.off()
  }else{
    warning('What filetype is that? Please specify either png, jpg or pdf')
  }
}

### ==== WKL Parser and batchinfo Generator ====
#' Inputs a dataframe directly imported from a .wkl file and returns a dataframe with filename, column_number and type of sample.
#'
#' @param wkl Dataframe. This is the imported wkl file in .wkl format, fresh from the MassHunter!
#'
#' @return data.frame object matching file names to column types
#' @export
#'
#' @examples
#' wkl_parser(wkl) (It can't get simpler than this)
wkl_parser <- function(wkl){
  
  require (tidyverse)
  
  names(wkl) <- 'Yay'
  wkl <- wkl[which(grepl('AcqMethod|DataFileName', wkl$Yay) | grepl('<Name>', wkl$Yay)),]
  
  ### Assuming the first 4 lines are always present (wkl[c(1:4),])
  wkl <- wkl[-c(1:4)]
  col_type <- wkl[which(grepl('.m<', wkl))]
  data_file <- wkl[which(grepl('.d<', wkl))]
  sample_pool <- wkl[which(grepl('<Name>', wkl))]
  
  for (i in 1:length(col_type)){
    col <- col_type[i]
    col_start <- str_locate(col, '.m</')[1]-2
    col <- substr(col, col_start+1, col_start + 1)
    col_type[i] <- col
    
    rm(col, col_start, i)
  }
  
  for (j in 1:length(data_file)){
    filename <- data_file[j]
    file_start <- as.numeric(regexpr("\\\\[^\\]*$", filename)) + 1
    file_end <- as.numeric(gregexpr('.d</', filename)) + 1
    filename <- filename %>% substr(file_start, file_end)
    data_file[j] <- filename
    
    rm(filename, file_start, file_end, j)
  }
  
  for (k in 1:length(sample_pool)){
    type <- sample_pool[k]
    type <- type %>% str_remove('.*?>') %>% str_remove('<.*')
    if(grepl('blank|pool', type) | grepl('N/A', type)){
      sample_pool[k] <- type
    }else{
      sample_pool[k] <- 'Sample'
    }
  }
  
  df <- data.frame(data_file, col_type, sample_pool)
  names(df) <- c('Filename', 'ColumnType', 'Type')
  
  return(df)
}


#' Reads all .wkl files from MassHunter in 'Input WKL Files' folder and outputs the batchinfo file used as the input for MRMKit.
#' Please ensure that you have your files in the 'Input WKL Files' folder
#'
#' @param cohort_name Just how you want to name your file. 
#'
#' @return A .txt file consolidating the columns used for .d file in a cohort.
#' @export
#'
#' @examples
#' batchinfo_generator('Cohort1') Yay.
batchinfo_generator <- function(cohort_name = 'Cohort'){
  
  require (tidyverse)
  
  master_df <- mastered_df <- data.frame(Filename = character(), ColumnType = numeric(),
                                         Type = character(), Batch = numeric())
  
  ## Reading the input files. Please ensure that you have the files in 'Input WKL Files' folder
  batch <- 1
  for (filename in list.files('./Input WKL Files')){
    wkl <- read.delim(paste0('./Input WKL Files/', filename), as.is = T, check.names = F)
    df <- wkl_parser(wkl)
    
    df$Batch <- as.numeric(df$ColumnType) + (batch-1)*2
    
    master_df <- master_df %>% rbind(df)
    
    print(batch)
    batch <- batch + 1
    
    rm(df, wkl, filename)
  }
  
  rm(batch)
  
  ## Sometimes Dr Leroy puts 'WorklistData', you can append to this list if required
  if(any(grepl('WorklistData', master_df$Filename))){
    master_df <- master_df[-which(grepl('WorklistData', master_df$Filename)),]
  }
  
  ## Some ordering and cleaning (Separate by length of characters in filename - differing by number of digits of file number [eg. 10 vs 100 vs 1000])
  nametypes <- table(nchar(master_df$Filename))
  if(length(nametypes) > 1){
    for(i in 1:length(nametypes)){
      length_filename <- names(nametypes)[i]
      df_current <- master_df[which(nchar(master_df$Filename) == length_filename),]
      df_current <- df_current[order(df_current$Filename),]
      mastered_df <- mastered_df %>% rbind(df_current)
      rm(length_filename, df_current, i)
    }
  }
  
  master_df <- mastered_df
  rm(mastered_df, nametypes)
  
  if(any(duplicated(master_df$Filename))){ # Wait there shouldn't be. :O
    master_df <- master_df[-which(duplicated(master_df$Filename)),]
    warning('There are duplicate names D:. Check the wkl files for dupes. I have removed the duplicates using duplicated() but it may still be kind of iffy.')
  }
  
  ### Final Arrangement and Cleaning
  master_df <- master_df[,c(1, 4, 3)] 
  names(master_df) <- c('file', 'batch', 'type')
  
  master_df <- master_df[-which(grepl('blank|N/A', master_df$type)),] # Removes Blanks
  master_df <- master_df[order(master_df$batch),] # Orders the Batches
  master_df$type[which(grepl('pool', master_df$type))] <- 'BQC' # Names 'pool' to 'BQC'
  master_df$file <- master_df$file %>% str_replace_all('.d', '.mzML') # .d to .mzML in filename
  
  if(!file.exists('BatchInfo (Column-Corrected)')){
    dir.create('BatchInfo (Column-Corrected)')
  }  
  
  ## Printing
  write.table(master_df, paste0('./BatchInfo (Column-Corrected)/batchinfo_', cohort_name, '.txt'), row.names = F, sep = '\t', quote = F)
}
