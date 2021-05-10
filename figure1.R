
# choose parameters for this pathway
# all functions return a grob plot object
# global parameters
min_expr_threshold = 0.35
pathway_index = 2
# p is the pathway index
pathway_name =  param_list[[pathway_index]][[1]]
min_genes_pathway = param_list[[pathway_index]][[2]]
# min_index for the min_expr value
silhoutte_max <- function(p = pathway_index,
                          min_expression = min_expr_threshold,
                          x_offset = 1, pathway_name = 'Bmp' ,
                          min.y = 0.1, max.y = 0.5){
    #p = 2 # which pathway
    # index x[[1]] means the ggplot object
    min_index = which(min_expr_list==min_expression)

    x = all_res_silh[[p]][[min_index]] # 4 here is index for the optimal min.expr
    # Extract data.frames from list and remove offset on the x-axis
    boot_df = x[[2]] %>% dplyr::filter(k>x_offset)
    boot_df_control = x[[3]] %>% dplyr::filter(k>x_offset)

    # Plot together with shaded area
    boot_df %>% ggplot(aes(x = k , y = m)) +
        geom_ribbon(aes(ymin = m -s, ymax = m+s) ,fill = 'lightblue', alpha = 0.2) +
        geom_line(color ='blue') + geom_ribbon(data = boot_df_control, aes(ymin = m -s, ymax = m+s), alpha = 0.2) +
        geom_line(data=boot_df_control)  + theme_pubr(base_size = 10 ) + ylab('Silhouette') +
        xlab('Number clusters') + ggtitle(pathway_name) + scale_y_continuous(limits=c(min.y, max.y))-> g
    # Compute z-score from randomized data
    diff = x[[2]]$m-x[[3]]$m;
    z_score = diff/x[[3]]$s
    max_z = max(z_score, na.rm = T)
    z_score = z_score / max_z
    df = data.frame(k=1:100, z = z_score)
    df %>% dplyr::filter(k>x_offset) -> df

    # Plot z-score as second plot on the bottom
    z_score_plot <- axis_canvas(g, axis = 'x') +
        geom_line(data = df, aes(x = k, y = z), color = 'red') +
        geom_vline(xintercept = which.max(z_score))
    combined_plot <- insert_xaxis_grob(g,z_score_plot, position = 'bottom ')


    #grid.arrange(grobs = list(combined_plot))
    return(combined_plot)
  }


# Run main pipeline with diversity controls

# choose optimal k from silhouette plot
optimal_k_pathway = 20
control_res  = fullControlPathway(this_pathway =pathway_name,
                                  k_pathway = optimal_k_pathway , filter_pathway = 'both',
                                  this_seurat = master_seurat,
                                  null_list = hvg_genes,
                                  n_samples = 100, filter_manual = T,
                                  min_expr_gene = min_expr_threshold,min_genes_ON = min_genes_pathway)

# this function takes the result from the full pipeline
# plot the distribution of diversity across all pathway clusters
# two controls show the upper and lower bounds on diversity epected from similar and random cell types.
violin_diversity <- function(control_res = list(), p_value = T,plot_null = F){
    my_comparisons = list( c('pathway','transcriptome'),c('pathway','pos control') )

    if(!plot_null){
      df = control_res$diversity %>% dplyr::filter(type %in% c('pathway','transcriptome','pos control'))

      lines_palette = c("#00AFBB", "#a8a8a8", "#000000") # color order: pathway, transcriptome, random
    }else{
      df = control_res$diversity %>% dplyr::filter(type %in% c('pathway','transcriptome','pos control','null dist'))
      lines_palette = c("#00AFBB", "#a8a8a8", "#000000",'#6f7e96') # color order: pathway, transcriptome, random
    }

    violin_plot <- ggviolin(df, x = 'type',y='d', fill = 'type',
                            palette = lines_palette,
                            add = 'boxplot',
                            add.params = list(fill ='white'),
                            font.label = list(size = 10, size = "plain"))
    if(p_value){
       violin_plot +  stat_compare_means(comparisons = my_comparisons) +
        stat_compare_means(label.y = control_res$diversity$d %>% max ) -> violin_plot
      }else{
        violin_plot + scale_y_continuous(trans='sqrt') +scale_x_discrete(labels = NULL) + coord_flip() -> violin_plot
      }

    violin_plot + ylab('Cell type diversity') +
        xlab('Cell type clustering') + theme_pubr(base_size = 10 ) -> violin_plot

        return(violin_plot)
}

# same as above but in ECDF format
ecdf_diversity <- function(control_res = list() ,plot_null=F, pathway_name = 'Bmp'){


  if(!plot_null){
    df = control_res$diversity %>% dplyr::filter(type %in% c('pathway','transcriptome','pos control'))

    lines_palette = c("#00AFBB", "#000000", "#a8a8a8") # color order: pathway, transcriptome, random
  }else{
    df = control_res$diversity %>% dplyr::filter(type %in% c('pathway','transcriptome','pos control','null dist'))
    lines_palette = c("#00AFBB", "#000000", "#a8a8a8",'#6f7e96') # color order: pathway, transcriptome, random
  }
  quantile_line = quantile(df %>% dplyr::filter(type=='transcriptome') %>% pull(d), 0.9)

    df  %>% ggplot(aes(x = d, col= type)) + stat_ecdf()  + ggtitle(pathway_name) +
      scale_x_continuous(trans='sqrt') +theme_pubr(base_size = 10) +
  		scale_color_manual(values =lines_palette) +
  		ylab('Fraction of clusters') + xlab('Cell type diversity') +
      geom_vline(xintercept=quantile_line,linetype="dashed", color = "lightgray", size = 1)-> ecdf_plot
}

# plot all profiles ranked by cell type diversity
diversity_plot <- function(stats_list = list() , title_main = 'Notch'){


		path_stats = stats_list$pathway
		global_stats = stats_list$global
    upper_quantile = 0.75

		g <-	path_stats %>% ggplot(aes(x = rank, y = diversity, size = round(n*0.5) )) + geom_point()  +
			    geom_hline(yintercept=quantile(global_stats$diversity, upper_quantile),linetype="dashed", color = "lightgray", size = 1) +
			    geom_hline(yintercept=quantile(global_stats$diversity, 1-upper_quantile),linetype="dashed", color = "lightgray", size = 1) +
			    ylab("Cell type diversity (PCA)")  + theme_pubr(base_size =10 )  +
			    xlab("Pathway profile")  + ggtitle(title_main ) +
			    annotate(geom="text", x=3, y=quantile(global_stats$diversity, upper_quantile) +2 , label="Expected",
			             color="black")

     return(g)
}

rank_diversity <- function(which_pathway ='Bmp', k_motifs = 20, min_expression = 0.2,
                          min_genes_pathway = 2, embedding_matrix = c() ,
                          global_dist_metric ='euclidean' , make_plot = T ){
    # takes the pathway name as input:
    # all parameters come from upstream functions
    res_list = quickPipeline(which_pathway = which_pathway ,
                        master_seurat = master_seurat,
                        k_final = k_motifs , min_expr = min_expression, min_genes_on = min_genes_pathway,
                        which_profiles = 'both')


    stats_list = global_clusterig(df_devel = res_list$df_devel,
                                master_seurat = master_seurat,
                                k_final = k_motifs,
                                manual_embedding = embedding_matrix ,
                                dist_metric =global_dist_metric,
                                pathway_genes = genesPathway(which_pathway ) )

    if(make_plot){
      return(diversity_plot(stats_list, which_pathway))
    }else{
      return(stats_list)
    }
}
