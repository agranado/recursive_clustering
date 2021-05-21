
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
  		ylab('Fraction of clusters') + xlab('Cell state dispersion') +
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

# latest version of this functino. takes the rank data.frame
# works better with the integrated pipeline.
diversity_plot2 <-function(rank_df = data.frame(), title_main='',control_res = list(),
                            upper_quantile=0.9){
  title_main=which_pathway

  global_df = control_res$diversity %>% dplyr::filter(type=='transcriptome') %>% rename(diversity=d)

  rank_df %>% ggplot(aes(x=rank,y=diversity,size=n *0.3)) + geom_point(alpha=0.4) +
      geom_hline(yintercept=quantile(global_df$diversity, upper_quantile),linetype="dashed", color = "lightgray", size = 1) +
      geom_hline(yintercept=quantile(global_df$diversity, 1-upper_quantile),linetype="dashed", color = "lightgray", size = 1) +
      ylab("Cell state dispersion (PCA)")  + theme_pubr(base_size =10 )  +
      xlab("Pathway profile")  + ggtitle(title_main ) +
      annotate(geom="text", x=5, y=quantile(global_df$diversity, upper_quantile) +2 , label="Expected",
               color="black")
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

# Find profiles in disperse cell states
# first we get the average profile matrix
# returns a matrix of average profiles grouped by their profile class
average_profileMatrix <-function(control_res =list(),
                                pathway_name=''){

      # 1. New: Select all profiles
    diverse_df <- control_res$profiles
    # 2. tidy data.frame to average gene expression within a pathway profile
    diverse_df %>% pivot_longer(cols = genesPathway(pathway_name),
                                names_to = 'gene', values_to = 'expression') %>%
        select(cell_id, cell_ontology_class, Tissue,
               Cell_class, gene, expression, class_label, rank, diversity, n) -> diverse_tidy
    # 3. average profile
    diverse_tidy %>% group_by(class_label,gene) %>%
        summarise(mean_expr = mean(expression),
                  rank = mean(rank),diversity = mean(diversity),
                  cell_types = mean(n)) -> diverse_tidy

    # 4. wide format
    diverse_tidy %>% pivot_wider(id_cols = c(class_label,rank, diversity, cell_types),
                                 names_from=gene,values_from=mean_expr) %>%
                     tibble::column_to_rownames('class_label')-> diverse_mat


    return(diverse_mat)


}

# May 20th 2021
# type: disperse or private
#
profileHeatmap <- function(diverse_mat=matrix(),
                            pathway_name = '',
                            control_res = list(),type='disperse',
                            save_dir = 'plots/figuresApril/fig4/',
                            sig_level = 0.05, plot_organ_heatmap = F , filter_tissues = c() ){


    # annotate with p-values to filter only significant profiles
    diverse_mat %>% tibble::rownames_to_column(var='label') %>%
    left_join(rank_df %>% select(label,p_value,mean_silh,p_value_private) %>%
                  mutate(label=as.character(label)),
              by ='label')  -> diverse_mat_ann

    if(type=='disperse'){
      # type I: disperse cell cell_types
      # profiles that express cell types with higher dispersion than 0.9 quantile of the
      # transcriptome + significan p_value
      control_res$diversity %>%
          dplyr::filter(type=='transcriptome') %>% # choose the null model
          pull(d) %>% quantile(0.90) -> divers_quantile

      diverse_mat_ann %>% dplyr::filter(diversity>divers_quantile & p_value <= sig_level) -> diverse_mat
    }else if(type=='private'){
      # type II: private
      # we filter for profiles that express in similar cell types
      # here we are not looking for profiles with -more similar cell types- than the transcriptomes
      # since the transcriptome imposes a lower bound in theory
      # here we are looking for profiles with diversity scores within the expected distribution
      # by default less than the median + significant p_value
      control_res$diversity %>%
          dplyr::filter(type=='transcriptome') %>% # choose the null model
          pull(d) %>% quantile(0.50) -> median_dispersion

      diverse_mat_ann %>% dplyr::filter(diversity<=median_dispersion & p_value_private <= sig_level) -> diverse_mat
    }

    row.names(diverse_mat) <- diverse_mat$label
    motif_heatmap <- superheat(diverse_mat[,genesPathway(pathway_name)],
                               pretty.order.rows = T,
                               heat.pal = blues_pal(10),
                               #order.rows = order(diverse_mat$rank),

                               #heat.pal.values = seq(0,1,0.1),
                               bottom.label.text.angle = 90,
                               yr = sqrt(diverse_mat$cell_types),
                               yr.plot.type='bar',
                               yr.axis.name = "N cell types",
                               row.title = paste(type, 'profiles'),
                               column.title = "Pathway genes",
                               bottom.label.size = 0.3,
                               grid.hline.col = "white",
                               grid.hline.size = 2,

                               left.label.text.size = 3,
                               bottom.label.text.size = 3

    )

    # plot organ heatmap
    if(plot_organ_heatmap){
      organDistribution_heatmap(control_res, diverse_mat, pathway_name, filter_tissues, save_dir)
    }


    # save to pdf
    gg = motif_heatmap$plot
    # given the font size and ggsave parameters
    # saving at width = 3.5 is optimal for 11 genes
    # height = 4.5 works for 11 profiles
    fig_width = 3.5/11 * length(genesPathway(pathway_name))
    # 1.5 in baseline for plot annotations
    fig_height = 0.37*4.5 + 3.5/11 * dim(diverse_mat)[1] # 37% of total height is for legend/ the rest is for the heatmap
    ggsave(paste(save_dir,
                  pathway_name,'_',type,'_heatmap_may.pdf',sep=''),
                  gg, height=fig_height, width = fig_width)

    # return the graphic object, the averaged matrix
    return(list(motif_heatmap, diverse_mat))
}

# for disperse profiles, we want to explore
# their distribution across different tissues in the organism
# for this, we create a table of organ vs profile
# also, we want to keep the same order as in the motif heatmap
#
organDistribution_heatmap <- function(control_res = list() ,
                                      diverse_mat  = data.frame(), pathway_name = '',
                                      filter_out_tissues = c(),
                                      save_dir = 'plots/figuresApril/fig4/') {

  diverse_df = control_res$profiles # extracte the pathway clustering
  diverse_df %>% dplyr::filter(class_label %in% diverse_mat$label ) %>%
    select(Tissue, class_label ) %>% dplyr::filter(!Tissue %in% filter_out_tissues) %>% group_by(class_label,Tissue)%>%
    count  %>%
    pivot_wider(id_cols=class_label,names_from = Tissue,
                values_from= n,values_fill = 0) -> tissue_distribution


  # color palette for the organ heatmap
  tissue_pal<-colorRampPalette(brewer.pal(n = 9, name = 'PuRd'))

  x = tissue_distribution %>% ungroup %>% select(-class_label) %>% as.matrix()
  row.names(x) <- tissue_distribution$class_label
  # make tissue names pretty
  colnames(x)<-str_replace(string = colnames(x),pattern ="_",replacement = ' ') %>% firstup()

  # use pheatmap
  # we want to plot the profiles in the same order as diverse mat heatmap
  #
  pheatmap(sqrt(x), treeheight_row = 20,treeheight_col = 20,
  		clustering_method = 'ward.D2',col = tissue_pal(100),
  		fontsize =12,angle_col = 45,
  		filename = paste(save_dir, pathway_name, '_Tissue_distribution_bmp.pdf'),
  		height =4, width = 6) # keep this size to make it similar to the motif heatmap


}
