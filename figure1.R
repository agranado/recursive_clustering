# May 22nd 2021
# choose parameters for this pathway
# all functions return a grob plot object

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
average_profileMatrix <-function(control_res =list(), rank_df = data.frame() ,
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
                     rename(label=class_label)-> diverse_mat


    # 5. label profiles as disperse of private
    #    left_join with p_value profiles.
    diverse_mat  %>%
    left_join(rank_df %>% select(label,p_value,mean_silh,p_value_private),  by ='label')  -> diverse_mat


    # type I: disperse cell cell_types
    # profiles that express cell types with higher dispersion than 0.9 quantile of the
    # transcriptome + significan p_value
    diverse_mat$profile_class ='profile'

    control_res$diversity %>%
        dplyr::filter(type=='transcriptome') %>% # choose the null model
        pull(d) %>% quantile(0.90) -> divers_quantile

    diverse_mat$adj_pval = p.adjust(diverse_mat$p_value, method ='BH')
    # filter for high dispersion profiles
    diverse_mat %>% mutate(profile_class =ifelse(diversity>divers_quantile & mean_silh >0 & cell_types> 2 & adj_pval <=sig_level, 'disperse',profile_class)) -> diverse_mat

    #diverse_mat_ann %>% dplyr::filter(diversity>divers_quantile & p_value <= sig_level) -> diverse_mat
    # type II: private
    # we filter for profiles that express in similar cell types
    # here we are not looking for profiles with -more similar cell types- than the transcriptomes
    # since the transcriptome imposes a lower bound in theory
    # here we are looking for profiles with diversity scores within the expected distribution
    # by default less than the median + significant p_value
    control_res$diversity %>%
        dplyr::filter(type=='transcriptome') %>% # choose the null model
        pull(d) %>% quantile(0.75) -> median_dispersion # profile falls within the distribution for similar cell types

    diverse_mat$adj_pval_private = p.adjust(diverse_mat$p_value_private, method ='BH')

    diverse_mat %>% mutate(profile_class = ifelse(diversity<=median_dispersion & cell_types>2 & adj_pval_private <= sig_level, 'private',profile_class)) -> diverse_mat


    return(diverse_mat %>% as.data.frame)


}

# May 20th 2021
# type: disperse or private
# diversE_mat: average profiles matrix with labels for disperse and private
profileHeatmap <- function(diverse_mat=data.frame(),
                            pathway_name = '',
                            control_res = list(),type='disperse',
                            save_dir = 'plots/figuresApril/fig4/',
                            sig_level = 0.05, plot_organ_heatmap = F , filter_tissues = c() ){


    # annotate with p-values to filter only significant profiles
    # diverse_mat %>% tibble::rownames_to_column(var='label') %>%
    # left_join(rank_df %>% select(label,p_value,mean_silh,p_value_private) %>%
    #               mutate(label=as.character(label)),
    #           by ='label')  -> diverse_mat_ann


    diverse_mat %>% dplyr::filter(profile_class==type) -> diverse_mat
    diverse_mat %>% as.data.frame -> diverse_mat
    row.names(diverse_mat) <- diverse_mat$label
    # we define the rows order by clustering independently so we can
    # keep the same ordering for the organ heatmap
    diverse_mat[,genesPathway(pathway_name)] %>% pheatmap(clustering_method = 'ward.D2', silent=T) -> p_clust


    motif_heatmap <- superheat(diverse_mat[,genesPathway(pathway_name)],
                               #pretty.order.rows = T,
                               order.rows = p_clust$tree_row$order ,
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
      organDistribution_heatmap(control_res, diverse_mat, pathway_name, filter_tissues, save_dir, row_order =  rev(p_clust$tree_row$order))
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
                                      save_dir = 'plots/figuresApril/fig4/',
                                      row_order = c() ) {

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
  pheatmap(sqrt(x[row_order, ]), treeheight_row = 20,treeheight_col = 20,
  		clustering_method = 'ward.D2',col = tissue_pal(100),
  		fontsize =12,angle_col = 45,cluster_rows = F,
  		filename = paste(save_dir, pathway_name, '_Tissue_distribution_bmp.pdf'),
  		height =4, width = 6) # keep this size to make it similar to the motif heatmap


}

# we can generate the palette once and use it for all plots
pathway_palette<- function(pathway_df = data.frame() ){
  # generate a palette for all profiles
  dot_colors = makeQualitativePal(length(unique(pathway_df$class_label)),
                                  rand_order = T, skip = skip_colors, tail_colors = F )
}

# umap distribution of diverse and private profiles
# plot all profiles together
umapDistribution <- function(pathway_df = data.frame(), # pathway clustering
                              dot_colors = c() # color array for each profile
                            ){

      # global embedding from seurat object
      umap_coords = Embeddings(master_seurat, 'umap') %>% as.data.frame
      umap_coords$cell_id = row.names(umap_coords)
      # assign profile class to those cell types expressing the pathway
      umap_coords %>% left_join(pathway_df %>% select(cell_id,class_label),by='cell_id') %>%
          replace_na(replace = list(class_label = 0)) -> umap_df

      dot_colors[1] ="#D3D3D3"
      umap_df$class_label <- as.character(umap_df$class_label)

      umap_df %>%  ggplot(aes(x = UMAP_1,y = UMAP_2,color = class_label)) +
          geom_point(size = 1) + scale_color_manual(values = dot_colors) +
          theme_minimal() + theme(text = element_text(size = 10)) +
          coord_cartesian(xlim= c(-5,7),ylim = c(-10,5))
}

# mini umap array
# plot multiple profiles in the umap projection as a single grid. one umap per profile
umapGridPlot <- function(diverse_mat = data.frame(), # data.frame with dispersion scores, p_values, etc
                          pathway_df = data.frame(), #pathway clustering
                            sig_level = 0.05, save_pdf = F,
                            user_embedding = c() ){
    # conditions for significant profiles
    # 1. silhouette score >0
    # 2. more than 2 cell types
    # 3. p_value (adjusted) < 0.05
    # 4. among the highest (lowest) dispersion
    # diverse_profiles <- rank_df %>% dplyr::filter(mean_silh >0 & n> 2 & adj_pval <=sig_level ) %>%
    #                   arrange( desc(diversity)) %>% top_n(n =5, wt = diversity) %>% pull(label )
    #
    # specific_profiles <- rank_df %>% dplyr::filter(mean_silh>0 & n>2 & p_value_private <=sig_level ) %>% # choose the null model
    #     top_n(n =5, wt = -diversity) %>% pull(label )


    diverse_profiles <- diverse_mat %>% dplyr::filter(profile_class=='disperse') %>%
                          mutate(total_diversity = diversity*cell_types) %>%
                          top_n(n =5, wt = total_diversity) %>% pull(label )
    specific_profiles <- diverse_mat %>% dplyr::filter(profile_class=='private') %>%
                          mutate(total_diversity = diversity*cell_types) %>%
                          top_n(n =5, wt = -diversity) %>% pull(label )


    plot_profiles = c(diverse_profiles, specific_profiles)

    # global embedding from seurat object
    if(length(user_embedding)==0){
      umap_coords = Embeddings(master_seurat, 'umap') %>% as.data.frame
    }else{
      umap_coords = user_embedding
    }
    umap_coords$cell_id = row.names(umap_coords)
    # assign profile class to those cell types expressing the pathway
    umap_coords %>% left_join(pathway_df %>% select(cell_id,class_label),by='cell_id') %>%
        replace_na(replace = list(class_label = 0)) -> umap_df



    mini_umap_list = list()
    for(i in 1:length(plot_profiles)){
        which_motif = plot_profiles[i]


        umap_df %>% dplyr::filter(class_label == which_motif ) -> motif_df

        umap_df %>%  ggplot(aes(x = UMAP_1,y = UMAP_2)) +
            geom_point(size = 0.3, color = "#D3D3D3")  +
            theme_minimal() +
            geom_point(data = motif_df, aes(x =UMAP_1, y= UMAP_2), color ='red', size = 0.5) +
            theme_pubr(base_size =8 ) + ggtitle(paste('Profile', which_motif, 'n=',dim(motif_df)[1])) +
            coord_cartesian(xlim= c(-5,5),ylim = c(-9,5)) +
            theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                  axis.ticks.y = element_blank(), axis.text.y  = element_blank(),line = element_blank(),
                  axis.title.y=element_blank(), axis.title.x=element_blank()) -> g_umap


        mini_umap_list[[i]] <- g_umap

    }

    if(!save_pdf){
      gg = grid.arrange(grobs= mini_umap_list,ncol = length(diverse_profiles))
    }else{
      gg = arrangeGrob(grobs= mini_umap_list,ncol = length(diverse_profiles))
    }

    return(gg)
}

# plot a dendrogram including all expressing cell types
# global metric can be computed directly from gene expression values of from a PCA embedding
globalDendrogramPlot <-function(control_res = list(), # main results object
                                seurat_obj = c() ,
                                hvg_genes = c(), # list of genes to create the global dendrogram
                                dist_metric ='cosine', # clustering distance in global space. applies for expression and PCA spaces
                                use_pca = F, # wether to use PCA coordinates instead of gene expression for global dendrogram
                                npcs = 100, # default for all seurat objects
                                save_dir = 'plots/figuresApril/fig4/'){

  pathway_df = control_res$profiles
  row.names(pathway_df) <- pathway_df$cell_id

  if(!use_pca){
    # use highly variable genes and compute distance directly on gene expression (log scale)
    hvg_df = makeMainDataFrame(hvg_genes, master_seurat = seurat_obj)
    hvg_genes = hvg_genes[hvg_genes %in% colnames(hvg_df )]
    # filter only for expresing cell types in this pathway
    hvg_df %>% dplyr::filter(cell_id %in% pathway_df$cell_id ) -> hvg_df
    row.names(hvg_df) <-hvg_df$cell_id

    x = hvg_df[, hvg_genes]
    x_scaled = scale(x ) # standardScaler
    row.names(x_scaled) <- hvg_df$cell_id
  }else{
    x_scaled = Embeddings(seurat_obj, reduction='pca')
    x_scaled = x_scaled[pathway_df$cell_id,1:npcs] # already set row names as cell id

  }
  # use cosine distance of euclidean. applies for gene expression and PCA spaces
  if(dist_metric=='cosine') clust_dist_metric=dist.cosine(x_scaled) else clust_dist_metric = dist(x_scaled)

  k_motifs = length(pathway_df$class_label %>% unique ) # pathway states are already clustered
  # colors for pathway classes
  colors_1206$class_label <- makeQualitativePal(k_motifs, glasbey_use = T, skip  =1) # skip white color
  names(colors_1206$class_label) <- 1:k_motifs

  pathway_df$class_label = pathway_df$class_label  %>% as.character()
  # Make heatmap
  pheatmap(x_scaled, annotation_row = pathway_df %>% select( Cell_class,age, class_label),
  				annotation_colors = colors_1206, show_colnames = F, show_rownames = F,
  				clustering_distance_rows = clust_dist_metric , treeheight_col = 0,
  				cutree_rows = 12 , fontsize = 10,
  				filename = paste(save_dir ,pathway_name,'global_dendrogram.pdf',sep=""),
  				height = 5, width =4)
}
