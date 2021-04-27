## functions Mar 2021
## Clustering and silhouette score statistics.
## Applies for any Seurat object

library(gridExtra)
library(stylo)
library(cluster)

#Apr 2021
# Quick pipeline functions
# Filters non-expressing cell types
# Filters also for devel or adult profiles


genesPathway <- function(which_pathway = 'Notch'){

  this_pathway = all_pathways %>% dplyr::filter(pathway ==which_pathway ) %>% pull(gene)
  this_pathway = this_pathway[which(this_pathway %in% row.names(master_seurat))]

  return(this_pathway)
}

# 1. create devel/adult data.frame
# 2. Filter cell types with at least 2 expressed genes
# 3. Cluster cosine + ward.D2
# 4. returns clustered data.frame
#
quickPipeline <- function(which_pathway = 'Notch', master_seurat = c(),
                        k_final = 25,
                        min_genes_on = 1, min_expr = 0.2,
                        which_profiles = 'devel',
                        rand_control = F,
                        silent_plot = T,
                        manual_filter_cells = c(),
                        verbose = F
                         ){
        if(!rand_control){
          this_pathway = all_pathways %>% dplyr::filter(pathway ==which_pathway ) %>% pull(gene)
        }else{
          this_pathway = which_pathway # random set directly as input
        }
        this_pathway = this_pathway[which(this_pathway %in% row.names(master_seurat))]



        df_devel <- normalizedDevel(this_pathway,sat_val = 0.99,
                fill_zero_rows = F, master_seurat = master_seurat ,
                which_datasets = which_profiles)

        df_devel %>% mutate(cell_type = paste(global_cluster, '--',Tissue,': ', cell_ontology_class,'-', age, sep="")) -> df_devel
        df_devel$genes_on = rowSums(df_devel[,this_pathway]>min_expr )

        row.names(df_devel) <- df_devel$cell_type

        # two types of filtering. Either by min expression of pathway genes. Default
        # or by user-specified list of cell types. For control purposes
        if(length(manual_filter_cells)==0){
          df_devel %>% dplyr::filter(genes_on>min_genes_on ) -> df_devel
        }else{
          df_devel %>% dplyr::filter(cell_id %in% manual_filter_cells) -> df_devel
          # minority of gene sets will have non-expressing cells
          # this happens very rarely and only 1 or 2 cells will have no expression, so we can remove them
          expressing_cells =df_devel$cell_id[rowSums(df_devel[,this_pathway])>0]
          df_devel %>% dplyr::filter(cell_id %in% expressing_cells) -> df_devel
        }

        # Heatmap
        p = pheatmap(df_devel[,this_pathway ],
                     annotation_row = df_devel %>% select(dataset, Cell_class),
                     annotation_colors = colors_1206, show_rownames = T, fontsize = 5,
                     cutree_rows = k_final, clustering_method = 'ward.D2',
                      clustering_distance_rows = dist.cosine(as.matrix(df_devel[,this_pathway])),
                     cluster_cols = F, silent = silent_plot)


        ####################
        cos_labels = cutree( p$tree_row, k = k_final)
        df_devel$class_label = cos_labels


        df_devel$class_label <- as.numeric(df_devel$class_label)

        df_devel %>% gather(key = 'gene', value ='expression', this_pathway ) %>%
            group_by(class_label, gene) %>%
            summarise(mean_expr = mean(expression), n_cell_types = n()) %>%
            spread(gene, mean_expr)  %>% tibble::column_to_rownames(var  ="class_label") -> x
        # if available
        return(list('matrix' = x, 'df_devel' = df_devel) )
}

# calls the above function
# clusters the pathway genes and makes the corresponding heatmap
quickHeatmap <- function(which_pathway = 'Bmp', master_seurat = c() , filter_profiles = 'adult', k = 36, silent = F,
                            min_genes_ = 1, min_expr_gene = 0.2, save_pdf = F,
                            pdf_file = ''){
    this_pathway = genesPathway(which_pathway )
    res_list = quickPipeline(which_pathway = which_pathway,
                             master_seurat = master_seurat,
                             which_profiles = filter_profiles,
                             k_final = k, manual_filter_cells = c() ,verbose = T,
                            min_genes_on = min_genes_, min_expr = min_expr_gene)

    df_devel = res_list$df_devel
    if(!save_pdf){
    pheatmap(df_devel[,this_pathway ],
             annotation_row = df_devel %>% select(dataset, Cell_class),
             annotation_colors = colors_1206, show_rownames = F, fontsize = 10,
             cutree_rows = k, clustering_method = 'ward.D2',
             clustering_distance_rows = dist.cosine(df_devel[,this_pathway ] %>% as.matrix ), silent = silent,
            main = paste(which_pathway, '--', dim(df_devel)[1],'cell types'), col = blues_pal(100))

    }else{
      pheatmap(df_devel[,this_pathway ] ,
               annotation_row = df_devel %>% select(dataset, Cell_class),
               annotation_colors = colors_1206, show_rownames = F, fontsize = 10,
               cutree_rows = k, clustering_method = 'ward.D2',
               clustering_distance_rows = dist.cosine(df_devel[,this_pathway ] %>% as.matrix ), silent = silent,
              main = paste(which_pathway, '--', dim(df_devel)[1],'cell types'), col = blues_pal(100),
              filename = pdf_file, height = 5, width = 5)

    }

    return(df_devel)
}

#
# data is already clustered by pathway
# only for the cell types in df_devel -- those with a pathway label
# Mar 25
# manual_embedding: a user-provided latent space representing all the cells in the dataset. A distance matrix will be computed from here
global_clusterig<- function(df_devel = data.frame() , master_seurat = c() ,
														k_final = 25,n_pcs = 100 ,manual_embedding = c()){
        	# cells are filtered
        	df_devel$class_label <- df_devel$class_label %>% as.character()

        	# pathway colors
        	colors_1206$class_label = makeQualitativePal(length(df_devel$class_label %>% unique))
        	names(colors_1206$class_label) <- df_devel$class_label %>% unique

        	# PCA embedding
        	# Number of principal components
          # umap_coords by default uses cell_id as the row.name
          if(length(manual_embedding)==0){
        	umap_coords<- Embeddings(master_seurat,reduction = 'pca')
        	umap_coords <- umap_coords[,1:n_pcs]
        }else{
          # this matrix should represent a latent space with row.names equal cell_id
          # when provided directly, the number of principal components is ignored
          umap_coords = manual_embedding
        }

        	# use only those cell that have a notch profile
        	scaled_data = umap_coords[df_devel$cell_id, ]
        	row.names(scaled_data)<-df_devel$cell_type # row name is the full cell type

        	p_global = pheatmap(scaled_data %>% t ,
        				annotation_col = df_devel %>% dplyr::select( age,dataset,Tissue,Cell_class, class_label),
        				annotation_colors = colors_1206 ,
        				clustering_method = 'ward.D2',
        				cutree_cols = k_final, show_rownames = F,
        				show_colnames = F, fontsize = 12,silent= T)

                # for each pathway label, we calculate the mean distance
         # in the PCA space for the cell types in that label
         k_final = length(df_devel$class_label %>% unique )
         sapply(1:k_final, function(x){
             class_cells = df_devel$cell_id[df_devel$class_label==x]
             if(length(class_cells)>1){
             umap_coords[class_cells,] %>% dist %>% as.matrix -> dist_class

             return(dist_class[upper.tri(dist_class)] %>% mean())

             }else{
                 return(0)
             }

         }) -> diversity_pathway

        #global stats // transcriptome
        # same idea, we use the global label -----
        global_labels = cutree(p_global$tree_col, k = k_final)
        df_devel$global_label = global_labels
        # singletons have diversity = 0 by definition
        sapply(1:k_final, function(x){
            class_cells = df_devel$cell_id[df_devel$global_label==x]
            if(length(class_cells)>1){
            umap_coords[class_cells,] %>% dist %>% as.matrix -> dist_class

            return(dist_class[upper.tri(dist_class)] %>% mean())
          }else{
            return(0)
          }


        }) -> diversity_global

        # transcriptome
        df_devel %>% group_by(global_label) %>% count %>% as.data.frame() -> global_stats
        global_stats$diversity = diversity_global

        # pathway
        df_devel %>% group_by(class_label) %>% count %>% as.data.frame() %>%
        	mutate(class_label = as.numeric(class_label)) %>%
        	arrange(class_label ) -> path_stats
        path_stats$diversity = diversity_pathway


        path_stats$type = "pathway"
        global_stats$type = "transcriptome"

        path_stats %>% rename(label = class_label) -> path_stats
        global_stats %>% rename(label = global_label) -> global_stats

        # rank the profiles by diversity
        path_stats$rank = rank(path_stats$diversity)

        return(list('pathway' = path_stats, 'global' = global_stats))
}

# Plot all motifs ranked by their diversity score
# Apr 19th -- migrated from Notion
diversity_plot <-function(stats_list = list() , title_main = 'Notch'){

		# data.frames with diversity scores for each profile in the pathway
		# we use global as a control for the expected diversity in similar cell types
		#
		path_stats = stats_list$pathway
		global_stats = stats_list$global
    upper_quantile = 0.80

		g <-	path_stats %>% ggplot(aes(x = rank, y = diversity, size = n )) + geom_point()  +
			    geom_hline(yintercept=quantile(global_stats$diversity, upper_quantile),linetype="dashed", color = "orange", size = 1) +
			    geom_hline(yintercept=quantile(global_stats$diversity,.10),linetype="dashed", color = "orange", size = 1) +
			    ylab("Cell type diversity (PCA)")  + theme_pubr(base_size =12 )  +
			    xlab("Pathway profile")  + ggtitle(title_main ) +
			    annotate(geom="text", x=3, y=quantile(global_stats$diversity, upper_quantile) +2 , label="Expected",
			             color="black")

     return(g)
}

# only for the cell types in df_devel -- those with a pathway label
# Mar 25
positive_control<- function(df_devel = data.frame() , master_seurat = c() ,
														n_pcs = 100 , n_random = 100, user_embedding =c() ){
        	# cells are filtered
        	df_devel$class_label <- df_devel$class_label %>% as.character()
          k_final = length(df_devel$class_label %>% unique )
        	# PCA embedding
        	# Number of principal components
          # umap_coords by default uses cell_id as the row.name
          if(length(user_embedding)==0){
          	umap_coords<- Embeddings(master_seurat,reduction = 'pca')
          	umap_coords <- umap_coords[,1:n_pcs]
          }else{
            # Here the latent space already has the row.names
            # dimensions should be specified before
            umap_coords = user_embedding
          }

          control_list = list()
          for( i in 1:n_random){
             df_devel$class_label = sample(df_devel$class_label)
             # for each pathway label, we calculate the mean distance
             # in the PCA space for the cell types in that label
             sapply(1:k_final, function(x){
                class_cells = df_devel$cell_id[df_devel$class_label==x]
                if(length(class_cells)>1){
                umap_coords[class_cells,] %>% dist %>% as.matrix -> dist_class

                return(dist_class[upper.tri(dist_class)] %>% mean())

              }else{
                return(0)
              }

            }) -> diversity_pathway


            # pathway
            df_devel %>% group_by(class_label) %>% count %>% as.data.frame() %>%
            	mutate(class_label = as.numeric(class_label)) %>%
            	arrange(class_label ) -> path_stats
            path_stats$diversity = diversity_pathway

            path_stats$type = "pathway"
            path_stats %>% rename(label = class_label) -> path_stats
            # rank the profiles by diversity
            path_stats$rank = rank(path_stats$diversity)

            control_list[[i]] = path_stats
        }
        return(control_list )
}



# null_dist: list of genes to sample from to generate the null distribution
# which_pathway: pathway name (used for number of genes)
# which_k: optimal k found for the pathway-- test random genes using this same k
# seurat_obj: fetchData from this seurat
# n_rand: how many random sets of genes
# # filter_profile: 'devel' or 'adult' to filter cell types before running the pipeline.

# filter_profile: 'devel' or 'adult' to filter cell types before running the pipeline.
controlDiversity<-function(null_dist = c(), which_pathway = c(), which_k = c(), seurat_obj= c() ,
													 n_rand=100, filter_profile = 'adult',
                            filter_manual_cells = c() , verbose = F, n_pcs = 100, manual_embedding = c()){
	control_stats = list()


	    this_pathway = genesPathway(which_pathway)

      rand_list = list()
      for(rr in 1:n_rand)
        rand_list[[rr]] = sample(null_dist, length(this_pathway))

      # Run control pipeline in parallel for the whole list of random pathways
      # specify the manual list of cell types: overrides the filtering by min.expression
      control_stats = mclapply(rand_list, parallel_pipeline, cut_k = which_k,
                            seurat_master = seurat_obj, manual_cell_types = filter_manual_cells,
                            profile_filter = filter_profile,global_npcs = n_pcs, user_embedding = manual_embedding, mc.cores = 15)

	    # res_list = quickPipeline(which_pathway = rand_pathway,
	    #                          master_seurat = seurat_obj,
	    #                          which_profiles =filter_profile,
	    #                          k_final = which_k,
	    #                          rand_control = T)
      #
	    # stats_list = global_clusterig(df_devel = res_list$df_devel,
	    #                               master_seurat = seurat_obj,
	    #                               k_final = which_k)
	    #control_stats[[i]] = stats_list



	# make control_df
	do.call(rbind,lapply(control_stats, function(x){rbind(x$pathway %>% select(-rank), x$global)})) -> control_df

  #control_df %>% dplyr::filter(type=='pathway') -> control_df
	return(control_df)

}

# Executes negative control in parallel
# Only works in linux/mac
# manual_cell_types: user-specified list of cell types to use for clustering
# this is useful for random sets of genes, to use the same cell types as the pathway such that
# calculations of diversity are comparable.
parallel_pipeline<- function(x, cut_k = c(), seurat_master = c(), profile_filter = 'adult' ,
                            manual_cell_types = c() , global_npcs = 100, user_embedding = c()){
    # quickPipeline takes manual_cell_types as the list with cell_ids to include
    # if empty, the function will filter cell_types based on the expression of the selected genes
    res_list = quickPipeline(which_pathway = x,
                         master_seurat = seurat_master,
                         which_profiles = profile_filter,
                         k_final = cut_k, manual_filter_cells = manual_cell_types,
                         rand_control = T)

    stats_list = global_clusterig(df_devel = res_list$df_devel,
                                  master_seurat = seurat_master,
                                  k_final = cut_k, n_pcs = global_npcs, manual_embedding = user_embedding)
   return(stats_list)
}

# Mar31 2021
# runs the above functions in order for a given pathway
# Apr 14th updates:
#  We want to specify the min.expr threshold for each gene to be considered ON
#  Also, since different pathways comprise different numbers of genes, we want a fraction of the pathway to be ON as opposed to the same number for all pathways
# Update Apr 16th: Now number of PCs for global distance is a parameter to the main function.
fullControlPathway <- function(this_pathway = c() ,
                               k_pathway = c(),
                               filter_pathway = 'adult',
                               this_seurat = c(),
                               null_list = c(),
                               n_samples = 100,
                               filter_manual = F,
                               verbose_output = F,
                               min_expr_gene = 0.2,
                               min_genes_ON = 1,
                               n_pcs_global = 100,
                               embedding_matrix = c()  # Number of principal components to consider for diversity metric
                                ){

     # 1. Run the quick pipeline for the real pathway
     #    with the specified k (found previously as optimal)
     res_list = quickPipeline(which_pathway = this_pathway,
  															master_seurat = this_seurat,
  															which_profiles = filter_pathway,
  															k_final = k_pathway, verbose = verbose_output,
                                min_genes_on = min_genes_ON, min_expr = min_expr_gene)
    if(verbose_output)
      print(paste(this_pathway,": ", dim(res_list$df_devel)[1]," expressing cells" ))



    # 2. Compute diversity for the cell types containing each of the pathway profiles
     stats_list = global_clusterig(df_devel = res_list$df_devel,
                                    master_seurat = this_seurat,
                                    k_final = k_pathway, n_pcs = n_pcs_global,manual_embedding = embedding_matrix)
    #print('Done clustering real pathway')
    # 3. Negative control: random sets of genes from null list
    # Note: here we want two types of filtering:
    #       current: filter cell types based on the expression of the random set
    #       probably better: filter cell types that express the real pathway so the diversity is comparable
    # new option to filter the same cell types as the real pathway
    if(filter_manual){
      manual_cell_types = res_list$df_devel$cell_id
    }else{
      manual_cell_types = c()
    }
    # add the filter using filter_manual_cells parameter
    control_df = controlDiversity(null_dist = null_list,
    															which_pathway =this_pathway,
    															which_k = k_pathway, seurat_obj = this_seurat,
    															n_rand = n_samples, filter_profile = filter_pathway,
                                  filter_manual_cells = manual_cell_types, verbose = verbose_output,
                                  n_pcs = n_pcs_global, manual_embedding = embedding_matrix)
    #print(paste('Done clustering ',n_samples,' random sets'))

    # 4. Positive control: re-distribute pathway profiles in randomly selected cell types.
    # returns diversity scores only for the fake pathway (no global calculation)
    pos_control = positive_control(res_list$df_devel, master_seurat = this_seurat,n_random = n_samples ,n_pcs = n_pcs_global,
                                      user_embedding = embedding_matrix)
      # merge all in data.frame
    pos_control = do.call(rbind, pos_control )
    # Split negative control results: pathway (random genes), transcriptome(global diversity)
    control_df_path = control_df %>% dplyr::filter(type=='pathway')
    control_df_global = control_df %>% dplyr::filter(type=='transcriptome')

    df_diversity = rbind(data.frame(d =control_df_path$diversity, type ='null dist'), #random sets
    					 data.frame(d = control_df_global$diversity, type ='transcriptome dist'), #global dendrogram using same cells from random sets
               data.frame(d = stats_list$pathway$diversity, type = 'pathway'), # actual pathway
               data.frame(d = stats_list$global$diversity, type = 'transcriptome'), # global dendrogram using the same cells than the real pathway
               data.frame(d = pos_control$diversity, type ='pos control')) #randomized pathway profiles across cell types

    df_recurrence = rbind(data.frame(d =control_df_path$diversity*control_df_path$n, type ='null dist'),
           data.frame(d = control_df_global$diversity*control_df_global$n, type ='transcriptome dist'),
           data.frame(d = stats_list$pathway$diversity*stats_list$pathway$n, type = 'pathway'),
           data.frame(d = stats_list$global$diversity*stats_list$global$n, type = 'transcriptome'),
           data.frame(d = pos_control$diversity*pos_control$n, type ='pos control'))

    return(list(diversity= df_diversity, recurrece = df_recurrence ))
}

# Apr 1 2021
# plot histograms for the pathway diversity compared to the expected null model
fullControlPlots <- function(control_res = list(), this_pathway = c() , remove_transcriptome = F){

  # clusteringt the transcriptome for these same cells using the 100PCs
  # so far the right control is the random sets of HVGs
  # using 100PCs on the same cells tells us the expected diversity in the cells we are using
  if(remove_transcriptome)
    control_res$recurrece %>% dplyr::select(-transcriptome) -> control_res$recurrece

  control_res$recurrece %>% ggplot(aes(x = d, col= type)) + stat_ecdf() +
    scale_color_ggthemr_d() + ggtitle(this_pathway) +
    theme(text = element_text(size =20)) +
    scale_x_continuous(trans='sqrt')-> g1

    control_res$recurrece %>% dplyr::filter(type != 'transcriptome') %>%
   ggplot( aes(x = d, y = ..density.., fill = type )) +
   	geom_histogram(position = "identity", alpha = 0.2, bins = 50) +
   	theme(text = element_text(size = 20 )) + scale_x_continuous(trans='sqrt') -> g2

    return(list(g1,g2))

}


# Parameter screen plots
# test different null models
pctAbovePlot <- function(p,threshold, plot_margin = 0.1, metric_index = 'recurrece'){
  plot_colors = makeQualitativePal(5,glasbey_use = F,rand_order = T)

  # threshold = 0.75
  # p = 1

  for(n in 1:length(null_models)){
      lapply(null_res_list[[n]][[p]], function(x){
          # pct of profiles above the 0.9 quantile for the null distribution
          x[[metric_index]] %>% dplyr::filter(type=='null dist') %>% pull(d) %>% quantile(threshold) -> threshold_quant
          x[[metric_index]] %>% dplyr::filter(type=='pathway') %>% pull(d) -> pd
          sum(pd>=threshold_quant)

      }) %>% unlist -> n_above
      pct_above =n_above/k_var
      names(pct_above) = k_var
      if(n==1)
          plot(k_var, pct_above,type = "o", main = paste(test_pathways[p],'% Profiles above', threshold),
                            col = plot_colors[n],
                            ylim = c(0,(1-threshold + plot_margin) ), lwd = 1.5)

      lines(k_var, pct_above, col = plot_colors[n], lwd = 1.5)


  }
  abline(h = 1-threshold)
  #if(threshold< 0.5) y_legend = (1-threshold + plot_margin) else y_legend = plot_margin
  y_legend = (1-threshold + plot_margin)

  legend(14,y_legend,legend = c('HVG','HVG4k','Markers1','Markers2','Allgenes'),
         fill = plot_colors[1:length(null_models)],bg="transparent")

}

make_superheat<- function(x, which_pathway ='Notch'){

    this_pathway = genesPathway(which_pathway)
    # add text to mark expression values above a threshold
    text_ann = matrix("", nrow = dim(x)[1],ncol = length(this_pathway))

    text_ann[x[,this_pathway] >= 0.1] = "+"

    superheat(x[,this_pathway],
              pretty.order.rows = T,
              bottom.label.text.angle = 90,
              yr  = x$n_cell_types,
              yr.axis.name = 'N cell types',
              yr.point.size = 3,

              left.label.size = 0.4,left.label.text.size = 5,
              left.label.text.alignment = "left",
              X.text = text_ann,
              X.text.size = 4,
              X.text.col = 'white')
}

## end quick pipeline functions

###############################
## Silhouette score from scratch

# Mar 2021
# silhouette score scanning for different values of k for the full pathway list
# specify the number of random sets of genes to compute the expected distribution of silhouette scores
# Parameters: distance metric
#             clustering_method
#
silhPathwayControl <- function(x, cells_include = c(), k_max = 100 ,
                               clust_method = "complete", sat_val = 0.99,
                               dist_metric = 'euclidean', master_obj = c() ){
    this_pathway = x
    devel_adult <- normalizedDevel(this_pathway, sat_val,
                                   fill_zero_rows =T, master_seurat = master_obj)


    devel_adult %>% dplyr::filter( cell_id %in% cells_include) -> devel_adult
    row.names(devel_adult) <- devel_adult$global_cluster
    #
    x = devel_adult[,this_pathway] %>% as.matrix
    if(dist_metric =="euclidean"){
        dist_mat = dist(x)
    }else if(dist_metric=="cosine"){
        dist_mat = dist.cosine(x ) #matrix should not have zero rows
    }

    # Evaluate silhouette score for different values of k
    # cluster only once, then partition the tree using different values of k
    master_clustering = hclust(dist_mat, method =clust_method )
    ss = c()
    for(k in 2:k_max){
        clustering = cutree(master_clustering, k)
        devel_adult$motif_label = clustering

        master_clustered = devel_adult #this has motif label now
        s1<-cluster_silhouette(master_clustered, this_pathway = this_pathway,
                               dist = dist_metric ) # original clustering
        ss[k] = mean(s1$ms)
    }

    #umap_stats <-    makeUmapStats(master_clustered, s1,
    #																			this_pathway = this_pathway,
    #																			dist_method ='user',
    #																			user_dist_matrix = dist_pca)

    #umap_stats_cosine <-    makeUmapStats(master_clustered, s1,
    #																			this_pathway = this_pathway,
    #																			dist_method ='user',
    #																			user_dist_matrix = dist_cosine_pca)
    return(ss)
}

# Apr 2021
# Calculate optimal number of clusters for a given set of genes based on the silhouette score
# Performs bootstrap analysis to estimate standard deviation
silhPathwayBootstrap <- function(x, cells_include = c(), k_max = 100 ,
                               clust_method = "complete", sat_val = 0.99,
                               dist_metric = 'euclidean', master_obj = c() , pct_boots = 0.9,
                               n_boots = 100,
                               weighted_silhouette =F, control_silh = F){
    this_pathway = x
    devel_adult_main <- normalizedDevel(this_pathway, sat_val,
                                   fill_zero_rows =T, master_seurat = master_obj)

    # # Randomize columns -- destroy gene-gene correlations
    # if(control_silh)
    #   devel_adult_main = randomize_columns(devel_adult_main, this_pathway )

    boot_list = list()
    for(b in 1:n_boots){
      # Control 1: Bootstrap, sample the list of cell types
    if(!control_silh){
        # Filter out cell types that are not expressing.
        # but sample the list of cell types before to create a boostrap distribution
        cells_include_boots  = sample(cells_include, round(length(cells_include)*pct_boots) )
        devel_adult_main %>% dplyr::filter( cell_id %in% cells_include_boots) -> devel_adult
        row.names(devel_adult) <- devel_adult$global_cluster
    }else{
      # Control 2: re-shuffle the data to get a lower bound on silhouette score
      # by destroying gene-gene correlations but keeping the distribution fo each gene.
      # Before, we still need to filter for the cell types to include
      devel_adult_main %>% dplyr::filter( cell_id %in% cells_include) -> devel_adult_main
      # randomize expression within those cells! We don't want gene expression values from cell types not originally included
      devel_adult = randomize_columns(devel_adult_main, this_pathway )

    }

      #
      x = devel_adult[,this_pathway] %>% as.matrix
      if(dist_metric =="euclidean"){
          dist_mat = dist(x)
      }else if(dist_metric=="cosine"){
          dist_mat = dist.cosine(x ) #matrix should not have zero rows
      }

      # we compute the clustering only once
      main_clustering = hclust(dist_mat, method =clust_method )
      ss = c()
      for(k in 2:k_max){
          # cut the tree with different k
          clustering = cutree(main_clustering, k)
          devel_adult$motif_label = clustering

          master_clustered = devel_adult #this has motif label now
          s1<-cluster_silhouette(master_clustered, this_pathway = this_pathway,
                                 dist = dist_metric ) # original clustering
          if(!weighted_silhouette){ # normal silhouette score -- average across all clusters
            ss[k] = mean(s1$ms)
          }else{
            # A weighted average to consider the number of cell types within each cluster
            ss[k] = mean(s1$ms * s1$n)
          }
      }

      boot_list[[b]] = ss
    }

    return(boot_list )
}

# Negative control
# destroy gene-gene correlations while keeping individual distributions for each gene the same
randomize_columns<-function(df, this_pathway){
    for(i in 1:length(this_pathway))
        df[,this_pathway[i]] = sample(df[,this_pathway[i]])
    return(df )

}
# plot output of silhPathwayControl
# silhouette score for target pathway + expectation distribution as shaded area
ggplotPathway <- function(results,type = 'silh', pathway_name =''){

    pathway.data = results[[2]]
    rand.ensemble = results[[1]]


    if(type =='silh'){
        low.lim = 2
        high.lim = length(pathway.data)+1
        score.type ="Silhouette"
        legend.pos = 'none'
    }else if(type =="wss"){
        low.lim = 1
        high.lim = length(pathway.data)
        score.type = "WSS"
        legend.pos = 'none'
    }else if(type =="pca"){
        low.lim = 1
        # number of principal components
        high.lim = length(pathway.data)

        rand.ensemble = t(apply(results[[1]],1,cumsum))
        pathway.data = cumsum(pathway.data)
        score.type = "PCA"
        legend.pos = c(0.7,0.2 )
    }else if(type =='MI'){
        score.type = 'MI'

    }

    #repeats
    row.names(rand.ensemble)<-as.character(1:dim(rand.ensemble)[1])
    # K
    colnames(rand.ensemble)<-as.character(low.lim:high.lim)

    aa<-as.data.frame(rand.ensemble)
    bb= data.frame( k = low.lim:high.lim, m_score = pathway.data,sd_score = rep(0,length(pathway.data)),data = rep("Pathway",length(pathway.data)))

    control.df<-gather(aa, "k","score") %>% group_by(k) %>%
        summarise(m_score = mean(score),sd_score = sd(score)) %>% mutate(k = as.numeric(k)) %>%
        mutate(data = rep("Control",length(pathway.data)))

    control.df %>% rbind(bb) %>% ggplot(aes(x = k, y = m_score, color =data)) +
        geom_ribbon(data=control.df, aes(ymin = m_score - sd_score, ymax=m_score + sd_score),alpha = 0.2,color =NA) +
        geom_line(size = 0.5) + theme_classic() + scale_colour_manual(values=c("black", "deeppink3")) + theme(text =element_text(size =15)) +
        xlab("N clusters") + ylab(paste("Score (",score.type, ")",sep="")) +
        theme(legend.position= legend.pos) +
        theme(axis.text.x=element_text(size=10), axis.text.y = element_text(size =10)) +
        ggtitle(pathway_name)  -> p

    if(type=="wss"){
        p = p +   scale_y_continuous(trans='log10') +coord_cartesian(ylim=c(100,1000),xlim = c(2,20))
    }else if(type=="silh"){
        p = p + coord_cartesian(xlim=c(4,high.lim))
    }
    return(p)
}

# March 8th
# First function that will scan different k values
# run spectrum with those k and create a data.frame with information about silhouette and other stats
# From here we can choose the best k that will intitialize the pipeline
findOptimalK <- function(ks = c(50,60, 70,80,90,100,120), pathway = c(), top_k = 1){


  # does not run in parallel! but does not take too long
  silh_list = lapply(ks, silh_run_k, kernel ='stsc', this_pathway = pathway )

  for(i in 1:length(silh_list))
      silh_list[[i]] %>% mutate(rank = rank(desc(ms))) %>% mutate(k = ks[i] %>% as.character) -> silh_list[[i]]

  param.var.k <-do.call(rbind, silh_list)

  # choose best clustering based on number of cluster with higher silhouette
  # also maximize the number of cells that belong to good clusters.
  param.var.k %>% rename(n_cells = n) %>% dplyr::filter(ms >0.5) %>%
	   group_by(k) %>% summarise(tot_cells = sum(n_cells),n_clust  = n() ) %>%
	    mutate(k_ratio = n_clust /as.numeric(k)) %>%
	     arrange(desc(tot_cells), desc(k_ratio)) %>% top_n(top_k, wt = tot_cells) -> best_k

  return(best_k )

}
# For a given pathway run clustering step by step
# Complete pipeline
# Scan for top 5 k for initial clustering
fullSpectralRecursive <- function(this_pathway, master_seurat = c(),
                                 k_opt  = 20, min_silh = 0.5 ){# for first round ){

      # 1. Gene selection
    #this_pathway <- rand_path

    devel_adult <- makeMainDataFrame(this_pathway, master_seurat = master_seurat) #pink variables go to Shiny

    n_motifs = 29 # for final clustering of the recursive labels

    kernel_1st = 'stsc'
    kernel_2nd = 'stsc'


    p_list<-clusterPathway(
        devel_adult = devel_adult,
        which_pathway = this_pathway,
        max_satu =1, # not needed if using min.max
        min_expr = 0.0, #filte on rowsums
        k_opt = k_opt, # number of clusters Spectral
        spec_method = 1,
        spec_kernel = kernel_1st,
        unique_cell_types = F,
        min_max = T,
        sat_quantile = T, # saturate values at quantile
        sat_val = 0.99,
        min_expr_single_cell = 0.0 #min expr to consider a gene ON (after min.max)
    )

    # 3. Merge with umap coordinates and meta data
    # if no recursive clustering,then this is scatter_export
    master_clustered = makeMasterClustered(p_list, k_opt = k_opt)

    # 4. Calculate silhouette score in the UMAP space
    silh_scores  = cluster_silhouette(master_clustered, this_pathway )

    # 5. Or using PCA (user-specified matrix)
    umap_stats_pca <- makeUmapStats(master_clustered, silh_scores,
                                   this_pathway = this_pathway,
    															dist_method = 'user', user_dist_matrix = dist_pca)

    scatter_export = master_clustered


        # 6. Recursive pipeline
    recursive_res = recursiveSpectral(p_list, k_opt,n_motifs ,
    								silh_cutoff = min_silh , master_clustered,
    								spectrum_kernel = kernel_2nd, this_pathway=this_pathway )
    master_recursive = recursive_res$master # main data frame with meta data and labels
    scatter_export2= recursive_res$scatter_export # main data frame with umap coordinates


    # we use the recursive clustering silhouette scores
    # scatter_expor2 has the recursive labels for each profile
    umap_stats_recursive <- makeUmapStats(scatter_export2,
    													silh_scores = s2, dist_method = "user",
    													user_dist_matrix = dist_pca ,
    													this_pathway = this_pathway )


    return(list(recursive_res , master_recursive) )
}
