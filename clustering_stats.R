## functions Mar 2021
## Clustering and silhouette score statistics.
## Applies for any Seurat object

library(gridExtra)
library(stylo)
library(cluster)

#Apr 2021
# Quick pipeline
# Filters non-expressing cell types
# Filters also for devel or adult profiles

k_final = 25

genesPathway <- function(which_pathway = 'Notch'){

  this_pathway = all_pathways %>% dplyr::filter(pathway ==which_pathway ) %>% pull(gene)
  this_pathway = this_pathway[which(this_pathway %in% row.names(master_seurat))]

}

# 1. create devel/adult data.frame
# 2. Filter cell types with at least 2 expressed genes
# 3. Cluster cosine + ward.D2
# 4. returns clustered data.frame

quickPipeline <- function(which_pathway = 'Notch', master_seurat = c(),
                        k_final = 25,
                        min_genes_on = 1, min_expr = 0.2,
                        which_profiles = 'devel',
                        rand_control = F
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

        df_devel %>% dplyr::filter(genes_on>min_genes_on ) -> df_devel

        # Heatmap
        p = pheatmap(df_devel[,this_pathway ],
                     annotation_row = df_devel %>% select(dataset, Cell_class),
                     annotation_colors = colors_1206, show_rownames = T, fontsize = 5,
                     cutree_rows = k_final, clustering_method = 'ward.D2',
                      clustering_distance_rows = dist.cosine(as.matrix(df_devel[,this_pathway])),
                     cluster_cols = F, silent = T)


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

# data is already clustered by pathway
# only for the cell types in df_devel -- those with a pathway label
# Mar 25
global_clusterig<- function(df_devel = data.frame() , master_seurat = c() ,
														k_final = 25,n_pcs = 100 ){
        	# cells are filtered
        	df_devel$class_label <- df_devel$class_label %>% as.character()

        	# pathway colors
        	colors_1206$class_label = makeQualitativePal(length(df_devel$class_label %>% unique))
        	names(colors_1206$class_label) <- df_devel$class_label %>% unique

        	# PCA embedding
        	# Number of principal components
          # umap_coords by default uses cell_id as the row.name
        	umap_coords<- Embeddings(master_seurat,reduction = 'pca')
        	umap_coords <- umap_coords[,1:n_pcs]

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
         sapply(1:k_final, function(x){
            class_cells = df_devel$cell_id[df_devel$class_label==x]
            umap_coords[class_cells,] %>% dist %>% as.matrix -> dist_class

            dist_class[upper.tri(dist_class)] %>% mean()

        }) -> diversity_pathway

        #global stats // transcriptome
        # same idea, we use the global label -----
        global_labels = cutree(p_global$tree_col, k = k_final)
        df_devel$global_label = global_labels

        sapply(1:k_final, function(x){
            class_cells = df_devel$cell_id[df_devel$global_label==x]
            umap_coords[class_cells,] %>% dist %>% as.matrix -> dist_class

            dist_class[upper.tri(dist_class)] %>% mean()

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

    ss = c()
    for(k in 2:k_max){
        clustering = cutree(hclust(dist_mat, method =clust_method ), k)
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
        xlab("N clusters") + ylab(paste("Clustering Score (",score.type, ")",sep="")) +
        theme(legend.position= legend.pos) +
        theme(axis.text.x=element_text(size=15), axis.text.y = element_text(size =15)) +
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
