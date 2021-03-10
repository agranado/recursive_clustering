## functions Mar 2021
## Clustering and silhouette score statistics.
## Applies for any Seurat object

library(gridExtra)
library(stylo)
library(cluster)
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
