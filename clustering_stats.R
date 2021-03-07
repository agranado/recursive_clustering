## functions Mar 2021
## Clustering and silhouette score statistics.
## Applies for any Seurat object

library(gridExtra)
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
