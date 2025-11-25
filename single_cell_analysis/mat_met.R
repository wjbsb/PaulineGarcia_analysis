raw = qs::qread("/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/single_cell/run23052023_starsolo/Rds_Rdata/listCSO_soupX.qs")

lapply(raw, function(x) min(x@meta.data$nCount_RNA))
lapply(raw, function(x) max(x@meta.data$nCount_RNA))
lapply(raw, function(x) dim(x))

filtered = qs::qread("/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/single_cell/run23052023_starsolo/Rds_Rdata/listCSO_afterthresholds_soupX.qs")

lapply(filtered@filtered_data@list_cso, function(x) min(x@meta.data$nCount_RNA))
lapply(filtered@filtered_data@list_cso, function(x) max(x@meta.data$nCount_RNA))
lapply(filtered@filtered_data@list_cso, function(x) dim(x))

doublet = qs::qread("/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/single_cell/run23052023_starsolo/Rds_Rdata/listcso2integrated_soupX.qs")

lapply(doublet, function(x) table(x@meta.data$scDblFinder.class))

doublet_f=list()
for (d in names(doublet)) {
  
  Idents(doublet[[d]]) = "scDblFinder.class"
  tmp = subset(doublet[[d]], subset = scDblFinder.class == "singlet")
  
  doublet_f[[d]] = tmp
}

lapply(doublet_f, function(x) min(x@meta.data$nCount_RNA))
lapply(doublet_f, function(x) max(x@meta.data$nCount_RNA))
lapply(doublet_f, function(x) median(x@meta.data$nCount_RNA))
lapply(doublet_f, function(x) min(x@meta.data$nFeature_RNA))
lapply(doublet_f, function(x) max(x@meta.data$nFeature_RNA))
lapply(doublet_f, function(x) median(x@meta.data$nFeature_RNA))

lapply(doublet_f, function(x) min(x@meta.data$percent.protein_coding))
lapply(doublet_f, function(x) max(x@meta.data$percent.protein_coding))
lapply(doublet_f, function(x) median(x@meta.data$percent.protein_coding))

lapply(doublet_f, function(x) min(x@meta.data$percent.rRNA))
lapply(doublet_f, function(x) max(x@meta.data$percent.rRNA))
lapply(doublet_f, function(x) median(x@meta.data$percent.rRNA))

lapply(doublet_f, function(x) min(x@meta.data$percent.MT))
lapply(doublet_f, function(x) max(x@meta.data$percent.MT))
lapply(doublet_f, function(x) median(x@meta.data$percent.MT))

lapply(doublet_f, function(x) dim(x))
FeaturePlot(doublet_f$KO,features=c("nCount_RNA","nFeature_RNA","percent.protein_coding","percent.rRNA","percent.MT"))

integrated = qs::qread("/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/single_cell/run23052023_starsolo/Rds_Rdata/integrated.qs")

min(integrated@meta.data$nCount_RNA)
max(integrated@meta.data$nCount_RNA)
dim(integrated)

FeaturePlot(integrated,features=c("nCount_RNA","nFeature_RNA","percent.protein_coding","percent.rRNA","percent.MT"),cols=wesanderson::wes_palette("Zissou1"))

config = readRDS("/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/single_cell/run23052023_starsolo/Rds_Rdata/config.Rds")
source("/home/mna_bioinfo/Bureau/development/github/single_cell_RNA_seq/functions_v2.R")


cso = filtered@filtered_data@list_cso$KO
cols = "percent.protein_coding"
path_custom_thresholds = config@paths@input_thresholds

#Custom thresholds
if (path_custom_thresholds!="NULL") {
  custom_threshold_by_sample <- read.table(path_custom_thresholds, sep="\t", header = T) %>%
    filter(str_detect(string=cso@project.name, pattern=sample))
  mtdt <- cso@meta.data %>% dplyr::select(-sample)
  mtdt$cells_in_custom_thresholds <- ifelse(
    #threshold
    (mtdt$nCount_RNA >=  custom_threshold_by_sample$count_down & mtdt$nCount_RNA <= custom_threshold_by_sample$count_up) &
      (mtdt$nFeature_RNA >= custom_threshold_by_sample$feature_down & mtdt$nFeature_RNA <= custom_threshold_by_sample$feature_up) &
      mtdt$percent.mt <= custom_threshold_by_sample$percent.mt_up,
    #output
    "Retained","Don't Retained"
  )
  p2 <- ggplot(mtdt, aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
    geom_point() +
    scale_color_gradientn(colours=c("green","blue","orange"),
                          na.value="red",
                          breaks=round(as.vector(sort(c(0,25,50,70))),2),
                          values=scales::rescale(c(MIN,MID,MAX)),limits=c(MIN,MAX)) +
    #ecriture scientifique
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    facet_wrap(~cells_in_custom_thresholds,labeller=as_labeller(labeller_modif(mtdt,"cells_in_custom_thresholds"))) +
    theme_legrand2 +
    labs(color="% MT")
  
  bcrank <- DropletUtils::barcodeRanks(counts(as.SingleCellExperiment(cso)), fit.bounds = c(custom_threshold_by_sample$count_down, custom_threshold_by_sample$count_up))
  knee <- bcrank@metadata$knee
  inflection <- bcrank@metadata$inflection
  
  
  yintercept_dot <- c("Knee"=knee,
                      "Inflexion"=inflection,
                      "Up custom thresholds"=custom_threshold_by_sample$count_down,
                      "Down custom thresholds"=custom_threshold_by_sample$count_up)
  yintercept_colors <- c("red", "blue", "black", "green")
  yintercept <- data.frame(yintercept_dot, yintercept_colors)
  
  
  tmp <- filter(mtdt,
                between(nCount_RNA,custom_threshold_by_sample$count_down,custom_threshold_by_sample$count_up) &
                  between(nFeature_RNA,custom_threshold_by_sample$feature_down,custom_threshold_by_sample$feature_up) &
                  percent.mt <= custom_threshold_by_sample$percent.mt_up
  )
  
  mtdt$complexityValue <- log10(mtdt$nFeature_RNA)/log10(mtdt$nCount_RNA)
  mad_metrics_custom_f <- setNames(lapply(colnames(mtdt), function(y) outlierHistogram(data=tmp, x=y)), colnames(mtdt))
} else {
  bcrank <- DropletUtils::barcodeRanks(counts(as.SingleCellExperiment(cso)))
  knee <- bcrank@metadata$knee
  inflection <- bcrank@metadata$inflection
  yintercept_dot <- c("Knee"=knee,
                      "Inflexion"=inflection)
  yintercept_colors <- c("red", "blue")
  yintercept <- data.frame(yintercept_dot, yintercept_colors)
  
  mtdt <- cso@meta.data %>% dplyr::select(-sample)
  
  mtdt$complexityValue <- log10(mtdt$nFeature_RNA)/log10(mtdt$nCount_RNA)
  # mad_metrics_raw <- setNames(lapply(colnames(mtdt), function(y) outlierHistogram(data=mtdt, x=y)), colnames(mtdt))
}
library(patchwork)
library(ggridges)
p=ggplot(cso@meta.data, aes(x=nCount_RNA,y=1)) +  
  ggridges::geom_density_ridges(jittered_points = T,
                                alpha = 0.4,
                                scale=1,
                                rel_min_height=.01,
                                point_shape="|",
                                point_size=3,size=.25,
                                position=position_points_jitter(height=0)) +
  geom_vline(data=make_position_around_best_max(cso$nCount_RNA,logtransform = T)[[2]],aes(xintercept=10^x,color=factor(val))) +
  geom_vline(xintercept = c(quantile(cso$nCount_RNA,0.1),
                            median(cso$nCount_RNA),
                            quantile(cso$nCount_RNA,0.9)),linetype="dashed",color="black") +
  # geom_vline(xintercept = c(1,
  #                           10),linetype="dashed",color="darkred")
  scale_x_log10() +
  ggplot(cso@meta.data, aes(x=nFeature_RNA,y=1)) + 
  ggridges::geom_density_ridges(jittered_points = T, 
                                alpha = 0.4,
                                scale=1,rel_min_height=.01,
                                point_shape="|",
                                point_size=3,size=.25,
                                position=position_points_jitter(height=0)) +
  geom_vline(data=make_position_around_best_max(cso$nFeature_RNA,logtransform = T)[[2]],aes(xintercept=10^x,color=factor(val))) +
  geom_vline(xintercept = c(quantile(cso$nFeature_RNA,0.1),
                            median(cso$nFeature_RNA),
                            quantile(cso$nFeature_RNA,0.9)),linetype="dashed",color="black") +
  scale_x_log10() +
  ggplot(cso@meta.data, aes(x=nCount_RNA,y=nFeature_RNA)) +
  geom_point(alpha=.3) +
  scale_x_log10() +
  scale_y_log10() +
  # geom_hex(bindwidth=c(1,1)) +
  # scale_fill_gradient(low="darkblue",high="darkorange")
  geom_density2d(adjust=.1) 

sample = cso@meta.data$sample %>% unique()
ggsave(filename = paste0("view_metrics_",sample,".pdf"),
       plot = p,device = cairo_pdf,
       path = paste0(configuration_scrnaseq@paths@output_plots,"/DataExploration"),width = 15,height = 10,dpi=1000)
#  
# library(mclust)
#  mixmodel = Mclust(log10(cso$nCount_RNA),model="V")
#  mu = mixmodel$par$mean
#  sigma = sqrt(mixmodel$par$variance$sigmasq) 
#  # Graphique
#  mclust1Dplot(log10(cso$nCount_RNA), parameters=mixmodel$parameters, z=mixmodel$z, what="density", col='red', new=TRUE)
#  hist(log10(cso$nCount_RNA), freq=FALSE, breaks=100, add=TRUE)
#  abline(v=mu, col='blue')
#  abline(v=c(mu+sigma,mu-sigma), col='blue', lty=2)

print("calculate limits")

ncount_limits <- make_position_around_best_max(cso$nCount_RNA,T)[[1]]
nfeature_limits <- make_position_around_best_max(x = cso$nFeature_RNA,logtransform = T)[[1]]

if (any(c(is.null(ncount_limits)), is.null(nfeature_limits))) {
  ncount_limits <- make_position_around_best_max(cso$nCount_RNA,F)[[1]]
  nfeature_limits <- make_position_around_best_max(x = cso$nFeature_RNA,logtransform = F)[[1]]
}

if (is.null(nfeature_limits)) {
  nfeature_limits = data.frame(
    yintercept_dot=c(min(cso$nFeature_RNA),max(cso$nFeature_RNA)),
    yintercept_clors=NA,
    row.names=c("min_elbow","max_elbow")
  )
}

# if (is.na(ncount_limits$yintercept_dot[3])) {ncount_limits$yintercept_dot[3]=ncount_limits$yintercept_dot[2]}
# if (is.na(nfeature_limits$yintercept_dot[3])) {nfeature_limits$yintercept_dot[3]=nfeature_limits$yintercept_dot[2]}


# p1 <- ggplot(bcrank %>% data.frame() %>% arrange(desc(rank)), aes(rank, total)) +
#   geom_point() + 
#   scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
#                 labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#   scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
#                 labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#   geom_hline(yintercept = c(metadata(bcrank)$knee,
#                             metadata(bcrank)$inflection), 
#              color = c("blue","green")) +
#   theme_legrand2 + 
#   labs(y="nCount_RNA",x="Barcodes rank")
df <- cso@meta.data %>%
  mutate(cells_pass = case_when(
    
    nCount_RNA <=  metadata(bcrank)$inflection ~ "Not Available - Inflexion",
    nCount_RNA <= ncount_limits$yintercept_dot[1] & nFeature_RNA <= nfeature_limits$yintercept_dot[1] ~ "Low",
    between(nCount_RNA,ncount_limits$yintercept_dot[1],ncount_limits$yintercept_dot[2]) & 
      between(nFeature_RNA,nfeature_limits$yintercept_dot[1],nfeature_limits$yintercept_dot[2]) ~ "Medium",
    nCount_RNA > ncount_limits$yintercept_dot[2] & nFeature_RNA > nfeature_limits$yintercept_dot[2] ~ "High",
    nCount_RNA > ncount_limits$yintercept_dot[2] & nFeature_RNA < nfeature_limits$yintercept_dot[2] ~ "Medium+",
    between(x = nCount_RNA,left = ncount_limits$yintercept_dot[1],right = ncount_limits$yintercept_dot[2]) &
      nFeature_RNA >= nfeature_limits$yintercept_dot[2] ~ "Medium++",  
    nCount_RNA > ncount_limits$yintercept_dot[1] & nFeature_RNA < nfeature_limits$yintercept_dot[1] ~ "Medium-",
    .default = "Not Available - Default"
    
  ))

table(df$cells_pass)

p3 <- ggplot(df, aes(nCount_RNA, nFeature_RNA, color=cells_pass)) +
  geom_point() +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  geom_vline(xintercept = c(metadata(bcrank)$knee,
                            metadata(bcrank)$inflection), 
             color = c("blue","green")) +
  geom_vline(xintercept = ncount_limits$yintercept_dot, color = "black",linetype="dashed") +
  geom_hline(yintercept = nfeature_limits$yintercept_dot, color = "black",linetype="dashed") +
  theme_legrand2 
  

p4 <- ggplot(merge(bcrank %>% data.frame(), df,by=0), aes(rank, total,color=cells_pass)) +
  geom_point() +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  geom_hline(yintercept = c(metadata(bcrank)$knee,
                            metadata(bcrank)$inflection,
                            ncount_limits$yintercept_dot), 
             color = c("blue","green","black","black")) +
  theme_legrand2 + 
  labs(y="nCount_RNA",x="Barcodes rank")


ggsave(filename = paste0("view_filtering_",sample,".pdf"),
       plot = p3 + p4,device = cairo_pdf,
       path = paste0(configuration_scrnaseq@paths@output_plots,"/DataExploration"),width = 15,height = 10,dpi=1000)



if (old_thresholds) {
  #old thresholds 
  mtdt_f <- filter(df, cells_pass == "Medium") %>% dplyr::select(-sample, -cells_pass)
} else {
  #new thresholds
  mtdt_f <- filter(df, cells_pass %in% c("Low","Medium-","Medium","Medium+","Medium++","High")) %>% dplyr::select(-sample, -cells_pass)
}

if (!features_filtering_mincount) {
  # mahalanobis_raw <- mahalanobis_plot(dataset = mtdt_f,cols = c("percent.protein_coding","nCount_RNA","nFeature_RNA"))
  mahalanobis_raw <- mahalanobis_plot(dataset = mtdt_f,cols = c("percent.protein_coding","nCount_RNA","nFeature_RNA"))
} else {
  mahalanobis_raw <- mahalanobis_plot(dataset = mtdt_f,cols = c(paste0("percent.",cols),"nCount_RNA","nFeature_RNA"))
  
}

mahalanobis_raw$data2plot <- mahalanobis_raw$data2plot %>%
  mutate(state = case_when(
    pval >= 0.001 & outliers_mahalanobis == "Retained\n" ~ "Retained",
    pval >= 0.001 & outliers_mahalanobis == "Don't retained\n" ~ "Don't retained 1",
    pval < 0.001 & outliers_mahalanobis == "Retained\n" ~ "Don't retained 2",
    pval < 0.001 & outliers_mahalanobis == "Don't retained\n" ~ "Don't retained 3"
  ))

if (configuration_scrnaseq@parameters@species=="Dme") {
  method_single_transcriptomique = "sn"
} else {
  method_single_transcriptomique = "sc"
}


if (method_single_transcriptomique == "sc") {
  print("single cell")
  MIN <- min(na.omit(mahalanobis_raw$data2plot$mahalanobis)) %>% log10() %>% round(digits=2)
  MID <- quantile(na.omit(mahalanobis_raw$data2plot$mahalanobis), .5) %>% log10() %>%  round(digits=2)
  MAX <- max(na.omit(mahalanobis_raw$data2plot$mahalanobis)) %>% log10() %>% round(digits=2)
  
  g1 <- ggplot(mahalanobis_raw$data2plot, aes(log10(mahalanobis), fill=state)) +
    geom_density(alpha=0.5) + 
    geom_vline(xintercept = 10,linetype="dashed",color="red") +
    theme_bw() +
    labs(x=TeX("$Log_{10}(Mahalanobis\\ Distance)$"),y="Density")
  g2 <- ggplot(mahalanobis_raw$data2plot, aes(x=Dim.1, y=Dim.2, color=log10(mahalanobis))) +
    geom_point(size=0.7,alpha=0.5) +
    scale_color_gradientn(colours=c("green","blue","orange"),
                          na.value="red",
                          breaks=round(as.vector(sort(c(MIN,MID,MAX))),2),
                          values=scales::rescale(c(MIN,MID,MAX)),limits=c(MIN,MAX)) +
    geom_hline(yintercept = 0,linetype="dashed",color="black") +
    geom_vline(xintercept = 0,linetype="dashed",color="black") +
    theme_bw() +
    facet_wrap(~state, labeller=as_labeller(labeller_modif(mahalanobis_raw$data2plot,"state"))) +
    labs(
      x=paste0("PCA_1 (",
               round(mahalanobis_raw$data2pca$eig[1,2],2),"%)"),
      y=paste0("PCA_2 (",
               round(mahalanobis_raw$data2pca$eig[2,2],2),"%)"),
      color=TeX("$Log_{10}(Mahalanobis\\ Distance)$"))
} else {
  print("single nuclei")
  MIN <- min(na.omit(mahalanobis_raw$data2plot$mahalanobis)) %>% round(digits=2)
  MID <- quantile(na.omit(mahalanobis_raw$data2plot$mahalanobis), .5) %>% round(digits=2)
  MAX <- max(na.omit(mahalanobis_raw$data2plot$mahalanobis)) %>% round(digits=2)
  
  g1 <- ggplot(mahalanobis_raw$data2plot, aes(mahalanobis, fill=state)) +
    geom_density(alpha=0.5) + 
    geom_vline(xintercept = 10,linetype="dashed",color="red") +
    theme_bw() +
    labs(x=TeX("$Mahalanobis\\ Distance$"),y="Density")
  g2 <- ggplot(mahalanobis_raw$data2plot, aes(x=Dim.1, y=Dim.2, color=mahalanobis)) +
    geom_point(size=0.7,alpha=0.5) +
    scale_color_gradientn(colours=c("green","blue","orange"),
                          na.value="red",
                          breaks=round(as.vector(sort(c(MIN,MID,MAX))),2),
                          values=scales::rescale(c(MIN,MID,MAX)),limits=c(MIN,MAX)) +
    geom_hline(yintercept = 0,linetype="dashed",color="black") +
    geom_vline(xintercept = 0,linetype="dashed",color="black") +
    theme_bw() +
    facet_wrap(~state, labeller=as_labeller(labeller_modif(mahalanobis_raw$data2plot,"state"))) +
    labs(
      x=paste0("PCA_1 (",
               round(mahalanobis_raw$data2pca$eig[1,2],2),"%)"),
      y=paste0("PCA_2 (",
               round(mahalanobis_raw$data2pca$eig[2,2],2),"%)"),
      color=TeX("$Mahalanobis\\ Distance$"))
}


require(patchwork)
g_one <- g1 / g2

barcode <- list(
  raw_retained = rownames(mahalanobis_raw$data2plot[mahalanobis_raw$data2plot$state=="Retained",])
)

save_g_one=T
if (save_g_one) {
  qsave(g_one,file = paste0(configuration_scrnaseq@paths@output_plots,"/DataExploration/",cso$sample %>% unique(),"_mahalnobis_droplet.qs"))  
}

ggsave(path = paste0(configuration_scrnaseq@paths@output_plots,"/DataExploration"),
       filename=paste0(cso$sample %>% unique(),"_mahalnobis_droplet.pdf"),
       device = cairo_pdf,
       dpi=1000,width = 20,height = 10,
       plot=g_one)
cso_final <- subset(cso, cells = barcode$raw_retained)

cso_final = AddMetaData(cso_final,mahalanobis_raw$data2plot[,c("mahalanobis","pval")])

return(cso_final) 
}
calculate_ranks_add_thresholds <- function(cso, path_custom_thresholds) {
  
  cso@meta.data <- cso@meta.data %>%
    mutate(percent.mt = as.numeric(unlist(PercentageFeatureSet(cso, pattern = "^MT-|^mt"))),
           countPerfeature = nCount_RNA/nFeature_RNA,
           rank_nCount_RNA = rank(nCount_RNA),
           rank_nFeature_RNA = rank(nFeature_RNA),
           rank_percent.mt = rank(-percent.mt),
           rank_countPerfeature = rank(countPerfeature))
  cso@meta.data$mean_ranks <- rowMeans(cso@meta.data[, c("rank_nCount_RNA", "rank_nFeature_RNA", "rank_percent.mt", "rank_countPerfeature")])
  
  summary_ranks <- summary(cso@meta.data$mean_ranks)
  
  cso@meta.data$cells_in_threshold25_75 <- ifelse(cso@meta.data$mean_ranks > summary_ranks[2] & cso@meta.data$mean_ranks < summary_ranks[5], T, F)
  
  custom_threshold_by_sample <- read.table(path_custom_thresholds, sep="\t", header = T) %>%
    filter(str_detect(string=cso@project.name, pattern=sample))
  
  mad_ncount_rna <- mad(cso@meta.data$nCount_RNA)
  mad_nfeature_rna <- mad(cso@meta.data$nFeature_RNA)
  
  cso@meta.data$cells_in_custom_thresholds <- ifelse(
    #threshold
    (cso@meta.data$nCount_RNA >=  custom_threshold_by_sample$count_down & cso@meta.data$nCount_RNA <= custom_threshold_by_sample$count_up) &
      (cso@meta.data$nFeature_RNA >= custom_threshold_by_sample$feature_down & cso@meta.data$nFeature_RNA <= custom_threshold_by_sample$feature_up) &
      cso@meta.data$percent.mt <= custom_threshold_by_sample$percent.mt_up, #& cso@meta.data$countPerfeature >=3,
    #output
    T,F
  )
  cso@meta.data$count <- ifelse(cso@meta.data$nCount_RNA >= custom_threshold_by_sample$count_down & cso@meta.data$nCount_RNA <= custom_threshold_by_sample$count_up, T, F)
  cso@meta.data$feature <- ifelse(cso@meta.data$nFeature_RNA >=  custom_threshold_by_sample$feature_down & cso@meta.data$nFeature_RNA <= custom_threshold_by_sample$feature_up, T, F)
  cso@meta.data$mt <- ifelse(cso@meta.data$percent.mt <= custom_threshold_by_sample$percent.mt_up, T, F)
  cso@meta.data$complexity <- ifelse(cso@meta.data$countPerfeature >= 3, T, F)
  return(cso)
}