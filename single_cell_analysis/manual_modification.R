
integrated =   qread(file = paste0(configuration_scrnaseq@paths@output_archives,"/integrated_ccbr_annotation_hallmark.qs"))
no_integrated =   qread(file = paste0(configuration_scrnaseq@paths@output_archives,"/no_integrated_ccbr_annotation_hallmark.qs"))


DimPlot(integrated,reduction="fitsne",group.by=c("Original_Louvain","predicted.id"),label=T,split.by = "hypergroup")

DimPlot(no_integrated[[1]],reduction="fitsne",group.by=c("Original_Louvain","predicted.id"),label=T,split.by = "hypergroup")
DimPlot(no_integrated[[2]],reduction="fitsne",group.by=c("Original_Louvain","predicted.id"),label=T,split.by = "hypergroup")



library(scales)
DimPlot(no_integrated[[1]],reduction="fitsne",group.by=c("Original_Louvain","predicted.id"),label=T,split.by = "hypergroup",ncol = 2)

cluster_prediction = function(data) {
  
  mtdt = extract_metadata_reduction_seurat(data)
  return(ggplot(mtdt, aes(x=predicted.id,fill=Original_Louvain)) + geom_bar(aes(y=(..count..)/sum(..count..))) + coord_flip() +
           facet_wrap(~Original_Louvain,scales = "free") +
           scale_y_continuous(labels = percent))
}

barplot_cluster_prediction = lapply(c(integrated,no_integrated), function(x) cluster_prediction(x))

integrated_m =integrated
integrated_m@meta.data = integrated_m@meta.data %>%
  mutate(
    hypergroup = as.character(hypergroup),
    hypergroup = case_when(
      Original_Louvain %in% c(0,1,2,7,8,9,10,16,17,18,22) ~ "Immune Cells",
      Original_Louvain %in% c(21) ~ "Neural Cells",
      Original_Louvain %in% c(14) ~ "Tenocytes",
      Original_Louvain %in% c(13,23) ~ "Mural Cells",
      Original_Louvain %in% c(11,12,3) ~ "Endothelial Cells",
      Original_Louvain %in% c(5,6,19) ~ "FAPs",
      Original_Louvain %in% c(4,20,15,24) ~ "MuSCs",
    
    ))

hypergroup_prediction = function(data) {
  
  mtdt = extract_metadata_reduction_seurat(data)
  return(ggplot(mtdt, aes(x=Original_Louvain,fill=hypergroup)) + geom_bar(aes(y=(..count..)/sum(..count..))) + coord_flip() +
           facet_wrap(~hypergroup,scales = "free") +
           scale_y_continuous(labels = percent))
}
cluster_prediction(integrated_m)
hypergroup_prediction(integrated_m) + DimPlot(integrated,label=T,reduction="fitsne")

qsave(integrated_m,file = paste0(configuration_scrnaseq@paths@output_archives,"/integrated_manual_addingcorrections.qs"))

genes = c("Pecam1","Pax7","Myod1","Myog","Tnmd","Cd14","Col1a1","Gzma","Xcl1",
          "Rgs5","Irf7","Cd74","Itgam","Ptprc","Ctsa")

mtdt = do.call("rbind",lapply(no_integrated, extract_metadata_reduction_seurat))
genesExpression = do.call("rbind",lapply(no_integrated, function(x) FetchData(x,vars=genes)))

tmp = merge(mtdt,genesExpression,by=0) %>% pivot_longer(cols=genes,names_to = "gene",values_to="expression") %>% 
  data.frame() %>%
  mutate(gene = factor(gene,levels=genes))

ggplot() +
  geom_point(data = filter(tmp,sample=="KO"), aes(FItSNE_1,FItSNE_2),color="grey80") + 
  geom_point(data = filter(tmp,sample=="KO" & expression>0) %>% arrange(expression), aes(x=FItSNE_1,FItSNE_2,color=expression)) +
  scale_color_gradient2(low="grey60",high="darkblue",n.breaks=5) +
  facet_wrap(~gene,ncol = 4) +
  theme_legrand2

ggplot() +
  geom_point(data = filter(tmp,sample=="WT"), aes(FItSNE_1,FItSNE_2),color="grey80") + 
  geom_point(data = filter(tmp,sample=="WT" & expression>0) %>% arrange(expression), aes(x=FItSNE_1,FItSNE_2,color=expression)) +
  scale_color_gradient2(low="grey60",high="darkblue",n.breaks=5,na.value = "grey80") +
  facet_wrap(~gene,ncol = 4) +
  theme_legrand2
