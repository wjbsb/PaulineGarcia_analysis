plot_ok4presentation = function(plot,output=NA,ncol,nrow) {
  # made 1 page with 2 rows&columns
  # plot are in [1,1]. Legend are in [1,2]
  
  tmp = plot + 
    theme(legend.position="none")
  legend <- get_legend(plot + 
                         guides(color = guide_legend(override.aes = list(size = 5),ncol=ncol,nrow = nrow)) + 
                         labs(color="")) 
  # tmp$layers[[2]] = NULL
  l = list(tmp,legend)
  p=marrangeGrob(grobs = l, nrow=1, ncol=2,top="")
  if (!is.na(output)) {
    ggsave(output, 
           p,width = 15, height = 7, dpi=1000,units="in",
           device = cairo_pdf)
  }
  
  return(p)
}

integrated = qread(file=paste0(configuration_scrnaseq@paths@output_archives,"/integrated_ccbr_annotation_hallmark.qs"),nthreads=20)

integrated@meta.data = integrated@meta.data %>% 
  mutate(new_annotation = case_when(
    Original_Louvain == 0 ~ "Trem2+ Macro",
    Original_Louvain == 1 ~ "H2-Ab1 Macro",
    Original_Louvain == 2 ~ "Monocytes",
    Original_Louvain == 3 ~ "Capillary ECs",
    Original_Louvain == 4 ~ "MuSCs",
    Original_Louvain == 5 ~ "FAPs_2",
    Original_Louvain == 6 ~ "FAPs_1",
    Original_Louvain == 7 ~ "CD14+ Macro",
    Original_Louvain == 8 ~ "SPARC+ Macro",
    Original_Louvain == 9 ~ "Lars2+ macro",
    Original_Louvain == 10 ~ "H2-Ab1_2",
    Original_Louvain == 11 ~ "EC cap 2",
    Original_Louvain == 12 ~ "EC vessel",
    Original_Louvain == 13 ~ "Mural",
    Original_Louvain == 14 ~ "Tenocytes",
    Original_Louvain == 15 ~ "Cycling MuSCs",
    Original_Louvain == 16 ~ "Ly6c2+ macro",
    Original_Louvain == 17 ~ "T cells",
    Original_Louvain == 18 ~ "Dendritic cells",
    Original_Louvain == 19 ~ "Cycling FAPs",
    Original_Louvain == 20 ~ "Dying cells",
    Original_Louvain == 21 ~ "Schwann cells",
    Original_Louvain == 22 ~ "NK",
    Original_Louvain == 23 ~ "Mural 2",
    Original_Louvain == 24 ~ "Myonuc")) %>% 
  mutate(new_subtype2 = paste0(Original_Louvain," : ",new_annotation))

col_order = dplyr::select(integrated@meta.data,Original_Louvain,new_subtype2) %>% 
  unique() %>% arrange(Original_Louvain) %>% pull(new_subtype2)

integrated@meta.data = integrated@meta.data %>% mutate(new_subtype2 = factor(new_subtype2,levels = col_order))


integrated@meta.data = integrated@meta.data %>% mutate(
  sample = factor(sample, 
                  levels=c("KO","WT"),
                  labels=c("Setdb1^{MuSC-mutate(new_subtype = paste0(Original_Louvain," : ",subtype))KO} ","Control "))
)

predictions <- integrated@meta.data[,c("sample","condition","lineage","Original_Louvain",
                                       colnames(integrated@meta.data)[grep("prediction",colnames(integrated@meta.data))])] %>%
  rownames_to_column(var="barcode") %>%
  pivot_longer(cols=colnames(integrated@meta.data)[grep("prediction",colnames(integrated@meta.data))],
               names_to="prediction",values_to="score_prediction")

integrated$predicted.id[integrated$prediction.score.max<0.5] = "low_prediction"

p1 = make_plots(data = integrated,
                expressed_genes_list = c("Setdb1"),
                idents_data = c("predicted.id","Original_Louvain","hypergroup"),
                groupby = "predicted.id")

tmp = extract_metadata_reduction_seurat(integrated) %>%
  mutate(sample = ifelse(sample == "KO", "Setdb1^{MuSC-KO}","Control") %>% factor())

alphaP=1
sizeP=2
sizePtext=0.36*20
sizetext=15
family="Helvetica"

palettes <- randomcoloR::distinctColorPalette(length(unique(tmp$new_subtype2)))
p_cluster_integrated = 
  ggplot(tmp, aes(x=FItSNE_1,y=FItSNE_2,color=new_subtype2)) +
  geom_point(alpha=alphaP,size=sizeP,show.legend = T) +
  theme(legend.position = "right") +
  guides(color=guide_legend(ncol=1,title.position = "top"),
         fill=guide_legend(override.aes = aes(label=""))) +
  scale_color_manual(values=palettes) +
  scale_fill_manual(values=palettes) +
  labs(color="Clusters", x="FItSNE_1", y="FItSNE_2") +
  theme_classic() +
  theme(
    text=element_text(family = family),
    strip.text.x = element_text(face="bold.italic",size=sizetext,family = family),
    strip.background = element_blank(),
    legend.title = element_text(face="bold"),
    axis.title.x = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  )

plot_ok4presentation(p_cluster_integrated,
                     output = "/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier/integrated_louvain_v2.pdf",
                     ncol = 5,nrow = 5)

p_view_sample_unified = ggplot(tmp) +
  geom_point(aes(x=FItSNE_1,y=FItSNE_2,color=sample),size=2) +
  scale_color_manual(values=c("blue","red"),
                     labels=expression(Control,
                                       Setdb1^{MuSC-KO})) +
  theme_legrand2 +
  labs(color="Sample")


p_view_sample=ggplot() +
  geom_point(data=dplyr::select(tmp,-sample), aes(FItSNE_1,FItSNE_2),color="grey80") +
  geom_point(data=tmp,aes(x=FItSNE_1,y=FItSNE_2,color=sample),size=2) +
  scale_color_manual(values=c("blue","red"),
                     labels=expression(Control,
                                       Setdb1^{MuSC-KO})) +
  labs(color="Sample") +
  facet_wrap(~sample,labeller = as_labeller(labeller_modif(tmp,"sample"),label_parsed)) +
    theme_legrand2

ggsave(plot = p_view_sample_unified, path = "/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier/",
       filename="fig_samples_unified.pdf",width = 10, height = 7, dpi=1000,units="in",
       device = cairo_pdf)
ggsave(plot = p_view_sample, path = "/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier/",
       filename="fig_samples_separated.pdf",width = 10, height = 7, dpi=1000,units="in",
       device = cairo_pdf)
no_integrated = qread(file=paste0(configuration_scrnaseq@paths@output_archives,"/no_integrated_ccbr_annotation.qs"),nthreads=20)
names(no_integrated) = c("KO","WT")

no_integrated$WT@meta.data = no_integrated$WT@meta.data %>%
  mutate(
    subtype = case_when(
    
    Original_Louvain == "0" ~ "Endothelial Cells (capillary)",
    Original_Louvain == "1" ~ "Monocytes",
    Original_Louvain == "2" ~ "MØ (H2-Ab1+)",
    Original_Louvain == "3" ~ "MØ (Trem2+)",
    Original_Louvain == "4" ~ "FAPs & Tenocytes",
    Original_Louvain == "5" ~ "MØ (CD14+)",
    Original_Louvain == "6" ~ "MuSCs",
    Original_Louvain == "7" ~ "MØ (Lars2+)",
    Original_Louvain == "8" ~ "Dendritic Cells",
    Original_Louvain == "9" ~ "MØ (Ly6c2+)",
    Original_Louvain == "10" ~ "Mural Cells",
    Original_Louvain == "11" ~ "Endothelial Cells (vessel)",
    Original_Louvain == "12" ~ "Cycling MuSCs",
    Original_Louvain == "13" ~ "T Lymphocytes",
    Original_Louvain == "14" ~ "Myonuclei",
    Original_Louvain == "15" ~ "Natural Killer Cells",
    Original_Louvain == "16" ~ "Cycling FAPs",
    Original_Louvain == "17" ~ "Schwann Cells"
  ),
  new_hypergroup = case_when(
    Original_Louvain %in% c("0","11") ~ "Endothelial Cells",
    Original_Louvain %in% c("4","16") ~ "FAPs & Tenocytes",
    Original_Louvain == "10" ~ "Mural Cells",
    Original_Louvain %in% c("1","2","3","5","7","8","9") ~ "Macrophages",
    Original_Louvain %in% c("6","12","14") ~ "MuSCs & Myonuclei",
    Original_Louvain %in% c("13","15") ~ "Lymphocytes",
    Original_Louvain == "17" ~ "Schwann Cells"
  ))

DimPlot(no_integrated$KO, group.by = "Original_Louvain",reduction="fitsne_by_pca",label=T)

no_integrated$KO@meta.data = no_integrated$KO@meta.data %>%
  mutate(
    subtype = case_when(
      
      Original_Louvain == "0" ~ "MØ (Sparc+)",
      Original_Louvain == "1" ~ "MØ (H2-Ab1+)",
      Original_Louvain == "2" ~ "MØ (Trem2+)",
      Original_Louvain == "3" ~ "Monocytes",
      Original_Louvain == "4" ~ "IntMono/MØ",
      Original_Louvain == "5" ~ "Endothelial Cells (capillary)",
      Original_Louvain == "6" ~ "MØ (CD14+)",
      Original_Louvain == "7" ~ "FAPs",
      Original_Louvain == "8" ~ "MuSCs",
      Original_Louvain == "9" ~ "MØ (Lars2+)",
      Original_Louvain == "10" ~ "Myonuclei",
      Original_Louvain == "11" ~ "Tenocytes",
      Original_Louvain == "12" ~ "Dendritic Cells",
      Original_Louvain == "13" ~ "Endothelial Cells (vessel)",
      Original_Louvain == "14" ~ "T Lymphocytes",
      Original_Louvain == "15" ~ "Mural Cells",
      Original_Louvain == "16" ~ "MØ (Sparc+)",
      Original_Louvain == "17" ~ "Schwann Cells"
    ),
    new_hypergroup = case_when(
      Original_Louvain %in% c("5","13") ~ "Endothelial Cells",
      Original_Louvain %in% c("7","11") ~ "FAPs & Tenocytes",
      Original_Louvain == "15" ~ "Mural Cells",
      Original_Louvain %in% c("0","1","2","3","4","6","9","12") ~ "Macrophages",
      Original_Louvain %in% c("8","10") ~ "MuSCs & Myonuclei",
      Original_Louvain == c("14") ~ "Lymphocytes",
      Original_Louvain == "17" ~ "Schwann Cells"
    ),
    new_Original_Louvain = case_when(
      Original_Louvain == "0" ~ "0",
      Original_Louvain == "1" ~ "1",
      Original_Louvain == "2" ~ "2",
      Original_Louvain == "3" ~ "3",
      Original_Louvain == "4" ~ "4",
      Original_Louvain == "5" ~ "5",
      Original_Louvain == "6" ~ "6",
      Original_Louvain == "7" ~ "7",
      Original_Louvain == "8" ~ "8",
      Original_Louvain == "9" ~ "9",
      Original_Louvain == "10" ~ "10",
      Original_Louvain == "11" ~ "11",
      Original_Louvain == "12" ~ "12",
      Original_Louvain == "13" ~ "13",
      Original_Louvain == "14" ~ "14",
      Original_Louvain == "15" ~ "15",
      Original_Louvain == "16" ~ "0",
      Original_Louvain == "17" ~ "16"))
no_integrated$KO@meta.data = no_integrated$KO@meta.data %>% mutate(new_subtype = paste0(Original_Louvain," : ",subtype))


qsave(integrated,file=paste0(configuration_scrnaseq@paths@output_archives,"/integrated_manual_addingcorrectionns.qs"),nthreads = 20)
qsave(no_integrated,file=paste0(configuration_scrnaseq@paths@output_archives,"/no_integrated_manual_addingcorrectionns.qs"),nthreads = 20)

no_integrated=qread(file=paste0(configuration_scrnaseq@paths@output_archives,"/no_integrated_manual_addingcorrectionns.qs"),nthreads = 20)

integrated=qread(file=paste0(configuration_scrnaseq@paths@output_archives,"/integrated_manual_addingcorrectionns.qs"),nthreads = 20)




#######WT
#EC
#FAPs
#musc
#mural
#lymphocyte
#macrom
#schwann
genes_subtype = list(WT=c("Pecam1","Aqp1","Selp","Lrg1",
                          "Pdgfra","Dcn","Ccna2","Cenpf",
                          "Myog","Pax7",
                          "Rgs5","Myh11",
                          "Gzma","Nkg7","Ccl5","Klrd1",
                          "Folr2","Ccl8","H2-Ab1","Cd74","Trem2","Scamp1","Cd14","Gdf15","Lars2","Zeb2","Cd209a","Plc8","Ly6c2","S100a9",
                          "Plp1","Kcna1"),
                     KO=c("Pecam1","Aqp1","Selp","Lrg1",
                          "Pdgfra","Dcn","Tnmd","Thbs4","Mki67","Apod",
                          "Myog","Gsta4","Acta1","Tnnc2",
                          "Rgs5","Myh11",
                          "Gzma","Nkg7",
                          "Postn","Fbxo38","H2-Ab1","Cd74","Trem2","Spp1","Folr2","Ccl8","Casp8","Tor1b","Rgs5","Cd14","Gdf15","Cd209a","Plac8",
                          "Plp1","Kcna1"))
#######KO
#EC
#FAPs
#musc
#mural
#lymphocyte
#macro
#schwann
Idents(no_integrated$WT) = "subtype"
Idents(no_integrated$KO) = "subtype"
make_bubbleplot = function(data,genelist) {
  # data=no_integrated$KO
  # genelist=genes_subtype$KO
  mtdt = extract_metadata_reduction_seurat(data) %>% dplyr::select(subtype,new_hypergroup) %>% unique()
  DEA = FindAllMarkers(data,assay="RNA",features = genelist,test.use = "wilcox")

  genelist = genelist[unique(DEA$gene) %in% genelist]
  
  MIN=min(na.omit(DEA$avg_log2FC))
  MID=quantile(na.omit(DEA$avg_log2FC),.5)
  MAX=max(na.omit(DEA$avg_log2FC))
  
  DEA$cluster = as.character(DEA$cluster)
  
  DEA = inner_join(DEA,mtdt,by=c("cluster"="subtype")) %>% arrange(new_hypergroup)
  
  lvl = dplyr::select(DEA, cluster, new_hypergroup) %>% unique() %>% arrange(new_hypergroup)
  subtype_order = lvl %>% pull(cluster)
  hypergroup_order = lvl %>% pull(new_hypergroup) %>% unique()
  
  library(ggh4x)
  stripx <- strip_themed(background_x = elem_list_rect(fill = rainbow(2)))
  
  DEA = DEA %>% mutate(cluster = factor(cluster,levels=subtype_order),
                       new_hypergroup = factor(new_hypergroup, levels=hypergroup_order),
                       typegene = "Subtype Markers",
                       gene = factor(gene,levels=genelist))
  
  
  p = make_plots(data,expressed_genes_list = "Setdb1",idents_data = c("subtype","new_hypergroup"))
  p$palettes$new_hypergroup = p$palettes$new_hypergroup[levels(data$new_hypergroup)]
  
  l = lapply(p$palettes$new_hypergroup, function(x) element_rect(fill=x))
  
  
  return(ggplot(DEA, aes(x=gene,y=cluster)) + 
    # geom_point(data = filter(DEA,p_val_adj<0.05), aes(color=avg_log2FC,size=pct.1*100)) +
    geom_point(data = filter(DEA,p_val_adj<0.05), aes(color=avg_log2FC,size=pct.1*100)) +
    geom_point(data = filter(DEA,p_val_adj>=0.05), color="grey80",size=5) +
    scale_size_continuous(labels = scales::percent_format(scale = 1)) +
    scale_color_gradientn(colours=c("blue","white","red"),
                          na.value="grey80",
                          values=scales::rescale(c(MIN,MID,MAX)),limits=c(MIN,MAX)) +
    labs(x="",y="",color="Average Log2 Fold Change",size="% expression in cluster") +
    theme_classic() +
    facet_grid(~typegene,scales="free",space="free",shrink=T,switch="both") +
    facet_grid2(new_hypergroup~typegene,scales = "free",space = "free",shrink = T,switch = "both",
                strip = strip_themed(
                  background_y  = l)) +
    theme(axis.text.x = element_text(angle=45,hjust = 1,family = "Helvetica",size=11,face = "bold.italic"),
          axis.text.y = element_text(family="Helvetica",size=11,face="bold.italic"),
          strip.text.x.bottom = element_text(family = "Helvetica",size=11,face="bold.italic"),
          strip.background.x = element_rect(),
          strip.text.y.left = element_text(angle=0,family="Helvetica",size=11,face="bold.italic"),
          legend.position = "bottom")
  )
}

pko = make_bubbleplot(data=no_integrated$KO,genelist=genes_subtype$KO)
pwt = make_bubbleplot(data=no_integrated$WT,genelist=genes_subtype$WT)


no_integrated$KO@meta.data = no_integrated$KO@meta.data %>%
  mutate(new_subtype = paste0(new_Original_Louvain," : ",subtype),
         new_subtype = factor(new_subtype, levels = c(
           "0 : MØ (Sparc+)","1 : MØ (H2-Ab1+)","2 : MØ (Trem2+)","3 : Monocytes","4 : IntMono/MØ",
           "5 : Endothelial Cells (capillary)", "6 : MØ (CD14+)","7 : FAPs","8 : MuSCs","9 : MØ (Lars2+)",
           "10 : Myonuclei","11 : Tenocytes","12 : Dendritic Cells","13 : Endothelial Cells (vessel)",
           "14 : T Lymphocytes","15 : Mural Cells","16 : Schwann Cells"
         )))
no_integrated$WT@meta.data = no_integrated$WT@meta.data %>%
  mutate(new_subtype = paste0(Original_Louvain," : ",subtype),
         new_Original_Louvain = Original_Louvain,
         new_subtype = factor(new_subtype, levels = c(
           "0 : Endothelial Cells (capillary)","1 : Monocytes","2 : MØ (H2-Ab1+)","3 : MØ (Trem2+)","4 : FAPs & Tenocytes",
           "5 : MØ (CD14+)","6 : MuSCs","7 : MØ (Lars2+)","8 : Dendritic Cells","9 : MØ (Ly6c2+)","10 : Mural Cells",
           "11 : Endothelial Cells (vessel)","12 : Cycling MuSCs","13 : T Lymphocytes","14 : Myonuclei",
           "15 : Natural Killer Cells","16 : Cycling FAPs","17 : Schwann Cells"
         )))
make_fig5 = function(data,pathoutput) {
  family="helvetica"
  p1 = make_plots(data = data,  
                  expressed_genes_list = c("Setdb1"),
                  idents_data = c("predicted.id","hypergroup","subtype","new_subtype"),
                  groupby = "predicted.id")
  
  mtdt = extract_metadata_reduction_seurat(data)
  
  fig = ggplot(mtdt, aes(x=FItSNE_1,y=FItSNE_2,color=new_subtype)) +
    geom_point(alpha=alphaP,size=sizeP,show.legend = T) +
    # ggrepel::geom_text_repel(data = mtdt_label_position(mtdt,"new_Original_Louvain","new_Original_Louvain","fitsne"),
    #                          aes(x=medX,y=medY,
    #                              label=new_Original_Louvain),
    #                          color="black",
    #                          fontface = "bold",
    #                          max.iter = 3e3,
    #                          seed = 42,
    #                          force=1) +
    theme(legend.position = "right") +
    guides(color=guide_legend(ncol=1,title.position = "top"),
           fill=guide_legend(override.aes = aes(label=""))) +
    scale_color_manual(values=p1$palettes$new_subtype) +
    scale_fill_manual(values=p1$palettes$new_subtype) +
    labs(color="") +
    theme_classic() +
    theme(
      #theme elements
      text=element_text(family = family),
      #plot elements
      strip.text.x = element_text(face="bold.italic",size=sizetext,family = family),
      strip.background = element_blank(),
      #facetting elements
      #legend elements
      legend.title = element_text(face="bold"),
      #axis elements
      axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold")
      #panel elements
    )
  require(gginnards)
  
  output= plot_ok4presentation(
    fig,
    pathoutput,
    ncol=3,
    nrow=6
  )
  return(list(output=output,palette = p1$palettes))
}

fig5d = make_fig5(no_integrated$KO,"/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier/fig5d.pdf")
fig5e = make_fig5(no_integrated$WT,"/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier/fig5e.pdf")


lgenes=c("Pecam1","Pax7","Myod1","Myog","Tnmd","Cd14","Col1a1",
        "Gzma","Xcl1","Nlrp3","Rgs5","Irf7","Cd74","Itgam","Ptprc","Ctsa")

make_figS5 = function(data,lgenes,log=T) {
  
  fgenes = merge(FetchData(data,vars = lgenes),extract_metadata_reduction_seurat(data),by=0) %>%
    pivot_longer(cols = lgenes,names_to="genes",values_to="expression") %>%
    mutate(genes=factor(genes, levels=lgenes))
  if (log) {
    p=ggplot() + 
      geom_point(data = fgenes%>% filter(expression<=1), 
                  # aes(FItSNE_1,FItSNE_2),col="grey70",size=0.5) +
                 aes(FItSNE_1,FItSNE_2),col="grey70",size=sizeP) +
      geom_point(data = fgenes %>% arrange(expression) %>% filter(Expression>1), 
                 aes(FItSNE_1,FItSNE_2,col=log10(expression)),size=sizeP) + 
      facet_wrap(~genes) +
      theme_legrand2 + 
      labs(color="Expression") +
      colorspace::scale_color_continuous_sequential(palette = "Blues")
  } else {
    p=ggplot() + 
      geom_point(data = fgenes%>% filter(expression<=1), 
                 aes(FItSNE_1,FItSNE_2),col="grey70",size=0.5) +
      geom_point(data = fgenes %>% arrange(expression) %>% filter(expression>1), 
                 aes(FItSNE_1,FItSNE_2,col=expression),size=sizeP) + 
      facet_wrap(~genes) +
      theme_legrand2 +
      labs(color="Expression") +
      colorspace::scale_color_continuous_sequential(palette = "Blues")
  }
  return(p)
}

figS5a = make_figS5(no_integrated$KO,lgenes,log = F)
figS5b = make_figS5(no_integrated$WT,lgenes,log=F)
  
ggsave(plot = figS5a, path = "/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier/",
       filename="figS5a",width = 20, height = 10, dpi=1000,units="in",
       device = cairo_pdf)
ggsave(plot=figS5b, path = "/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier/",
       filename="figS5b",width = 20, height = 10, dpi=1000,units="in",
       device = cairo_pdf)



trajectories = qs::qread(
      # file = paste0(configuration_scrnaseq@paths@output_archives,"/",hypergroup,"/trajectories_analysis.qs")
      file = "/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier/trajectories_macrophages.qs",
)

ggsave(
  filename="all_trajectories_Macrophages_selected.pdf",device = cairo_pdf,
  # path = configuration_scrnaseq@paths@output_plots,
  path = "/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier",
  width = 10,height = 5, 
  plot = trajectories$topologies_analysis$plot$pbytrajectories
)

immune_cells=qread(file = paste0(configuration_scrnaseq@paths@output_archives,"/Immune Cells/integrated.qs"))
remove=read.csv(file="/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/single_cell/run23052023_starsolo/Rds_Rdata/Immune Cells/metadata2remove.csv",header=T,sep=",")


immune_cells@meta.data = immune_cells@meta.data %>%
  mutate(type=case_when(
    
    Original_Louvain %in% remove$Original_Louvain ~ "Lymphocytes",
    .default = "Macrophages"
    
  ))
Idents(immune_cells)="type"

macrophages = subset(immune_cells,idents = "Macrophages")
DimPlot(macrophages,reduction="fitsne_by_harmony",group.by="Original_Louvain")




genes_macro = c("Folr2","Ccl8","H2-Ab1","Cd74","Trem2","Scamp1",
                "Cd14","Gdf15","Lars2","Zeb2","Cd209a","Plac8",
                "Ly6c2","S100a9","Postn","Fbxo38","Casp8","Tor1b",
                "Rgs5","Sparc")

macrophages@meta.data = macrophages@meta.data %>%
  mutate(
    subtype = case_when(
      Original_Louvain == "0" ~ "MØ (H2-Ab1+)",
      Original_Louvain == "1" ~ "Monocytes",
      Original_Louvain == "2" ~ "MØ (Trem2+ A)",
      Original_Louvain == "3" ~ "MØ (Trem3+ B)",
      Original_Louvain == "4" ~ "MØ (CD14+)",
      Original_Louvain == "5" ~ "MØ (Lars+2 A)",
      Original_Louvain == "6" ~ "MØ (Sparc+ A)",
      Original_Louvain == "7" ~ "T_Cells A",
      Original_Louvain == "8" ~ "intMono/MØ",
      Original_Louvain == "9" ~ "MØ (Ly6c2+)",
      Original_Louvain == "10" ~ "MØ (Sparc+ B)",
      Original_Louvain == "11" ~ "Dendritic ?",
      Original_Louvain == "12" ~ "NK_Cells B",
      Original_Louvain == "13" ~ "MØ (Phagocytic)",
      Original_Louvain == "14" ~ "Pericytes",
      Original_Louvain == "15" ~ "T_Cells B",
      Original_Louvain == "16" ~ "NK_Cells B",
      Original_Louvain == "17" ~ "Neutrophils" 
    ), 
    new_Original_Louvain = Original_Louvain,
    new_subtype = paste0(Original_Louvain," : ",subtype),
         new_subtype = factor(new_subtype, levels = c(
           "0 : MØ (H2-Ab1+)","1 : Monocytes","2 : MØ (Trem2+ A)","3 : MØ (Trem3+ B)","4 : MØ (CD14+)",
           "5 : MØ (Lars+2 A)", "6 : MØ (Sparc+ A)","7 : T_Cells A","8 : intMono/MØ","9 : MØ (Ly6c2+)",
           "10 : MØ (Sparc+ B)","11 : Dendritic ?","12 : NK_Cells B","13 : MØ (Phagocytic)",
           "14 : Pericytes","15 : T_Cells B","16 : NK_Cells B","17 : Neutrophils"
         )))

qsave(x=macrophages,file="/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier/macrophages_f.qs")
macrophages <- qs::qread(file="/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier/macrophages_f.qs")

fig5c = make_fig5(macrophages,"/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier/fig5c.pdf")
fig_macro = make_figS5(macrophages,genes_macro,log = T)

ggsave(plot = fig_macro, path = "/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier/",
       filename="figS5ab_ter.pdf",width = 20, height = 10, dpi=1000,units="in",
       device = cairo_pdf)

vlnplot_genes_macro = lapply(genes_macro,function(x) VlnPlot(macrophages, features=x,group.by = "new_subtype"))

names(macrophages@reductions) = c("pca","harmony","umap","fitsne")
library(wesanderson)
fig5f = plot_imbalance_score_v3(data_seurat = macrophages,reduceDimension_selected = "fitsne",cl = "sample",genes = F,split = F)

fig5f_final = ggplot(fig5f$imbalance_scoring %>% arrange(scores_kNN), aes(x=FItSNE_1,y=FItSNE_2,color=scores_kNN)) +
  geom_point() +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1"),
                        breaks = c(quantile(fig5f$imbalance_scoring$scores_kNN,c(0,.9,1))),
                        limits = c(0,max(fig5f$imbalance_scoring$scores_kNN))) +
  labs(color="Imbalance Score - Smooth kNN") +
  theme_legrand2 +
  theme(legend.position = "bottom")
ggsave(plot=fig5f_final, path = "/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier/",
       filename="fig5f.pdf",width = 10, height = 10, dpi=1000,units="in",
       device = cairo_pdf)

##after slingshot_v2
# plot_imbalance_score_gene(merged_seurat,"fitsne","Fcer1g") + p1$canonical_plots$clustering$Original_Louvain
# plot_imbalance_score_gene(merged_seurat,"fitsne","Ncr1") + p1$canonical_plots$clustering$Original_Louvain
# plot_imbalance_score_gene(merged_seurat,"fitsne","Ms4a1",smooth = 20) + p1$canonical_plots$clustering$Original_Louvain
# plot_imbalance_score_gene(merged_seurat,"fitsne","Itgam") + p1$canonical_plots$clustering$Original_Louvain
# plot_imbalance_score_gene(merged_seurat,"fitsne","Adgre1") + p1$canonical_plots$clustering$Original_Louvain
# plot_imbalance_score_gene(merged_seurat,"fitsne","Ptprc") + p1$canonical_plots$clustering$Original_Louvain
# 
# plot_imbalance_score_v2(merged_seurat,"fitsne",x = "sample",split = F) + theme_legrand2 +
#   + p1$canonical_plots$clustering_splitted$predicted.id
# 
# FeaturePlot(merged_seurat,c("Lars2","H2-Ab1","Trem2","Sparc","Cd14","Ly6c2"),cols=c("blue","orange","red"),pt.size=1,order=T,reduction="fitsne") |
#   p1$canonical_plots$clustering$Original_Louvain

mtdt <- extract_metadata_reduction_seurat(macrophages)[,c(91:140,8,164)] %>%
  arrange(sample,new_subtype) %>%
  rename("Sample"="sample",
         "Clusters"="new_subtype")
colnames(mtdt) = str_replace_all(colnames(mtdt),"_kNN","")

mtdt = mutate(mtdt,Clusters=factor(Clusters))

ca = columnAnnotation(df =  mtdt[,c("Clusters","Sample")],
                      annotation_legend_param = list(
                        Sample = list(
                          labels=expression(Setdb1^{MuSC-KO},
                                            Control)
                        )
                      ),
                      col = list(
                        Sample = c("KO"="red","WT"="blue"),
                        Clusters = setNames(fig5c$palette$new_subtype,levels(mtdt$Clusters))
                      )
                    )


ht = Heatmap(t(mtdt[,c(1:50)]),
        name = "UCell scores - Smooth KNN",
        # col = circlize::colorRamp2(c(0,.05,.4), c("blue", "white", "red")),
        use_raster = F,
        top_annotation = ca,
        show_row_names = T,
        row_names_gp = gpar(fontface="bold.italic",fontsize="11",family="Helvetica"),
        cluster_rows = T,
        show_row_dend = T,
        row_km = 3,
        column_km = 3,
        cluster_columns = F,
        show_column_dend = F,
        show_column_names = F,
        heatmap_legend_param = list(
          direction="horizontal",nrow=1
        ))
pdf(file = "/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier/figS5f.pdf",width = 10,height = 10)
draw(ht)
dev.off()

tmp = extract_metadata_reduction_seurat(macrophages) %>%
  mutate(sample = ifelse(sample == "KO", "Setdb1^{MuSC-KO}","Control") %>% factor())

alphaP=1
sizeP=2
sizePtext=0.36*20
sizetext=15
family="Helvetica"

p=ggplot() +
  geom_point(data=dplyr::select(tmp,-sample), aes(FItSNE_1,FItSNE_2),color="grey80") +
  geom_point(data=tmp,aes(x=FItSNE_1,y=FItSNE_2,color=sample),size=2) +
  scale_color_manual(values=c("blue","red"),
                     labels=expression(Control,
                                       Setdb1^{MuSC-KO})) +
  labs(color="Sample") +
  facet_wrap(~sample,labeller = as_labeller(labeller_modif(tmp,"sample"),label_parsed)) +
  theme_legrand2

ggsave(plot = p, path = "/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier/",
       filename="fig_samples_macrophages.pdf",width = 10, height = 7, dpi=1000,units="in",
       device = cairo_pdf)

palettes = qs::qread("/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier/v3.2.3/palette.qs")
colSums_total = dim(extract_metadata_reduction_seurat(macrophages))[1]
a =  extract_metadata_reduction_seurat(macrophages) %>% group_by(Original_Louvain) %>% summarise(count=n()) %>% arrange(count) %>% pull(Original_Louvain)

colsums_percent <- function(x,total) {
  return(x/total*100)
}
p1=DimPlot(macrophages,reduction="fitsne",group.by="new_subtype",split.by = "sample",label=T)
fig2c = extract_metadata_reduction_seurat(macrophages) %>% 
  group_by(Original_Louvain,sample) %>% summarise(count=n()) %>%
  mutate(sample = factor(sample, 
                         levels = unique(macrophages$sample)),
         Original_Louvain = factor(Original_Louvain, levels=a)) %>% 
  pivot_wider(id_cols="Original_Louvain",names_from = "sample",values_from = "count") %>%
  summarise(across(where(is.numeric),~colsums_percent(.x,total=colSums_total))) %>%
  pivot_longer(cols = where(is.numeric),names_to = "sample",values_to = "vRatio")

write.csv(fig2c,file = "/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier/proportioncells_by_clusters_macrophages.csv",row.names = F)

maxvratio <- fig2c %>% dplyr::select(Original_Louvain,vRatio) %>% group_by(Original_Louvain) %>% summarise(max=sum(vRatio)) %>% pull(max) %>% max() %>% round(-1)
label_ratio = setNames(paste0(seq(0,maxvratio,10),"%"),seq(0,maxvratio,10))
label_x = extract_metadata_reduction_seurat(macrophages) %>% 
  dplyr::select(Original_Louvain,new_subtype) %>% unique() %>% 
  mutate(Original_Louvain = factor(Original_Louvain, levels=levels(fig2c$Original_Louvain))) %>%
  arrange(Original_Louvain) %>% 
  mutate(new_subtype=as.character(new_subtype)) %>% pull(new_subtype)
fig2c = fig2c %>%
  arrange(Original_Louvain) %>%
  mutate(Original_Louvain = factor(Original_Louvain,labels=label_x)) %>%
  ggplot(aes(x=Original_Louvain,y=vRatio,fill=sample)) + 
  geom_col() +
  scale_fill_manual(values=c("blue","red"),
                     labels=expression(Control,
                                       Setdb1^{MuSC-KO})) +
  # scale_fill_manual(values=palettes$sample,labels=scales::parse_format()) +
  # scale_y_continuous(limits = c(0,30),labels=seq(0,max(label_ratio))) +
  scale_y_continuous(limits = as.numeric(c(min(names(label_ratio)),max(names(label_ratio)))),
                     breaks=as.numeric(names(label_ratio)), 
                     labels = label_ratio) +
  labs(x="",y="Relative Cells\nAbundance",fill="Sample") +
  theme_legrand2 +
  theme(axis.text = element_text(face="bold.italic",family="Helvetica",size=11)) +
  coord_flip() +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(ncol=2))
ggsave(plot = fig2c, path = "/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/figures_papier/",
       filename="barplot_samples_macrophages_by_subtype.pdf",width = 10, height = 7, dpi=1000,units="in",
       device = cairo_pdf)
