# Load libraries 
library(dplyr)
library(tidyverse)
library(edgeR)
library(lemon)
library(ggvenn)
library(ComplexHeatmap)
library(viridis)
library(ggpubr)
library(gridExtra)

setwd("/home/alejandra/Documentos/Cranberry_transcriptomics/") #Set working directory

count_data <- read.table("Counts/counts.cranberry.txt", header = T) #Load counts table
genetic_data <- count_data %>% dplyr::select(Geneid:Length) #Keep genetic information

rownames(count_data) <- NULL
count_data <- count_data %>% 
  dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>% #Selecting only sample columns 
  mutate(Geneid=gsub("gene:|_.*", "",Geneid)) %>% #Shorting gene names
  column_to_rownames("Geneid") %>% #Gene names are now the rownames of the count table 
  dplyr::rename("C1.R1.3" = R1.3.2, "C1.R4.6"=R4.6.2, "C1.R7.9"=R7.9.2, #Rename samples name
                "C2.R1.3" = R1.3.4, "C2.R4.6"=R4.6.4, "C2.R7.9"=R7.9.4,
                "C3.R1.3" = R1.3.6, "C3.R4.6"=R4.6.6, "C3.R7.9"=R7.9.6, 
                "C1.V1.3" = VI.3.2, "C1.V4.6"=V4.6.2, "C1.V7.9"=V7.9.2,
                "C2.V1.3" = V1.3.4, "C2.V4.6"=V4.6.4, "C2.V7.9"=V7.9.4,
                "C3.V1.3" = V1.3.6, "C3.V4.6"=V4.6.6, "C3.V7.9"=V7.9.6) %>%
  dplyr::select(matches("C1.R"), matches("C2.R"), matches("C3.R"), #Reordering columns 
                matches("C1.V"), matches("C2.V"), matches("C3.V"))
head(count_data)

dim(count_data) #Of 22836 annotated gene models
table(rowSums(count_data) == 0) #21705 have counts in at least one sample.


## EdgeR analysis ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

group <- gsub("(R|V).*$", "\\1", colnames(count_data)) #Grouping factor

y <- DGEList(counts=count_data, group = group, genes = rownames(count_data)) #Store data in a DGEList object
y$samples
y$samples$group

keep <- filterByExpr(y) #filter low expressed genes
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y) #19075 genes after filtering low expressed genes 

y <- calcNormFactors(y) #Calculate normalization factors
y <- estimateDisp(y) #Estimate the dispersion and tagwise dispersion 

#Explore data
plotMDS(y, gene.selection = "common")
# plotMDS(y, pch = 16, gene.selection = "common")

design <- model.matrix(~0+group, data=y$genes)
colnames(design) <- levels(y$samples$group)
fit <- glmQLFit(y, design)

results <- list() #List in blank to store the contrasts results 

for (i in c("C1","C2","C3")) {
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#          Comparisons between R and V bud tissue at the same stage         :
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   ## Make contrasts
      a <-paste(i,"R", sep = ".")
      b <- paste(i, "V", sep = ".")
      contrast <- makeContrasts(paste(a, b, sep = "-"), levels=design)
      results[[paste(i, "RvsV", sep = "_")]] <-  glmQLFTest(fit, contrast=contrast)
      results[[paste(i, "RvsV", sep = "_")]]$table <- 
        results[[paste(i, "RvsV", sep = "_")]]$table %>% 
        dplyr::mutate(FDR=p.adjust(PValue, method="BH"))
      
   ## How many DEG are in each comparison (pad < 0.05 and log Fold Change value > 1)
      print(paste(
      paste(i, "RvsV DEG", sep = " "),  nrow(results[[paste(i,"RvsV", sep = "_")]]$table %>%
                                               dplyr::filter(FDR < 0.05 & abs(logFC) > 1))
      )) 
}

for (i in c("R","V")) {
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#         Comparisons between across stages on same bud tissue     :
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
     ## Make contrasts
        #C1 vs C2 contrast
        temp <- paste(paste("C1",i, sep = "."), paste("C2", i, sep = "."), sep = "-")
        contrast <- makeContrasts(temp, levels=design)
        results[[paste(i,"C1vsC2", sep = "_")]] <-  glmQLFTest(fit, contrast=contrast)
        
        results[[paste(i, "C1vsC2", sep = "_")]]$table <- 
          results[[paste(i, "C1vsC2", sep = "_")]]$table %>% 
          dplyr::mutate(FDR=p.adjust(PValue, method="BH"))

        #C2 vs C3 contrast
        temp <- paste(paste("C2",i, sep = "."), paste("C3", i, sep = "."), sep = "-")
        contrast <- makeContrasts(temp,levels=design)
        results[[paste(i,"C2vsC3", sep = "_")]] <-  glmQLFTest(fit, contrast=contrast)
        
        results[[paste(i, "C2vsC3", sep = "_")]]$table <- 
          results[[paste(i, "C2vsC3", sep = "_")]]$table %>% 
          dplyr::mutate(FDR=p.adjust(PValue, method="BH"))
    
        #C1 vs C2 contrast
        temp <- paste(paste("C1",i, sep = "."), paste("C3", i, sep = "."), sep = "-")
        contrast <- makeContrasts(temp, levels=design)
        results[[paste(i,"C1vsC3", sep = "_")]] <-  glmQLFTest(fit, contrast=contrast)
        
        results[[paste(i, "C1vsC3", sep = "_")]]$table <- 
          results[[paste(i, "C1vsC3", sep = "_")]]$table %>% 
          dplyr::mutate(FDR=p.adjust(PValue, method="BH"))
        
     ## How many DEG are in each comparison (pad < 0.05 and log Fold Change value > 1)
        print(paste( 
          paste(i, "C1 vs C2 DEG", sep = " "), nrow(results[[paste(i,"C1vsC2", sep = "_")]]$table %>%
                                                      dplyr::filter(FDR < 0.05 & abs(logFC) > 1)), sep = ":"
        )) #How many DEG are in C1 vs C2
        
        print(paste(
          paste(i, "C2 vs C3 DEG", sep = " "), nrow(results[[paste(i,"C2vsC3", sep = "_")]]$table %>%
                                                      dplyr::filter(FDR < 0.05 & abs(logFC) > 1)), sep = ":"
        )) #How many DEG are in C2 vs C3
        
        print(paste(
          paste(i, "C1 vs C3 DEG", sep = " "), nrow(results[[paste(i,"C1vsC3", sep = "_")]]$table %>%
                                                      dplyr::filter(FDR < 0.05 & abs(logFC) > 1)), sep = ":"
        )) #How many DEG are in C1 vs C3
}

for (i in names(results)) {
#::::::::::::::::::::::::::
#         Plot data       :
#::::::::::::::::::::::::::
    results[[i]]$table <- 
      results[[i]]$table %>% 
      dplyr::mutate(status=ifelse(FDR < 0.05 & abs(logFC) > 1, "Significant DEG","No significant DEG")) 
 ## MA plots
    assign(paste("ma",i, sep = "_"), 
           ggplot(results[[i]]$table, aes(y=logFC, x=logCPM, color=as.factor(status))) +
             geom_point() +
             scale_color_manual(name = "",
                                values = c("No significant DEG"="#929ea1","Significant DEG"= "darkred")) +
             geom_hline(yintercept= 0, linetype="dashed", color = "black") +
             labs(title = gsub("_"," ",i), sep = " ")  + 
             theme(plot.tag = element_text(face =  "bold")) + theme_classic()
           )
 ## Volcano plots
    assign(paste("vol", i, sep = "_"), 
           ggplot(results[[i]]$table, aes(x=logFC, y=-log10(FDR), color=as.factor(status))) +
             geom_point() +
             scale_color_manual(name = "",
                                values = c("No significant DEG"="#929ea1","Significant DEG"= "darkred")) +
             geom_hline(yintercept= -log10(0.05), linetype="dashed", color = "black") +
             geom_vline(xintercept= c(1,-1), linetype="dashed", color = "black") +
             labs(title = gsub("_"," ",i), sep = " ") + 
             theme(plot.tag = element_text(face =  "bold")) + theme_classic()
           )
  ## Write results
    # write.table(results[[i]], paste("DEG_results/",".txt", sep = i) , row.names = F, col.names = T)
    
}

## Arrange plots ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Same stage between reproductive and vegetative buds
grid_arrange_shared_legend(vol_C1_RvsV + labs(tag = "A"), ma_C1_RvsV + labs(tag = "B"),
                           vol_C2_RvsV + labs(tag = "C"), ma_C2_RvsV + labs(tag = "D"),
                           vol_C3_RvsV + labs(tag = "E"), ma_C3_RvsV + labs(tag = "F"),
                           ncol=2, nrow=3)

#Between samples of reproductive buds
grid_arrange_shared_legend(vol_R_C1vsC2 + labs(tag = "A"), ma_R_C1vsC2 + labs(tag = "B"),
                           vol_R_C2vsC3 + labs(tag = "C"), ma_R_C2vsC3 + labs(tag = "D"),
                           vol_R_C1vsC3 + labs(tag = "E"), ma_R_C1vsC3 + labs(tag = "F"),
                           ncol = 2, nrow = 3)

#Between samples of vegetative buds
grid_arrange_shared_legend(vol_V_C1vsC2 + labs(tag = "A"), ma_V_C1vsC2 + labs(tag = "B"),
                           vol_V_C2vsC3 + labs(tag = "C"), ma_V_C2vsC3 + labs(tag = "D"),
                           vol_V_C1vsC3 + labs(tag = "E"), ma_V_C1vsC3 + labs(tag = "F"),
                           ncol = 2, nrow = 3)

## Bar plot ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

for (i in c("C1","C2","C3")) {
  assign(i,results[[paste(i, "RvsV", sep = "_")]]$table %>% 
           dplyr::filter(FDR < 0.05 & abs(logFC) > 1) %>%
           dplyr::mutate(deg=ifelse(FDR < 0.05 & logFC > 1, "Up","Down")) %>% 
           dplyr::mutate(stage = paste(i, "RvsV", sep = "_")) %>% dplyr::select(deg,stage))
}

temp <- nrow(C1) + nrow(C2) + nrow(C3) 
temp2 <- c(rep("C1_RvsV", times=nrow(C1)), rep("C2_RvsV", times=nrow(C2)),
           rep("C3_RvsV", times=nrow(C3)))

all <- data.frame(deg=rep("All", times=temp), stage=temp2)

df_res <- rbind(C1, C2, C3, all)

df_res$stage <- factor(df_res$stage, levels = c("C1_RvsV","C2_RvsV","C3_RvsV"))
df_res$deg <- factor(df_res$deg, levels = c("Up","Down","All"))

df <- as.data.frame(table(df_res))

ggplot(data = df, aes(x = stage, y = Freq, fill = deg)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75) +
  geom_text(aes(label = Freq), vjust = -0.5, #fontface = "bold"
            position = position_dodge(.9), size = 4) + 
  labs(x = "", y = "DEG gene counts\n", ) +
  theme_classic() + 
  theme(legend.title = element_blank(), legend.position = c(0.1,0.9), 
        legend.background = element_rect(linetype = "solid", color = "black",
                                         size=0.5), 
        legend.text = element_text(size = 10),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=11, face = "bold"),
        axis.text.y = element_text(size=11))

## Venn diagram :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Create logical table
venn_df<- tibble(count_data %>% 
                   rownames_to_column(., "geneid")) %>%
                   dplyr::select("geneid") %>%
                   dplyr::mutate(C1_RvsV= geneid %in% rownames(C1)) %>%
                   dplyr::mutate(C2_RvsV= geneid %in% rownames(C2)) %>%
                   dplyr::mutate(C3_RvsV= geneid %in% rownames(C3))


ggplot(venn_df) +
  geom_venn(aes(A = C1_RvsV, B = C2_RvsV, C= C3_RvsV), 
            fill_color = c("#F8766D", "#00BA38", "#619CFF"), 
            stroke_color = "#302e2e", 
            stroke_size = 0.5,
            fill_alpha = 0.6, 
            set_name_size = 6, 
            text_size = 6, 
            text_color = "#1a1918", 
            set_names = c("C1 RvsV","C2 RvsV", "C3 RvsV")) +
  coord_fixed(xlim = c(-1.6,1.6), ylim = c(-1.8,1.8)) +
  theme_void() +
  theme(plot.tag = element_text(size = 15, face = "bold", vjust = -7, hjust = 1),
        plot.title=element_text(size=15,hjust=0.5, vjust = -5)) 

## Heatmap
deg_genes <- venn_df %>% dplyr::filter(C1_RvsV == T | C2_RvsV == T | C3_RvsV == T) %>% dplyr::select(geneid) 
nrow(deg_genes) #2548 DEG in at least one comparison.

#Normalization
#Calculate cpm counts
dge <- DGEList(count_data)
dge <- calcNormFactors(dge, method="TMM")
counts_cpm <- cpm(dge, normalized.lib.sizes = T,log=FALSE) 
#check normalization
gather_raw <- gather(log2(count_data+1))
gather_raw <- gather_raw %>% mutate(Sample=gsub("(R|V).*", "\\1",key))

gather_cpm <- gather(log2(counts_cpm %>% as.data.frame()+1))
gather_cpm <- gather_cpm %>% mutate(Sample=gsub("(R|V).*", "\\1",key))

list_gather <- list(raw=gather_raw, cpm=gather_cpm)

for (i in c("raw","cpm")) {
  #Density plot
  assign(
    paste("density", i, sep = "_"),
    ggplot(list_gather[[i]]) +
      aes(x=value, fill=key) +
      geom_density(alpha=.3) + 
      labs(y="Density", x=paste(paste("log2(", i, sep =  ""), "counts +1)", sep = " "), fill="")
          )
  assign(
    paste("R", i, sep = "_"),
    ggplot(list_gather[[i]] %>% dplyr::filter(endsWith(Sample, 'R')), 
           aes(x=value, fill=Sample)) + geom_density(alpha=.3) + 
      labs(y="Density", x=paste(paste("log2(", i, sep =  ""), "counts +1)", sep = " "), fill="")
      )
  
  assign(
    paste("V", i, sep = "_"),
    ggplot(list_gather[[i]] %>% dplyr::filter(endsWith(Sample, 'V')), 
           aes(x=value, fill=Sample)) + geom_density(alpha=.3) + 
      labs(y="Density", x=paste(paste("log2(", i, sep =  ""), "counts +1)", sep = " "), fill="")
  )
  
  # Blox plot
  assign(
    paste("vp", i, sep = "_"),
    ggplot(list_gather[[i]], aes(x=key, y=value, fill=Sample)) + 
      geom_boxplot() +
      labs(y= "log2(cpm counts +1)", x="", title = paste(i,"counts", sep = " ")) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
         )
  
}

ggarrange(density_raw, density_cpm, ncol=2, common.legend = T)
ggarrange(R_raw, V_raw, ncol = 2, labels = c("A","B"))
ggarrange(R_cpm, V_cpm, ncol = 2,  labels = c("A","B"))
ggarrange(vp_raw, vp_cpm, ncol = 2, common.legend = T)


#Filter only DEG genes in at least one contrast
deg_cpm_counts <- counts_cpm %>% 
  as.data.frame() %>%
  dplyr::filter(rownames(.) %in% deg_genes$geneid) 

dim(deg_cpm_counts) #2548 
mat_scaled = t(scale(t(deg_cpm_counts)))
colnames(mat_scaled)

ann <- data.frame(Stage=rep(c("C1","C2","C3"), each=3), 
                  Buds = rep(c("Reproductive","Vegetative"), each=9))
ann$Buds <- factor(ann$Buds, levels = c("Reproductive","Vegetative"))
ann <- ann[,c("Buds", "Stage")]

colAnn <- HeatmapAnnotation(df=ann,
                            show_annotation_name = F,
                            which = 'col',
                            na_col = 'black',
                            col = list(Buds=c(Reproductive="cadetblue4",Vegetative="cadetblue3"),
                                       Stage=c(C1="burlywood4",C2="darkgoldenrod3",C3="mediumorchid4")),
                            annotation_height = 0.6,
                            height = 2)
meth_col_fun = magma(256)

hm <-  Heatmap(mat_scaled,
               cluster_row_slices = F,
               name="Expression",
               col = meth_col_fun,
               # split = split,
               na_col = "black",
               row_gap = unit(1.6, "mm"),
               cluster_rows = T,
               show_row_dend = T,
               show_row_names = F,
               border = T,
               width = unit(9, "cm"), height = unit(16, "cm"),
               cluster_columns = F,
               show_column_dend = F,
               row_names_gp = gpar(fontsize = 8),
               top_annotation = colAnn) 
draw(hm)

#Save session information
sink("session.info.txt")
sessionInfo()
sink()
