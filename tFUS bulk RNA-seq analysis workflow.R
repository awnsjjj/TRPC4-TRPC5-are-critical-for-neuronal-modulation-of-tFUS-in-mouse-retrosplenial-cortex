Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)


###Package installation
install.packages("BiocManager")
install.packages("remotes")
library(remotes)
library(BiocManager)
install.packages("org.Mm.eg.db")
install.packages("ggplot2")
install.packages("vctrs")
install.packages("UpSetR")
install.packages("reshape2")
install.packages("data.table")
BiocManager::install(c("DEGseq", "qvalue"))
install_github("yanlinlin82/ggvenn")

#####Differential analysis
library(qvalue)
library(DEGseq)
library(data.table)
geneExpFile <-data.frame(fread("gene_count_matrix.csv",header=T))
geneExpFile$name<-make.unique(geneExpFile$name)
row.names(geneExpFile) <-geneExpFile$name
geneExpFile1<-geneExpFile[3:14]
geneExpFile1<-geneExpFile1+1
write.csv(geneExpFile1,file="gene+count+matrix.csv")
treatment_positive_df <-readGeneExp(file = "gene+count+matrix.csv",geneCol = 1,valCol = c(2,4,6),sep = ",")
control_positive_df <-readGeneExp(file = "gene+count+matrix.csv",geneCol = 1,valCol = c(8,10,12),sep = ",")
treatment_df <-readGeneExp(file = "gene+count+matrix.csv",geneCol = 1,valCol = c(3,5,7),sep = ",")
control_df <-readGeneExp(file = "gene+count+matrix.csv",geneCol = 1,valCol = c(9,11,13),sep = ",")
treatment_positive_depth <-c(48359508,52611832,50916518)
control_positive_depth <-c(41714122,48820574,49144424)
treatment_depth<-c(49698164,42923092,63185328)
control_depth <-c(50013002,47967490,47360832)
#options(digits = 22)
#Comparison within the control group
DEGexp(geneExpMatrix2 = control_positive_df, geneCol1 = 1, expCol1 = 2:4, depth1 = control_positive_depth, groupLabel1 = "Tre",
       geneExpMatrix1 = control_df, geneCol2 = 1, expCol2 = 2:4, depth2 = control_depth, groupLabel2 = "Con",
       method = "LRT", normalMethod = "median", outputDir = "tmp")

#Comparison within the tFus group
DEGexp(geneExpMatrix2 = treatment_positive_df, geneCol1 = 1, expCol1 = 2:4, depth1 = treatment_positive_depth, groupLabel1 = "Tre",
       geneExpMatrix1 = treatment_df, geneCol2 = 1, expCol2 = 2:4, depth2 = treatment_depth, groupLabel2 = "Con",
       method = "LRT", normalMethod = "median", outputDir = "tmp")
#Comparison between the control and the tFus group
DEGexp(geneExpMatrix2 = treatment_positive_df, geneCol1 = 1, expCol1 = 2:4, depth1 = treatment_positive_depth, groupLabel1 = "Tre",
       geneExpMatrix1 = control_positive_df, geneCol2 = 1, expCol2 = 2:4, depth2 = control_positive_depth, groupLabel2 = "Con",
       method = "LRT", normalMethod = "median", outputDir = "tmp")
#Rename each output txt file and saved as csv format
########vocalno painting
library(ggplot2)

contr_compare<-data.frame(fread("control+G vs -G adj.csv",header=T))
Tfus_compare<-data.frame(fread("tfus+G vs -G adj.csv",header=T))
Tfus_contr_compare<-data.frame(fread("tfus+G vs control +G adj.csv",header=T))
M<-contr_compare
M<-Tfus_compare
M<-Tfus_contr_compare
logFC_cutoff <-2
Pvalue_cutoff<-0.001
#Let each comparison data equals to M and repeat the following codes
M$change[M$log2.Fold_change.>=logFC_cutoff & M$q.value.Benjamini.et.al..1995.<Pvalue_cutoff] <-"Down"
M$change[(M$log2.Fold_change.<logFC_cutoff &
            M$log2.Fold_change.>-logFC_cutoff)
         | M$q.value.Benjamini.et.al..1995.>Pvalue_cutoff] <-"No"
M$change[M$log2.Fold_change.<=-logFC_cutoff & M$q.value.Benjamini.et.al..1995.<Pvalue_cutoff] <-"Up"
table(M$change)
M$change <-as.factor(M$change)
p<-ggplot(M,aes(x=log2.Fold_change.,y=-log10(q.value.Benjamini.et.al..1995.),color=change))+
  geom_point(alpha=0.4,size=1)+
  scale_color_manual(values = c("Down"='#006699',"No"='#bebebe',"Up"='#ffad21'))+
  geom_vline(xintercept = c(-logFC_cutoff,logFC_cutoff),linetype="dashed",color="black",linewidth=1)+
  geom_hline(yintercept=-log10(Pvalue_cutoff),linetype="dashed",color="black",linewidth=1)+
  labs(x="log2(Fold Change)",y="-log10(P Value)")+
  theme_bw()+
  theme(legend.position = "right")
p1<-p+scale_y_continuous(limits = c(0,20))
ggsave('vocalno.pdf',plot = p1,width = 8,height = 6)
contr_compare<-merge(M,geneExpFile1[,7:12],by.x = 1,by.y = 0)
Tfus_compare<-merge(M,geneExpFile1[,1:6],by.x = 1,by.y = 0)
Tfus_contr_compare<-merge(M,geneExpFile1[,c(1,3,5,7,9,11)],by.x = 1,by.y = 0)
write.csv(contr_compare,file="contr_compare.csv")
write.csv(Tfus_compare,file="Tfus_compare.csv")
write.csv(Tfus_contr_compare,file="Tfus_contr_compare.csv")

############Venn painting
library(grid)
library(vctrs)
library(ggvenn)
conUp<-data.frame(fread("conUp.csv",header=T))
tfusUp<-data.frame(fread("tfusUp.csv",header=T))
tvcUp<-data.frame(fread("tvcUp.csv",header=T))

x <-list('Control'=conUp$GeneNames,
         'Tfus vs Control'=tvcUp$GeneNames,
         'Tfus'=tfusUp$GeneNames)
Up<-ggvenn(x,
           show_percentage = F,
           stroke_color = "white",
           fill_color = c("#b2e7cb","#b2d4ec","#ffb2b2"),
           set_name_color = c("#4a9b83","#1d6295","#ff0000"))
ggsave('venn_Up.pdf',plot = Up,width = 8,height = 6)

conDown<-data.frame(fread("conDown.csv",header=T))
tfusDown<-data.frame(fread("tfusDown.csv",header=T))
tvcDown<-data.frame(fread("tvcDown.csv",header=T))
x <-list('Control'=conDown$GeneNames,
         'Tfus vs Control'=tvcDown$GeneNames,
         'Tfus'=tfusDown$GeneNames)
Down<-ggvenn(x,
             show_percentage = F,
             stroke_color = "white",
             fill_color = c("#b2e7cb","#b2d4ec","#ffb2b2"),
             set_name_color = c("#4a9b83","#1d6295","#ff0000"))
ggsave('venn_Down.pdf',plot = Down,width = 8,height = 6)

#########intersection filtering
df1<-tvcUp$GeneNames
df2<-conUp$GeneNames
df3<-tfusUp$GeneNames
tvcup_re<-setdiff(df1,df2)
tfusup_re<-setdiff(df3,df2)
tvc_plus_tfus<-intersect(tvcup_re,tfusup_re)

#############GO analysis
library(org.Mm.eg.db)
library(clusterProfiler)

gene1<-tvc_plus_tfus
gene_ENTREZID <- unlist(na.omit(mapIds(x = org.Mm.eg.db,
                                       keys =  gene1,
                                       keytype = "SYMBOL",
                                       column = "ENTREZID",
                                       multiVals = "first")))

go_enrich_results_ALL <- enrichGO(gene = gene_ENTREZID,
                                  OrgDb = "org.Mm.eg.db",
                                  ont   = "ALL"  ,    
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.05,
                                  readable      = TRUE)
write.csv(go_enrich_results_ALL@result, 'GO_gene_ALL_enrichresults.csv') 

#####KEGG analysis
kegg_enrich_results <- enrichKEGG(gene  = gene_ENTREZID,
                                  organism  = "mmu",
                                  keyType = "kegg",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.05)
kk_read <- DOSE::setReadable(kegg_enrich_results, 
                             OrgDb="org.Mm.eg.db", 
                            keyType='ENTREZID')#ENTREZID to gene Symbol
write.csv(kk_read@result,'KEGG_gene_enrichresults.csv') 

########GO terms painting
library(enrichplot)

bp <-data.frame(fread("GO_gene_ALL_enrichresults_tfusplustvcBP.csv",header = T))
mf <-data.frame(fread("GO_gene_ALL_enrichresults_tfusplustvcMF.csv",header = T))
cc <-data.frame(fread("GO_gene_ALL_enrichresults_tfusplustvcCC.csv",header = T))

d<-bp
d<-mf
d<-cc
#Let each GO ontology equals to d and repeat the following codes
d<-d[order(-d$Count),]
d <-d[1:13,]
d$Description <- factor(d$Description,levels=d$Description)

mytheme <- theme(axis.title=element_text(face="bold", size=14,colour = 'black'), 
                 axis.text.y =element_text(face="bold", size=14,colour = 'black'), 
                 axis.text.x=element_text(size=8),
                 axis.line = element_line(linewidth=0.5, colour = 'black'), 
                 panel.background = element_rect(color='black'), 
                 legend.key = element_blank() 
)
p <- ggplot(d,aes(x=Count,y=Description,colour=-1*log10(p.adjust),size=Count))+
  geom_point()+
  scale_size(range=c(2, 8))+
  scale_colour_gradient(low = "blue",high = "red")+
  theme_bw()+
  ylab("GO_BP Pathway Terms")+
  xlab("Gene numbers")+
  labs(color=expression(-log[10](PValue)))+mytheme
ggsave('GO_BP.pdf',plot = p,width = 10,height = 6)

#########KEGG terms painting
kegg_up <-data.frame(fread("KEGG_gene_enrichresults_tfusplustvc.csv",header = T))
kegg_down <-data.frame(fread("KEGG_gene_enrichresults_tfusdowntvc.csv",header = T))
d<-kegg_up
d<-d[order(-d$Count),]
d <-d[1:10,]
d$Description <- factor(d$Description,levels=d$Description)
kegg_up<-d
d<-kegg_down
d<-d[order(-d$Count),]
d <-d[1:10,]
d$Description <- factor(d$Description,levels=d$Description)
kegg_down<-d
kegg<-rbind(kegg_up,kegg_down)
kegg$number <- factor(rev(1:nrow(kegg)))
kegg$type<-factor(c(rep("Up", 10),rep("Down", 10)),levels=c("Up", "Down"))

p <- ggplot(data=kegg, aes(x=number, y=Count, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  scale_fill_manual(values = c('#FD8D62',"#8DA1CB")) + theme_test() + 
  scale_x_discrete(labels=kegg$Description) +
  xlab("KEGG term") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched KEGG Terms")
ggsave('kegg_all.pdf',plot = p,width = 10,height = 8)

################preparation file for cytoscape and Upset of GO terms
library(UpSetR)

go_results<- data.frame(fread("pathway_select.csv",header = T))
colnames(go_results)[colnames(go_results)=="Type"]<-"Ontology"

nodes_list <- list()
edges_list <- list()
for (i in 1:nrow(go_results)) {
  pathway <- go_results[i, ]
  genes <- unlist(strsplit(pathway$geneID, "/"))
  pathway_node <- data.frame(
    ID = pathway$ID,
    name = pathway$Description,
    Type = "Pathway",
    Ontology = pathway$Ontology,
    stringsAsFactors = FALSE
  )
  nodes_list[[i]] <- pathway_node
  
  gene_nodes <- data.frame(
    ID = genes,
    name = genes,
    Type = "Gene",
    Ontology = NA,  
    stringsAsFactors = FALSE
  )
  nodes_list[[i]] <- rbind(nodes_list[[i]], gene_nodes)
  
  pathway_edges <- data.frame(
    fromNode_ID = rep(pathway$ID, length(genes)),
    fromNode_name = rep(pathway$Description, length(genes)),
    toNode_ID = genes,
    toNode_name = genes,
    Relations_Type = rep("Path-gene", length(genes)),
    Ontology = rep(pathway$Ontology, length(genes)),
    stringsAsFactors = FALSE
  )
  edges_list[[i]] <- pathway_edges
}

nodes <- do.call(rbind, nodes_list)
edges <- do.call(rbind, edges_list)

gene_freq<-table(edges$toNode_name)
gene_freq_df <- as.data.frame(gene_freq)
colnames(gene_freq_df) <- c("Gene", "Frequency")
gene_freq_df <- gene_freq_df[order(-gene_freq_df$Frequency), ]
select<-gene_freq_df[gene_freq_df$Frequency>=5,]
p<-ggplot(select, aes(x = Gene, y = Frequency)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Gene Frequency in GO Pathways", x = "Gene", y = "Frequency")
ggsave('gene_frequency_GO.pdf',plot =p,width = 10,height=6 )
write.table(nodes, "nodes.txt", sep = "\t", row.names = FALSE)
write.table(edges, "edges.txt", sep = "\t", row.names = FALSE)

x <-list('GO:0042391'=edges_list[[1]]$toNode_ID,
         'GO:1990351'=edges_list[[2]]$toNode_ID,
         'GO:1902495'=edges_list[[3]]$toNode_ID,
         'GO:0022804'=edges_list[[4]]$toNode_ID,
         'GO:0046873'=edges_list[[5]]$toNode_ID,
         'GO:0005216'=edges_list[[6]]$toNode_ID)
p<-upset(fromList(x), 
         nsets = 6,  
         order.by = "freq", 
         mainbar.y.label = "Intersection size", 
         sets.x.label = "Set size", 
         sets.bar.color = c("#b2e7cb","#ff0000","#ff0000","#b2d4ec","#b2d4ec","#ff0000"
         ),  
         matrix.color = "black", 
         main.bar.color = "black",  
         text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1),  
         shade.color = "gray88") 
pdf("upset_plot.pdf", width = 10, height = 8)
print(p)
dev.off()

#####heatmap painting
library(ggplot2)
library(reshape2)
select <- data.frame(fread("filtered_expr.csv",header = T))
row.names(select)<-select$GeneNames
select<-select[,11:16]
select<-log2(select+0.01)

scaled_matrix <- scale(t(select))
row.names(scaled_matrix)<-c('tFUS-1','tFUS-2','tFUS-3','Ctrl-1','Ctrl-2','Ctrl-3')
melted_matrix <- melt(scaled_matrix)

p<-ggplot(melted_matrix, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#008B8B", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x =element_text(angle = 60, vjust = 0.5, hjust = 0.3)) 
ggsave('heat_map_select.pdf',plot = p,width = 12,height = 3)