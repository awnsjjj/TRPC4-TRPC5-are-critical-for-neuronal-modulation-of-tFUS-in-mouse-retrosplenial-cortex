# TRPC4-TRPC5-are-critical-for-neuronal-modulation-of-tFUS-in-mouse-retrosplenial-cortex
Transcranial focused ultrasound (tFUS) non-invasively regulates neuronal activity, yet its underlying molecular mechanisms remain largely unknown. This study identifies TRPC4 and TRPC5 as critical mediators of tFUS-induced neuronal modulation. 
Project Overview
This project performs a bioinformatics analysis pipeline for RNA sequencing data of mouse samples. The pipeline includes  differential expression analysis, functional enrichment analysis, and paiting. The analysis focuses on identifying differentially expressed mRNAs and their functional roles using Gene Ontology (GO) and KEGG pathway enrichment.

Software and Versions
R: Version 4.2.1, for statistical analysis and visualization.
Cytoscape: Version 3.9.1, for network visualization.

Dependencies
R packages: Ensure the following are installed for visualization and analysis:
install.packages(c("qvalue","data.table","ggplot2","clusterProfiler","DEGseq","grid","vctrs","devtools","clusterProfiler","enrichplot","UpSetR","reshape2"))
devtools::install_github("yanlinlin82/ggvenn")
BiocManager::install("org.Mm.eg.db")
R script file:
Bulk RNA-seq analysis workflow.R
External databases:
Gene Ontology (GO): http://www.geneontology.org
KEGG: https://www.genome.jp/kegg/
Cytoscape: Download from https://cytoscape.org/.

Pipeline Steps
1. Differential analysis
Tool: R (Version 4.2.1)
R packages: qvalue, DEGseq, data.table
Input: gene_count_matrix.csv
Description:
Differential expression analysis was performed using the expression matrix to identify significantly differentially expressed mRNAs.
Criteria for significance:
|log2(fold change)| > 2
Adjusted P-value (P_adj) < 0.001
Output: List of differentially expressed mRNAs.

2. Vocalno painting  
Tool: R (Version 4.2.1)
R packages: ggplot2
Input: control+G vs -G adj.csv"/"tfus+G vs -G adj.csv/tfus+G vs control +G adj.csv
Description: 
Make vaocalno painting of differentially expressed mRNAs.
Output: Vocalno plots.

3. Venn painting and intersction filering
Tool: R (Version 4.2.1)
R packages: vctrs, devtools, ggvenn
Input: conUp.csv/tfusUp.csv/tvcUp.csv; conDown/tfusDown/tvcDown
Description: 
Make vaocalno painting of Up regulated or Down regulated DEG. Picking out the intersection gene names of different DEGs.
Output: Venn plots and intersction gene names.

4. Functional Enrichment Analysis
Tool: R (version 4.2.1)
R packages: org.Mm.eg.db, clusterProfiler
Input: intersction gene names obtained from step 3.
Description:
Gene Ontology (GO) and KEGG pathway enrichment analyses were performed.
Output:
Enrichment analysis results (tables).


5. Go terms and KEGG terms visualization
Tools: R (version 4.2.1)
R packages: enrichplot, ggplot2
Input: GO_gene_ALL_enrichresults_tfusplustvcBP.csv/GO_gene_ALL_enrichresults_tfusplustvcMF.csv/GO_gene_ALL_enrichresults_tfusplustvcCC.csv/KEGG_gene_enrichresults_tfusplustvc.csv/KEGG_gene_enrichresults_tfusdowntvc.csv
Description: 
Functional enrichment analysis results visualization.
Output: Go terms plots and KEGG terms plots.

6. Preparation file for cytoscape and Upset of GO terms
Tool: R (version 4.2.1)，Cytoscape (Version 3.9.1)
R packages: UpSetR
Input: pathway_select.csv
Description: 
The gene-pathway network visualization and intersection analysis visualization of intersted genes and pathways and the preparation process of Cytoscape input file. 
Output: Plots of intersection analysis and gene-pathway network.

7. Heatmap paiting
Tool: R (version 4.2.1)
R packages: ggplot2, reshape2
Input: filtered_expr.csv
Description:
Making heatmap of the interested genes from the intersection analysis.
Output:
Heatmap plots.

Usage
1. Prepare Input Data:
File required for each step has been list in the pipeline and can be easily found in the same filefold.

2. Run the Pipeline:
Follow the steps above in sequence. The total running duration is about half an hour.

3. Analyze in R:
Load the expression matrix into R and run differential expression analysis:
library(DEGseq)
#Comparison within the control group
DEGexp(geneExpMatrix2 = control_positive_df, geneCol1 = 1, expCol1 = 2:4, depth1 = control_positive_depth, groupLabel1 = "Tre",
       geneExpMatrix1 = control_df, geneCol2 = 1, expCol2 = 2:4, depth2 = control_depth, groupLabel2 = "Con",
       method = "LRT", normalMethod = "median", outputDir = "tmp")

4. Visualize Networks:
Open Cytoscape and import the Gene-Pathway network data for visualization.


License
This project is licensed under the GNU General Public License (GPL) version 2 or later, consistent with the licensing of R (version 4.2.1).

Contact
For questions or contributions, please contact lixiangy@zju.edu.cn.
