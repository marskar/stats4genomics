Transcriptomics and methylomics of human Monocytes and T cells	
Characteristics	Descriptions
characteristics: raceGenderSite	"Variable representing participant race, gender, and study site"
characteristics: bcell	Enrichment score from GSEA for B-cell contamination
characteristics: mono	Enrichment score from GSEA for monocyte contamination
characteristics: nkcell	Enrichment score from GSEA for natural killer cell contamination
characteristics: neutro	Enrichment score from GSEA for neutrophil contamination
characteristics: Tcell	Enrichment score from GSEA for T cell contamination
Calculation of enrichment scores:	
"To estimate residual sample contamination, we generated separate enrichment scores for neutrophils, B cells, monocytes, T cells, and natural killer cells. We implemented a Gene Set Enrichment Analysis [1] to calculate the enrichment scores using the gene signature of each blood cell type in the ranked list of expression values for each MESA sample. The cell type-specific signature genes were selected from previously defined lists [2] and passed the following additional filters: at least four-fold more highly expressed in the targeted cell type than in other cell populations and low expression levels in the targeted cells. "	
	
Reference List	
[1]	"Subramanian,A. et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc. Natl. Acad. Sci. U. S. A 102, 15545-15550 (2005)."
[2]	"Abbas,A.R. et al. Immune response in silico (IRIS): immune-specific genes identified from a compendium of microarray expression data. Genes Immun. 6, 319-331 (2005)."
