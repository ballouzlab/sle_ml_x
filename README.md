##  X chromosome transcriptional signatures in systemic lupus erythematosus: a machine learning analysis

This is the repository for the work:  
*_Machine learning reveals X chromosome transcriptional signatures that predict systemic lupus erythematosus in females_* 

Lachlan G. Gray <sup>1,2</sup>, Aiden Telfser  <sup>1,2</sup>, Hannah Hu <sup>1,2,3</sup>, Alisa Kane <sup>3</sup>, Elissa K. Deenick  <sup>1,4</sup>, Tri Giang Phan <sup>1,2</sup> and Sara Ballouz <sup>5</sup>
1.	Precision Immunology Program, Garvan Institute of Medical Research, Sydney NSW, Australia
2.	St Vincent’s Healthcare Clinical Campus, School of Clinical Medicine, Faculty of Medicine and Health, UNSW Sydney, Sydney, NSW, Australia
3.	Department of Immunology, Allergy and HIV, St Vincent’s Hospital, Sydney, NSW, Australia
4.	Kirby Institute, Faculty of Medicine and Health, UNSW Sydney, Sydney, NSW, Australia
5.	School of Computer Science and Engineering, Faculty of Engineering, UNSW Sydney, Sydney, NSW, Australia
 

### Abstract
#### Background
Systemic lupus erythematosus (SLE) has a significant female bias; however, it remains unclear whether X-chromosomal dysregulation increases the risk in females. Furthermore, while the role of X-linked genes in SLE pathogenesis is recognised, the analysis of genetic risks in SLE has largely focussed on autosomal genes.
#### Methods
Here, we analysed a public single-cell RNA sequencing (scRNA-seq) dataset of peripheral blood mononuclear cells using machine learning to test whether X-linked gene expression predicts female SLE in a cell type-specific manner. Using pseudobulked expression profiles from 1.1 million cells across 228 female donors, we trained an ensemble of classifiers on all X-chromosomal genes or a literature-derived set of SLE-associated genes. 
#### Results
In CD4⁺ and CD8⁺ T cells, models based on X-chromosomal genes achieved classification performance comparable to SLE gene models and were significantly enriched for XCI escape genes, with IL2RG, CD40LG, and TKTL1 emerging as key predictors. When applied to an independent paediatric–adult SLE cohort, CD8⁺ T-cell models retained robust performance, indicating a stable and reproducible X-linked signature in this subset. Flow cytometry revealed a modest increase in IL2RG protein expression in SLE T cells, consistent with partial escape from XCI in these cells. 
#### Conclusions
These findings demonstrate that X-linked transcriptional programs, enriched for known XCI escape genes and concentrated in activated T cells, encode a reproducible transcriptional signal of female SLE and provide a framework for dissecting the contributions of sex chromosomes to autoimmunity.
