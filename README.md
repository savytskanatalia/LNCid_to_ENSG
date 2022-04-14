# LNCid_to_ENSG
Dictionary for translating LNC and other trivial gene IDs to Ensembl ENSG ID (H.sapiens)




R-script for getting Ensemble ID for genes, whose first annotation ID is different from Ensembl (frequent for lnRNAs, however some genes had such annotation too). In the "annotation" folder original annotation .gtf file and the resulting dictionary.


Info on .gtf:
##description: evidence-based annotation of the human genome (GRCh38), version 29 (Ensembl 94) with added genes from LNCipedia, version 5.2

##provider: GENCODE, LNCipedia

##contact: ehutchins@tgen.org

##format: gff-version 2

##source-version rtracklayer 1.38.3

##date 2018-11-14

Only "gene" (column 3) rows were retained for analysis. 
