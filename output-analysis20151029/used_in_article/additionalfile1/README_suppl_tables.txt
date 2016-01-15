Additional file 1 contains the following files: 

1. Suppl Chromosomal location: Chromosomal locations with many significant miRNAs
2. boxplots_table 7: Boxplots of the expression of miRNAs differentially expressed in both malignant transitions (miRNAs from Table 7)
3. sampleannotation_suppl: Sample names and grouping, a more complete annotation is found in the github repository; https://github.com/ous-uio-bioinfo-core/Eurocan_microRNA
4. Comparison to Tahiri (Table and Figure): Comparing diffrentially expressed miRNAs in this study and Tahiri et al, 2014 who compared normal - benign and benign - invasive. 
5. consensuslists: Differentially expressen miRNAs with data from both statistical approaches
6. difflists_iClust_UCAM_only: 

Explanation for the genelists (consensuslists and difflists_iClust_UCAM_only):
Differentially expressed miRNAs between the two groups given in the file name.
The adjusted P-values (adj.P.Val) listed are False Discovery Rates (FDR).
Curated miRNAS: White= confirmed miRNAs, purple=possible miRNAs, yellow= equivocal, rejected= do not represent real miRNAs. (Fromm et al, 2015). Be aware that the "rejected" are still included in these lists, but not included in tables where significant miRNAs are counted. They may still be differentially expressed, although they are not microRNAs.
Fold change is in log2, a positive value means that the first group (in the filename) has the highest mean expression, i.e is up-regulated compared to the second group.	

Further explanation for consensus lists:
Log2 fold change is given for the merged approach and the test statistic for the meta approach.	
miRNAs found significantly(FDR<0.05) differentially expressed in the same direction using both approaches are listed as consensus=TRUE.
If the miRNA is altered in the same direction in both analyses: same_dir=TRUE.


