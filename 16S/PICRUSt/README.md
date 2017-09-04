# Running PICRUSt on the UCT Hex cluster.

##This is based on using the OTU table produced by the standard cbio pipeline for de novo OTU picking
1. Run step one of uparse_closed_ref_from_de_novo_forPICRUSt.sh with the file 'otus_repsetOUT.fa' from the standard cbio pipeline as input
2. Download OTU table created with standard cbio-pipeline, import into R using phyloseq and save as a phyloseq object (.RData file)
3. You will need the .RData file from 2. for the R script WISH_de_novo_to_closed_ref_OTUs.R together with the .uc file created in 1. 
4. Upload filtere .txt file containing only OTUs that mapped to GreenGenes IDs to hex
5. Continue with remaining steps in uparse_closed_ref_from_de_novo_forPICRUSt.sh (create new .biom file for use in PICRUSt)
6. Run PICRUSt on output from 5.
