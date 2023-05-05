# Fig.3: Enhancers/promoters regulating transcriptional programs in cancer.

---

**Directory structure**

* `make_mutliome_objects` -- First step of making multiome cohort-level objects.

* `ACR-to-gene_links` -- Second step of computing ACR-to-gene links using the multiome objects and plotting results
  * `Compute_links_to_genes_v7.0.R` - calculates ACR-to-gene links per cohort object
  * `Compute_diffuse_links_to_genes_v7.0.R` - calculates 100kb genomic windows -to-gene links per cohort (aka diffused links)
  * `1_Annotate_links_by_enhancers_geneHancerInter_filter_diffuse_DACRs.ipynb` - annotate the resulting links by EpiMap, ENCODE and GeneHancer databases and filter links by diffuse links, add inferCNV annotation
  * `2_linked_genes_summary.ipynb` - plot summary plots of links: 3b, 3c, S6a
  * `3_Heatmaps_links_coveragePlots_primary_vs_normal.ipynb` - plot heatmaps with links 3d and S6e, plot coverage plots in 3e, S6b, S6d, S6f, S6g
