
<h1 style="color:blue;"> NOCICEPTRA Resource </h1>

<h2> Welcome to Git repository for Data Analysis </h2>

<ul>
<li> Vst and normalized count data generated using the R-script </li>
<li> WGCNA R script to generate Module assignments and Hub-genes </li>
<li> Python scripts </li>
<ul> 
  <li> Script for Supercluster detection </li>
  <li> miRNA Edge Database construction </li>
  <li> Hub-gene networks and miRNA:hub-gene target-spaces </li>
  <li> miRNA enrichment analysis and mRNA enrichment analysis </li>
</ul>
<li> <a href = "https://hub.docker.com/repository/docker/muiphysiologie/nociceptra_mui"> NOCICEPTRA Tool at Docker </li>
</ul>

<h2> Pipelines </h2>

<p> RNA and small RNA sequencing were aligned to the human reference genome GRCH.38.p13 (genecode) with the spice-aware STAR-aligner and the parameter as indicated in the Publication
and normalized using DESseq2. TPMs were generated with a second pipeline using CLC-workbench with the GRCh38.p7 reference genome from the NCBI webpage</p>
<p> For both small and long RNA sequencing only unqiuely matched reads were counted. Since small RNA sequencing has also a lot of multipmappers we have employed a second strategy counting multimappers using the CLC Genomic Workbench and annotated using the miRBase 21 release (human and mouse) as well as the Homo_sapiens.GRCh38.ncrna annotation file. This data are provided within the Github repository but were not used in the publication, since we focused only on miRNAs.

<h2> Guidelines </h2>

<p> Please add the right file path to the analysis scripts </p>
<p> Or just copy the files into the same folder as the analysis scripts </p>
<ul> To use the scripts: <ul>
  <li> Generate first the WGCNA affiliation by using the WGCNA_mRNA and WGCNA_miRNA script to obtain module affiliation of genes and miRNAs </li>
  <li> Use the filtered AD2,AD3 and 840 DESeq2 normalized count matrices as input or the whole matrices </li>
  <li> If you like to change the filtering process, the python script Cluster_analysis_mrna_sequencing.ipynb provides a mask function for read counts </li>
  <li> For the analysis the Cluster_analysis_mrna_sequencing.ipynb should be run first, althoug we provide the miRNA_edge database as downloadable resource which can be found in the our Docker Tools </li>
  </ul> </ul>
  
  
 <h2> Data section </h2>
 <p> In the Data section we included data necessary for the analysis, but not databases from miRTarbase StarBase and TarBase </p>
 <p> We added the Ion channel Table --> Ion_channels.xlsx </p>
 <p> Raw counts for miRNAs and mRNAs derived from htseq, vst counts, as well as filtered counts which are used for the WGCNA analysis </p>
 <p> We further added gene enrichment tables, gene tables obtained from Pain database and a Day2Day Variance analysis where you can check differentially expressed genes between celllines at different days </p>
  
