
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
<p> For both small and long RNA sequencing only unqiuely matched reads were counted.Multimappers for miRNAs and ncRNAs were added in the files multimapper Section and in the NOCICEPTRA Tool </p>

<h2>Pipelines Guidelines </h2>

<p> First RNA and small RNA sequencing were aligned with STAR and counted with HTseq, then DESeq2 was used to normalized the data for all three cell lines
  Data was then subsequently parsed to the WGCNA script to detect coexpressed genes and miRNAs and derive modules.
  The downstream analysis was performed GO-profiling, KEGG pathways, miRNA::mRNA interactions and more!
</p>
  
  
 <h2> Data section </h2>
 <p> In the Data section we included data necessary for the analysis, but not databases from miRTarbase StarBase and TarBase </p>
 <p> We added the Ion channel Table --> Ion_channels.xlsx </p>
 <p> Raw counts for miRNAs and mRNAs derived from htseq, vst counts, as well as filtered counts which are used for the WGCNA analysis </p>
 <p> We further added gene enrichment tables, gene tables obtained from Pain database and a Day2Day Variance analysis where you can check differentially expressed genes between celllines at different days </p>
 
 <p> Furthermore we added a less stringent miRNA analysis using the CLC workbench and adjust for multimappers. This data can be found in the files section in the nc_counts_with_multimappers folder, but only uniquely mapped reads for miRNAs were used for the manuscript </p>
  
