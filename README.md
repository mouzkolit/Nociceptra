
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

<h2> Guidelines </h2>

<p> Please add the right file path to the analysis scripts </p>
<p> Or just copy the files into the same folder as the analysis scripts </p>
<ul> To use the scripts <ul>
  <li> Generate first the WGCNA affiliation by using the WGCNA_mRNA and WGCNA_miRNA script </li>
  <li> Use the filtered AD2,AD3 and 840 DESeq2 normalized count matrices as input </li>
  <li> If you like to change the filtering process, the python script Cluster_analysis_mrna_sequencing.ipynb provides a mask function </li>
  <li> For the analysis the Cluster_analysis_mrna_sequencing.ipynb should be run first, althoug we provide the miRNA_edge database as downloadable resource using our docker tool </li>
