
<h1 style="color:blue;"> NOCICEPTRA Resource </h1>

<h2> Welcome to Git repository for Data Analysis </h2>

<ul>
<li> Vst and normalized count data generated using the DESeq2 </li>
<li> TPMs generate with a second pipeline using CLC Workbench </li>
<li> WGCNA to generate conserved Module assignments and Hub-genes </li>
<li> Python scripts </li>
<ul> 
  <li> Script for Supercluster detection </li>
  <li> miRNA Edge Database construction </li>
  <li> Hub-gene networks and miRNA:hub-gene target-spaces </li>
  <li> miRNA enrichment analysis and mRNA disease affiliation analysis </li>
  <li> Gene enrichment and disease analysis </li>
  
</ul>
<li> <a href = "https://hub.docker.com/r/muiphysiologie/nociceptra_mui"> NOCICEPTRA Tool at Docker </li>
</ul>

<h2> Pipelines </h2>

<p> RNA and small RNA sequencing were aligned to the human reference genome GRCH.38.p13 (genecode) with the spice-aware STAR-aligner and the parameter as indicated in the Publication
and normalized using DESseq2. TPMs were generated with a second pipeline using CLC-workbench with the GRCh38.p7 reference genome from the NCBI webpage</p>
<p> For both small and long RNA sequencing only unqiuely matched reads were counted.Multimappers for miRNAs and ncRNAs were added in the files multimapper section and in the NOCICEPTRA Tool </p>



* Analysis performed by mouzkolit (GitHub Account)
