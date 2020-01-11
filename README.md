# DiS-TSS
Diana Signal - TSS

## Abstract
Promoters are genomic loci that comprise of all the necessary DNA motifs responsible for transcription initiation, thus being critical gene expression regulators. The spread, distribution and utilization of transcription start sites (TSS) within promoters  are poorly understood. Cap Analysis of Gene Expression (CAGE) has emerged as a popular gene expression profiling protocol, able to quantitate TSS usage by recognizing the 5’ end of capped RNA molecules. However, there is an increasing volume of studies in the literature suggesting that CAGE can also detect 5’ capping events which are transcription byproducts. These findings highlight the need for computational methods that can effectively remove the excessive amount of noise from CAGE samples, leading to accurate TSS annotation and promoter usage quantification. In this study, we present an annotation agnostic computational framework, DIANA Signal-TSS (DiS-TSS), that exclusively utilizes CAGE signal features in the spatial and frequency domains to accurately distinguish between peaks related to real transcription initiation events and biological or protocol-induced noise. 
<br>
Examle usage of tool:<br>
python DiStss.py -i <input_file> [OPTIONS]<br>
<br>
-i&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp- CAGE-seq, bam file, (provide full path e.g. /home/user/data.bam)<br>
-clusters&nbsp&nbsp- Bed file containing clusters. If set to "0" the tool will auto estimate the clusters. Default is 0 (auto).<br>
-@\t\t\t\t&nbsp- Number of theads to use. Default is 1.<br>
-MAPQ      - Skip alignments with map quality smaller than <-MAPQ>. Default is 10.<br>
-MAPDist   - Minimum distance for CTSS to be merged. Default is 25.<br>
-tpm       - Clusters expression below <-tpm> are expelled. Default is 1.<br>
<br>
