# DiS-TSS (Diana Signal - TSS)

Promoters are genomic loci that comprise of all the necessary DNA motifs responsible for transcription initiation, thus being critical gene expression regulators. The spread, distribution and utilization of transcription start sites (TSS) within promoters  are poorly understood. Cap Analysis of Gene Expression (CAGE) has emerged as a popular gene expression profiling protocol, able to quantitate TSS usage by recognizing the 5’ end of capped RNA molecules. However, there is an increasing volume of studies in the literature suggesting that CAGE can also detect 5’ capping events which are transcription byproducts. These findings highlight the need for computational methods that can effectively remove the excessive amount of noise from CAGE samples, leading to accurate TSS annotation and promoter usage quantification. In this study, we present an annotation agnostic computational framework, DIANA Signal-TSS (DiS-TSS), that exclusively utilizes CAGE signal features in the spatial and frequency domains to accurately distinguish between peaks related to real transcription initiation events and biological or protocol-induced noise. 

### Prerequisites

In order to run the tool you need to have installed python version 3.<br>
Also the following libraries python must be met:
```
kurtosis and skew from scipy.stats
find_peaks from scipy.signal
trapz and array from numpy
pickle
pandas
subprocess
```
### Running the algorithm
```
Examle usage of tool:
python DiStss.py -i <input_file> [OPTIONS]

-i         - CAGE-seq, bam file, (provide full path e.g. /home/user/data.bam)
-clusters  - Bed file containing clusters. If set to "0" the tool will auto estimate the clusters. Default is 0 (auto).
-@         - Number of theads to use. Default is 1.
-MAPQ      - Skip alignments with map quality smaller than <-MAPQ>. Default is 10.
-MAPDist   - Minimum distance for CTSS to be merged. Default is 25.
-tpm       - Clusters expression below <-tpm> are expelled. Default is 1.
```
### Example run
Make sure you are running python version 3.

Basic run letting the algorithm calculate the clusters from the bam file. Always provide full path to bam file.
```
python Dis_main.py -i /home/user/cage.bam
```
Provide your own custom clusters file. Clusters file must be in BED format.
```
python Dis_main.py -i /home/user/cage.bam -clusters home/clusters.bed
```


### Authors
Dimitris Grigoriadis, Nikos Perdikopanis, Georgios K Georgakilas and Artemis Hatzigeorgiou.

### Please cite
TBC...

### License
Free License.

### Acknowledgments
##### Contributions
A.H. supervised the study. N.P. and D.G. envisioned and designed the algorithm and performed the unsupervised analysis. D.G. developed the algorithm. N.P performed all comparisons. A.H., D.G. and N.P. made the figures and wrote the manuscript.
