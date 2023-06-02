# DavalosHeydari_scANNA
Scripts for scRNAseq preprocessing, data analysis, and figure generation for Davalos &amp; Heydari et al. (2023)

[![DOI:10.1101/2023.05.29.542760](http://img.shields.io/badge/DOI-10.1101/2023.05.29.542760-B31B1B.svg)](https://doi.org/10.1101/2023.05.29.542760)


Package repository can be found at https://github.com/SindiLab/scANNA

[scANNA tutorials](https://github.com/SindiLab/Tutorials/tree/main/scANNA%20single-cell%20ANalysis%20using%20Neural%20Attention) can be found at https://github.com/SindiLab/Tutorials


Repo structure is as follows 

```
.
├── Figure_Scripts
├── GlobalMarkerSelection/
│   └── Triku
├── Preprocessing/
│   └── TrainTestSplitting
└── SupervisedClassification/
    ├── CHETAH
    ├── DataSplitting
    ├── RandomForest
    ├── SingleCellNet
    ├── scClassify
    └── scPred/
        ├── ModelGeneration
        └── Evaluation
```


## Citation
If you found our work useful for your ressearch, please cite our preprint:
```
@article {Davalos2023.05.29.542760,
	author = {Oscar A. Davalos and A. Ali Heydari and Elana J Fertig and Suzanne Sindi and Katrina K Hoyer},
	title = {Boosting Single-Cell RNA Sequencing Analysis with Simple Neural Attention},
	elocation-id = {2023.05.29.542760},
	year = {2023},
	doi = {10.1101/2023.05.29.542760},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2023/06/01/2023.05.29.542760},
	eprint = {https://www.biorxiv.org/content/early/2023/06/01/2023.05.29.542760.full.pdf},
	journal = {bioRxiv}
}
```
