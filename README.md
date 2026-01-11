# Annapragada2025

This project is a workflowr repository containing code needed to
generate all figures and analyses in “Cell-free DNA fragmentomes for non-invasive detection of liver cirrhosis and other diseases”.


There are 5 folders of interest in this workflowr.

1)  analysis - contains code needed to generate each figure and
    supplementary figure of the paper.

``` r
library(workflowr)
wflow_build("analysis/*.Rmd")
```

2)  code - [empty - folder existence needed for workflowr structure]

3)  data - This contains raw data used to generate the figures and
    tables

4)  docs - This contains html of the markdown files in analysis, as well
    as the generated figures.

5)  output - See README in this folder for more details. Some required files for journal and convenience functions. Also contains code for model training. To use this code one must first generate ARTEMIS and DELFI features using the following code repositories: 
<https://github.com/cancer-genomics/delfi3>
<https://github.com/cancer-genomics/artemis_pipeline>

SessionInfo.Rmd in the analysis folder, and the corresponding html in
the docs folder document the versions of packages used to compile this
repo.

This repository is available on Github, and may be run as a workflowr
project to generate a webpage with all code and figures linked.
[workflowr](https://github.com/jdblischak/workflowr)

