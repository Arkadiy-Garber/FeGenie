#!/bin/bash

tar -xvf test_dataset.tar.gz

Rscript -e 'install.packages("grid", repos = "http://cran.us.r-project.org")'
Rscript -e 'install.packages("broom", repos = "http://cran.us.r-project.org”)'
Rscript -e 'install.packages("ggpubr", repos = "http://cran.us.r-project.org”)’
Rscript -e 'install.packages("ggplot2", repos = "http://cran.us.r-project.org")'
Rscript -e 'install.packages("reshape", repos = "http://cran.us.r-project.org")'
Rscript -e 'install.packages("reshape2", repos = "http://cran.us.r-project.org")'
Rscript -e 'install.packages("tidyverse", repos = "http://cran.us.r-project.org")'
Rscript -e 'install.packages("argparse", repos = "http://cran.us.r-project.org")'
Rscript -e 'install.packages("ggdendro", repos = "http://cran.us.r-project.org")'
Rscript -e 'install.packages("ggpubr", repos = "http://cran.us.r-project.org")'
Rscript -e 'install.packages("grid", repos = "http://cran.us.r-project.org")'
Rscript -e 'install.packages("pvclust", repos = "http://cran.us.r-project.org")'
