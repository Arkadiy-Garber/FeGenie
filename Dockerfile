FROM ncbi/blast

# https://stackoverflow.com/questions/44331836/apt-get-install-tzdata-noninteractive
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -y -m update && apt-get install -y --no-install-recommends hmmer prodigal python3
RUN apt-get install -y r-base
RUN which make
RUN which g++

RUN Rscript -e 'install.packages(c("remotes"), repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'remotes::install_version("cowplot", "0.9.2", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages(c("ggplot2", "reshape", "argparse", "ggpubr"), repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages(c("reshape2", "ggdendro", "grid", "pvclust"), repos = "http://cran.us.r-project.org")'

RUN apt-get autoclean && \
	apt-get -y autoremove libssl-dev && \
	rm -rf /var/lib/apt/lists/*

WORKDIR /data

CMD ["/bin/bash"]
