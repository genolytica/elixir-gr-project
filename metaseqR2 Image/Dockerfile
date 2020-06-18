FROM bioconductor/bioconductor_docker:RELEASE_3_11
# https://hub.docker.com/r/bioconductor/bioconductor_docker/tags

RUN apt-get update && apt-get install -y ca-certificates pandoc

# Install R packages
RUN R -e "install.packages(c('getopt','optparse','devtools','remotes'), repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('pmoulos/metaseqR2-local')"