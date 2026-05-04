FROM --platform=linux/amd64 rocker/tidyverse:4.3.1

COPY . /StratGWAS

RUN find /StratGWAS/src -name "*.so" -delete \
 && find /StratGWAS/src -name "*.o" -delete \
 && Rscript -e "install.packages(c('remotes', 'RcppParallel', 'RcppEigen'), repos='https://cloud.r-project.org')" \
 && R CMD INSTALL /StratGWAS