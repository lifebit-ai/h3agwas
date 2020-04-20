FROM r-base:3.6.0

RUN apt-get update && apt-get install -y \
    libssl-dev \
    libcurl4-openssl-dev

# Install R dependencies
RUN Rscript -e "install.packages('qqman')" && \
    Rscript -e "install.packages('manhattanly')" && \
    Rscript -e "install.packages('plotly')" && \
    Rscript -e "install.packages('tidyverse')" && \
    Rscript -e "install.packages('CMplot')"

# Install plink
RUN mkdir /plink && \
    cd /plink && \
    wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190617.zip && \
    unzip plink_linux_x86_64_20190617.zip && \
    mv plink /usr/local/bin/

RUN apt-get install -y procps && \
    wget https://github.com/jgm/pandoc/releases/download/2.7.3/pandoc-2.7.3-1-amd64.deb && \
    dpkg -i pandoc-2.7.3-1-amd64.deb