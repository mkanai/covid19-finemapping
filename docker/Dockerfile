FROM r-base:4.0.2

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Install additional components
RUN apt-get update --fix-missing && \
    apt-get install -y \
        bedtools \
        gawk \
        libcurl4-openssl-dev \
        libopenblas-dev \
        libssl-dev \
        libxml2-dev \
        python3-dev \
        python3-pip \
        && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# HTSLIB (for bgzip)
ENV HTSLIB_VERSION 1.10.2
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 -P ~ && \
    tar xjvf ~/htslib-${HTSLIB_VERSION}.tar.bz2 -C ~ && \
    cd ~/htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/usr --disable-bz2 --disable-lzma && \
    make && make install && \
    rm -r ~/htslib-${HTSLIB_VERSION}*

# SUSIE
# dev version: 0.9.1.0.ef8e7ad
ENV SUSIE_VERSION ef8e7ad
RUN Rscript -e 'install.packages(c("devtools", "argparse", "dplyr", "stringr", "data.table", "here", "R.utils", "Matrix"), noCache = TRUE, repos="https://cloud.r-project.org")' && \
    Rscript -e "Sys.setenv(TAR = \"/bin/tar\"); devtools::install_github(\"stephenslab/susieR@${SUSIE_VERSION}\")"

# for pysam
ENV HTSLIB_CONFIGURE_OPTIONS "--disable-bz2 --disable-lzma"
# pip
RUN pip install --no-cache cython && \
    pip install --no-cache numpy scipy pandas pybedtools

# Add scripts
COPY \
    python/make_finemap_inputs.py \
    R/run_abf.R \
    R/run_susieR.R \
    /usr/local/bin/
RUN chmod +x /usr/local/bin/*.py /usr/local/bin/*.R

VOLUME [ "/root/.config/gcloud" ]
CMD [ "/bin/bash" ]
