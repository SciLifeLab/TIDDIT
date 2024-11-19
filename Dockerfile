FROM condaforge/mambaforge:24.9.2-0

WORKDIR /app

## Set TIDDIT version
ARG TIDDIT_VERSION=3.9.0

## Add some info
LABEL base_image="python:3.8-slim"
LABEL software="TIDDIT.py"
LABEL software.version=${TIDDIT_VERSION}

## Download and extract
RUN conda install conda-forge::unzip
RUN conda install -c conda-forge pip gcc joblib
RUN conda install -c bioconda numpy cython pysam bwa

RUN wget https://github.com/SciLifeLab/TIDDIT/archive/TIDDIT-${TIDDIT_VERSION}.zip && \
    unzip TIDDIT-${TIDDIT_VERSION}.zip && \
    rm TIDDIT-${TIDDIT_VERSION}.zip

## Install
RUN cd TIDDIT-TIDDIT-${TIDDIT_VERSION}  && \
	pip install -e .

ENTRYPOINT ["tiddit"]
CMD ["--help"]
