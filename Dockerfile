FROM condaforge/mambaforge:24.9.2-0

WORKDIR /app

## Set TIDDIT version
ARG TIDDIT_VERSION=3.9.5

## Add some info
LABEL base_image="python:3.8-slim"
LABEL software="TIDDIT.py"
LABEL software.version=${TIDDIT_VERSION}

## Install dependencies
RUN conda install -c conda-forge pip gcc joblib
RUN conda install -c bioconda numpy cython pysam bwa

## Copy local source
COPY . /app/TIDDIT

## Install
RUN cd /app/TIDDIT && \
	pip install -e .

ENTRYPOINT ["tiddit"]
CMD ["--help"]
