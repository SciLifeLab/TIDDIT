FROM python:3.8-slim

ARG TIDDIT_VERSION=2.12.0

LABEL base_image="python:3.8-slim"
LABEL software="TIDDIT.py"
LABEL software.version=${TIDDIT_VERSION}

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
    build-essential \
    cmake \
    make \
    unzip \
    wget \
    zlib1g-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /app

RUN wget https://github.com/SciLifeLab/TIDDIT/archive/TIDDIT-${TIDDIT_VERSION}.zip && \
    unzip TIDDIT-${TIDDIT_VERSION}.zip && \
    rm TIDDIT-${TIDDIT_VERSION}.zip

RUN cd TIDDIT-TIDDIT-${TIDDIT_VERSION} && \
    ./INSTALL.sh && \
    chmod +x /app/TIDDIT-TIDDIT-${TIDDIT_VERSION}/TIDDIT.py && \
    ln -s /app/TIDDIT-TIDDIT-${TIDDIT_VERSION}/TIDDIT.py /usr/local/bin 

ENTRYPOINT ["TIDDIT.py"]
CMD ["--help"]