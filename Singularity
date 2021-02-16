BootStrap: debootstrap
OSVersion: trusty
MirrorURL: http://us.archive.ubuntu.com/ubuntu/

%environment
    SHELL=/bin/bash
    PATH=/opt/anaconda/bin:${PATH}


%runscript
    echo "This is what happens when you run the container..."
    PATH=/opt/anaconda/bin:${PATH}
    echo "This is what happens when you run the container..."

%post
    echo "Hello from inside the container"
    sed -i 's/$/ universe/' /etc/apt/sources.list
    apt-get update
    apt-get upgrade
    apt-get -y --force-yes install build-essential cmake make zlib1g-dev python python-dev python-setuptools git

    cd /root/ && wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    cd /root/ && chmod 700 ./Miniconda3-latest-Linux-x86_64.sh
    cd /root/ && bash ./Miniconda3-latest-Linux-x86_64.sh -b -p /opt/anaconda/ 

    export PATH=/opt/anaconda/bin:${PATH} 


    pip install numpy cython

    git clone https://github.com/SciLifeLab/TIDDIT.git
    mv TIDDIT/* /bin/
    cd /bin/ && ./INSTALL.sh
    chmod +x /bin/TIDDIT.py


