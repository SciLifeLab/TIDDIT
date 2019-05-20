BootStrap: debootstrap
OSVersion: trusty
MirrorURL: http://us.archive.ubuntu.com/ubuntu/


%runscript
    echo "This is what happens when you run the container..."


%post
    echo "Hello from inside the container"
    sed -i 's/$/ universe/' /etc/apt/sources.list
    apt-get update
    apt-get upgrade
    apt-get -y --force-yes install build-essential cmake make zlib1g-dev python python-dev python-setuptools git
    easy_install pip
    pip install numpy cython

    git clone https://github.com/SciLifeLab/TIDDIT.git
    mv TIDDIT/* /bin/
    cd /bin/ && ./INSTALL.sh
    chmod +x /bin/TIDDIT.py


