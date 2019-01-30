mkdir build
cd build
cmake ..
make 
cd ..
pip install numpy scipy cython pysam
cd src
python setup.py build_ext --inplace
