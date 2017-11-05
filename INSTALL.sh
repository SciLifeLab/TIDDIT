mkdir build
cd build
cmake ..
make 
cd ..
pip install numpy scipy cython
python setup.py build_ext --inplace
