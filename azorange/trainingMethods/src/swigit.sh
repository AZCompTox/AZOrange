#!/bin/tcsh

echo "============================= Compile randomForest.cc ======================================"
g++ -c -fPIC `pkg-config --cflags opencv` -o randomForest.o randomForest.cc

echo "============================== swig -c++ -python randomForest.i ======================================="
swig -c++ -python randomForest.i

echo "============================= Compile randomForest_wrap.cxx ======================================"
g++ -c -fPIC randomForest_wrap.cxx `pkg-config --cflags opencv`

echo "============================= Creating shared object _randomForest.so ================================="
g++ -shared `pkg-config --libs opencv` randomForest.o randomForest_wrap.o -o _randomForest.so

echo "============================= Compile PLS.cc and PlsAPI.cc ======================================"
pymake -nopython PLS.cc
pymake -nopython PlsAPI.cc

echo "============================== swig -c++ -python PlsAPI.i ======================================="
swig -c++ -python PlsAPI.i

echo "======================== pymake -so -nopython PlsAPI_wrap.cxx ==================================="
pymake -so -nopython PlsAPI_wrap.cxx

#echo "===================== cp OBJS/linux.../libPlsAPI_wrap.so ./_PlsAPI.so ==========================="
#cp OBJS/linux-i386___g++_dbg_double_nopython_blas_logging=dbg_numpy/libPlsAPI_wrap.so ./_pls.so


