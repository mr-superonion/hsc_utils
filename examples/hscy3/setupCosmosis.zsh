source cosmosis-configure

# mpi thread
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1

#cuda
export CUDA_VISIBLE_DEVICES=-1
export TF_CPP_MIN_LOG_LEVEL=2
export TF_ENABLE_ONEDNN_OPTS=0

export CC=gcc
export CXX=g++
export FC=gfortran
export MPIFC=mpif90

# Your cosmosis standard library
export cosmosis_dir=/work/xiangchong.li/cosmosis2-mpich/cosmosis-standard-library

# The additional cosmosis files
export cosmosis_utils=$PWD/cosmosis_utils
