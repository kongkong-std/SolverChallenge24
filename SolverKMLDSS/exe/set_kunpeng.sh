rm -rf ../src/kml_*
module purge
module use /workspace/public/software/modules
module load compilers/gcc/kunpenggcc/10.3.1/gcc10.3.1
module load kml-dss/kmldss/kmldss-gcc
module load kml-dss/openblas/openblas-gcc
export OMP_PROC_BIND=close

make all