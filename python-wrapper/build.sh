CONDA_ENV_PATH=/opt/anaconda/miniconda3/envs/test_cython/lib/
python setup.py build_ext --inplace --bitpit-path=/home/jc.giret/bitpit/ --madeleine-path=/home/jc.giret/SOURCES/madeleine_parallel/build/src/ --extensions-source=coupling.pyx  --metis-path=/home/jc.giret/.local/lib/ --lapack-path=$CONDA_ENV_PATH --mpi-include-path=$CONDA_ENV_PATH --petsc-path=$CONDA_ENV_PATH
