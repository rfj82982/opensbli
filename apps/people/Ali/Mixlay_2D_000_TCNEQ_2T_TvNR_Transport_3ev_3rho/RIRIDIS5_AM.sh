eval "$(conda shell.bash hook)"
shopt -s extglob                                              # In bash, to use rm -- !(file.txt), you must enable extglob:
rm !("TvNR.py"|"Mixlay_0D_000.py"|"export_IC.py"|"RPC_AM.sh"|"RIRIDIS5_AM.sh"|"Iridis_1gpu.slurm") -fr     # Remove everything but the initial files
#find . ! -name "Mixlay_2D_020.py" ! -name "Mixlay_2D_Plot.py" ! -name "Run_AM.sh" -type d -exec rm -f -r {} +
conda activate opensbli
python Mixlay_0D_000.py

#python export_IC.py
#python $OPS_TRANSLATOR/ops.py opensbli.cpp
#cp ~/app/opensbli_v1_github/apps/Makefile ./
#make clean
#make opensbli_cuda

# make opensbli_mpi
# mpirun -np 4 ./opensbli_mpi
