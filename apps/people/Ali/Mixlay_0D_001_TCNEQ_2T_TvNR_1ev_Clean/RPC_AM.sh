eval "$(conda shell.bash hook)"
shopt -s extglob                                              # In bash, to use rm -- !(file.txt), you must enable extglob:
rm !("TvNR.py"|"chemistry.py"|"Mixlay_2D_010.py"|"export_IC.py"|"RPC_AM.sh"|"RIRIDIS5_AM.sh"|"Iridis_1gpu.slurm") -fr     # Remove everything but the initial files
conda activate OpenSBLI
python Mixlay_2D_010.py

python export_IC.py
python TvNR.py

python $OPS_TRANSLATOR/ops.py opensbli.cpp
cp /home/ali/app/opensbli/apps/Makefile ./
make clean
make opensbli_mpi
mpirun -np 4 ./opensbli_mpi
