eval "$(conda shell.bash hook)"
shopt -s extglob                                              # In bash, to use rm -- !(file.txt), you must enable extglob:
rm !("Transport_Diffusion.txt"|"Transport_mu_3tc.txt"|"inf_to_zero_3rho.py"|"Transport_mu_1tc.txt"|"Transport.py"|"TvNR.py"|"chemistry.py"|"Mixlay_0D_000.py"|"export_IC.py"|"RPC_AM.sh"|"RIRIDIS5_AM.sh"|"Iridis_1gpu.slurm") -fr     # Remove everything but the initial files
conda activate OpenSBLI
python Mixlay_0D_000.py

python export_IC.py
python TvNR.py
python Transport.py
python inf_to_zero_3rho.py

python $OPS_TRANSLATOR/ops.py opensbli.cpp
cp /home/ali/app/opensbli/apps/Makefile ./
make clean
make opensbli_mpi
mpirun -np 4 ./opensbli_mpi
