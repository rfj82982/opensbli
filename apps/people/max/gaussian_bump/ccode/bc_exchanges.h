// Boundary condition exchange code on opensbliblock00 direction 2 left
ops_halo_group exchange75_block0 ;
{
int halo_iter[] = {block0np0 + 7, block0np1 + 7, 3};
int from_base[] = {-3, -3, 0};
int to_base[] = {-3, -3, block0np2};
int from_dir[] = {1, 2, 3};
int to_dir[] = {1, 2, 3};
ops_halo halo0 = ops_decl_halo(D00_B0, D00_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo halo1 = ops_decl_halo(D01_B0, D01_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo halo2 = ops_decl_halo(D10_B0, D10_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo halo3 = ops_decl_halo(D11_B0, D11_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo halo4 = ops_decl_halo(detJ_B0, detJ_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo grp[] = {halo0,halo1,halo2,halo3,halo4};
exchange75_block0 = ops_decl_halo_group(5,grp);
}
// Boundary condition exchange code on opensbliblock00 direction 2 right
ops_halo_group exchange76_block0 ;
{
int halo_iter[] = {block0np0 + 7, block0np1 + 7, 4};
int from_base[] = {-3, -3, block0np2 - 3};
int to_base[] = {-3, -3, -3};
int from_dir[] = {1, 2, 3};
int to_dir[] = {1, 2, 3};
ops_halo halo0 = ops_decl_halo(D00_B0, D00_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo halo1 = ops_decl_halo(D01_B0, D01_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo halo2 = ops_decl_halo(D10_B0, D10_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo halo3 = ops_decl_halo(D11_B0, D11_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo halo4 = ops_decl_halo(detJ_B0, detJ_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo grp[] = {halo0,halo1,halo2,halo3,halo4};
exchange76_block0 = ops_decl_halo_group(5,grp);
}
// Boundary condition exchange code on opensbliblock00 direction 2 left
ops_halo_group exchange65_block0 ;
{
int halo_iter[] = {block0np0 + 7, block0np1 + 7, 3};
int from_base[] = {-3, -3, 0};
int to_base[] = {-3, -3, block0np2};
int from_dir[] = {1, 2, 3};
int to_dir[] = {1, 2, 3};
ops_halo halo0 = ops_decl_halo(rho_B0, rho_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo halo1 = ops_decl_halo(rhou0_B0, rhou0_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo halo2 = ops_decl_halo(rhou1_B0, rhou1_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo halo3 = ops_decl_halo(rhou2_B0, rhou2_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo halo4 = ops_decl_halo(rhoE_B0, rhoE_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo grp[] = {halo0,halo1,halo2,halo3,halo4};
exchange65_block0 = ops_decl_halo_group(5,grp);
}
// Boundary condition exchange code on opensbliblock00 direction 2 right
ops_halo_group exchange66_block0 ;
{
int halo_iter[] = {block0np0 + 7, block0np1 + 7, 4};
int from_base[] = {-3, -3, block0np2 - 3};
int to_base[] = {-3, -3, -3};
int from_dir[] = {1, 2, 3};
int to_dir[] = {1, 2, 3};
ops_halo halo0 = ops_decl_halo(rho_B0, rho_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo halo1 = ops_decl_halo(rhou0_B0, rhou0_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo halo2 = ops_decl_halo(rhou1_B0, rhou1_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo halo3 = ops_decl_halo(rhou2_B0, rhou2_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo halo4 = ops_decl_halo(rhoE_B0, rhoE_B0, halo_iter, from_base, to_base, from_dir, to_dir);
ops_halo grp[] = {halo0,halo1,halo2,halo3,halo4};
exchange66_block0 = ops_decl_halo_group(5,grp);
}