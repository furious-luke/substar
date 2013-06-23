#ifndef bolshoi_types_hh
#define bolshoi_types_hh

#include <libhpc/libhpc.hh>

///
/// The Bolshoi halo data structure.
///
/// I've commented out a lot of values to reduce the size.
/// Unfortunately there is one enormous forest that consumes
/// an awful lot of memory.
///
struct bolshoi_halo_type
{
   // float     scale; // Scale factor of halo.
   long long  id; // ID of halo (unique across entire simulation).
   // float     desc_scale; // Scale of descendant halo, if applicable.
   long long  desc_id; // ID of descendant halo, if applicable.
   unsigned   num_prog; // Number of progenitors.
   long long  pid; // Host halo ID.
   // long long  upid; // Most massive host halo ID (only different from pid in
   //                  // cases of sub-subs or sub-sub-subs etc.
   // long long  desc_pid; // pid of descendant halo, if applicable.
   // int        phantom; // Nonzero for halos interpolated across timesteps.
   // float     sam_mvir; // Halo mass, smoothed across accretion history; always greater
   //                      // than sum of halo masses of contributing progenitors (Msun/h).
   //                      // Only for use with select semi-analytic models.
   float     mvir; // Halo mass (Msun/h).
   float     rvir; // Halo radius (kpc/h comoving).
   // float     rs; // Scale radius (kpc/h comoving).
   float     vrms; // Velocity dispersion (km/s physical).
   // int        mmp; // Whether the halo is the most massive progenitor or not.
   // float     scale_of_last_mm; // Scale factor of the last major merger (mass ratio > 0.3).
   float     vmax; // Maximum circular velocity (km/s physical).
   float     x, y, z; // Halo position (Mpc/h comoving).
   float     vx, vy, vz; // Halo velocity (km/s physical).
   float     jx, jy, jz; // Halo angular momenta ((Msub/h)*(Mpc/h)*km/s physical).
   float     spin; // Halo spin parameter.
   // unsigned   breadth_first_id; // Breadth-first ordering of halos within a tree.
   // unsigned   depth_first_id; // Depth-first ordering of halos withing a tree.
   // long long  tree_root_id; // ID of the halo at the last timestep in the tree.
   // long long  orig_halo_id; // Original halo ID from the halo finder.
   unsigned   snap_num; // Snapshot number from which halo originated.
   // unsigned   next_coprog_depthfirst_id; // Depth-first ID of next coprogenitor.
   // unsigned   last_coprog_depthfirst_id; // Depth-first ID of last progenitor.
   // float     rs_klypin; // Scale radius determined using Vmax and Mvir (see Rockstar paper).
   // float     mvir_all; // Mass enclosed within the specified overdensity, including unbound
   //                      // particles (Msun/h).
   // float     m200b_m2500c[4]; // Mass enclosed within specified overdensities (Msun/h).
   // float     x_offs; // Offset of density peak from average particle position (kpc/h comoving).
   // float     v_offs; // Offset of density peak from average particle velocity (km/s physical).
   // float     spin_bullock; // Bullock spin parameter (J/(sqrt(2)*GMVR)).
   // float     b_to_a, c_to_a; // Ration of second and third largest shape ellipsoid axes (B and C)
   //                            // to largest shape ellipsoid axis (A) (dimensionless).
   // float     ax, ay, az; // Largest shape ellipsoid axis (kpc/h comoving).
};

///
/// The SAGE compatible halo finder structure.
///
struct sage_halo_type
{
   // merger tree pointers 
   int descendant;
   int first_progenitor;
   int next_progenitor;
   int first_halo_in_fof_group;
   int next_halo_in_fof_group;

   // properties of halo 
   int Len;
   float M_Mean200, Mvir, M_TopHat;  // Mean 200 values (Mvir=M_Crit200)
   float Pos[3];
   float Vel[3];
   float VelDisp;
   float Vmax;
   float Spin[3];
   long long MostBoundID;

   // original position in subfind output 
   int SnapNum;
   int FileNr;
   int SubhaloIndex;
   float SubHalfMass;
};

#endif
