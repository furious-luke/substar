#include "types.hh"

using namespace hpc;

///
///
///
void
make_hdf5_types( h5::datatype& mem_type,
		 h5::datatype& file_type )
{
   // Create memory type.
   mem_type.compound( sizeof(bolshoi_halo_type) );
   mem_type.insert( long long, "globally unique halo ID", HOFFSET( bolshoi_halo_type, id ) );
   mem_type.insert( long long, "descendant ID", HOFFSET( bolshoi_halo_type, desc_id ) );
   mem_type.insert( unsigned, "number of progenitors", HOFFSET( bolshoi_halo_type, num_prog ) );
   mem_type.insert( long long, "host halod ID", HOFFSET( bolshoi_halo_type, pid ) );
   mem_type.insert( float, "virial mass", HOFFSET( bolshoi_halo_type, mvir ) );
   mem_type.insert( float, "virial radius", HOFFSET( bolshoi_halo_type, rvir ) );
   mem_type.insert( float, "vrms", HOFFSET( bolshoi_halo_type, vrms ) );
   mem_type.insert( float, "maximum velocity", HOFFSET( bolshoi_halo_type, vmax ) );
   mem_type.insert( float, "x position", HOFFSET( bolshoi_halo_type, x ) );
   mem_type.insert( float, "y position", HOFFSET( bolshoi_halo_type, y ) );
   mem_type.insert( float, "z position", HOFFSET( bolshoi_halo_type, z ) );
   mem_type.insert( float, "x velocity", HOFFSET( bolshoi_halo_type, vx ) );
   mem_type.insert( float, "y velocity", HOFFSET( bolshoi_halo_type, vy ) );
   mem_type.insert( float, "z velocity", HOFFSET( bolshoi_halo_type, vz ) );
   mem_type.insert( float, "x spin", HOFFSET( bolshoi_halo_type, jx ) );
   mem_type.insert( float, "y spin", HOFFSET( bolshoi_halo_type, jy ) );
   mem_type.insert( float, "z spin", HOFFSET( bolshoi_halo_type, jz ) );
   mem_type.insert( float, "some spin factor", HOFFSET( bolshoi_halo_type, spin ) );
   mem_type.insert( unsigned, "snapshot index", HOFFSET( bolshoi_halo_type, snap_num ) );

   // Create file type.
   file_type.compound( 88 );
   file_type.insert( h5::datatype::std_i64be, "globally unique halo ID", 0 );
   file_type.insert( h5::datatype::std_i64be, "descendant ID", 8 );
   file_type.insert( h5::datatype::std_i32be, "number of progenitors", 16 );
   file_type.insert( h5::datatype::std_i64be, "host halod ID", 20 );
   file_type.insert( h5::datatype::ieee_f32be, "virial mass", 28 );
   file_type.insert( h5::datatype::ieee_f32be, "virial radius", 32 );
   file_type.insert( h5::datatype::ieee_f32be, "vrms", 36 );
   file_type.insert( h5::datatype::ieee_f32be, "maximum velocity", 40 );
   file_type.insert( h5::datatype::ieee_f32be, "x position", 44 );
   file_type.insert( h5::datatype::ieee_f32be, "y position", 48 );
   file_type.insert( h5::datatype::ieee_f32be, "z position", 52 );
   file_type.insert( h5::datatype::ieee_f32be, "x velocity", 56 );
   file_type.insert( h5::datatype::ieee_f32be, "y velocity", 60 );
   file_type.insert( h5::datatype::ieee_f32be, "z velocity", 64 );
   file_type.insert( h5::datatype::ieee_f32be, "x spin", 68 );
   file_type.insert( h5::datatype::ieee_f32be, "y spin", 72 );
   file_type.insert( h5::datatype::ieee_f32be, "z spin", 76 );
   file_type.insert( h5::datatype::ieee_f32be, "some spin factor", 80 );
   file_type.insert( h5::datatype::std_i32be, "snapshot index", 84 );
}
