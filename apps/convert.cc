#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <libhpc/libhpc.hh>
#include "types.hh"
#include "exporter.hh"
#include "utils.hh"

using namespace hpc;
namespace fs = boost::filesystem;

extern const unsigned num_threads;
extern unsigned file_offsets[];

///
/// Tunable parameters.
///
const unsigned num_threads = 10;
const float    update_every = 10.0;

///
/// Constants.
///
const float particle_mass = 1.35e8;
const float newton_gravitation = 6.67428e-11;
const float solar_mass = 1.98892e30;
const float m_per_kpc = 3.08568025e19;

///
/// OpenMP shared variables.
///
unsigned file_offsets[num_threads];

///
/// Custom mass comparison. 
///
struct mass_compare_type
{
   mass_compare_type( const vector<bolshoi_halo_type>& forest )
      : forest( forest )
   {
   }

   bool
   operator()( unsigned op_a,
	       unsigned op_b )
   {
      return forest[op_a].mvir > forest[op_b].mvir;
   }

   const vector<bolshoi_halo_type>& forest;
};

// std::ostream&
// operator<<( std::ostream& strm,
//             const bolshoi_halo_type& obj )
// {
//    strm << "{scale: " << obj.scale;
//    strm << ", id: " << obj.id;
//    strm << ", desc_scale: " << obj.desc_scale;
//    strm << ", desc_id: " << obj.desc_id;
//    strm << ", num_prog: " << obj.num_prog;
//    strm << ", pid: " << obj.pid;
//    strm << ", upid: " << obj.upid;
//    strm << ", desc_pid: " << obj.desc_pid;
//    strm << ", phantom: " << obj.phantom;
//    strm << ", sam_mvir: " << obj.sam_mvir;
//    strm << ", mvir: " << obj.mvir;
//    strm << ", rs: " << obj.rs;
//    strm << ", vrms: " << obj.vrms;
//    strm << ", mmp: " << obj.mmp;
//    strm << ", scale_of_last_mm: " << obj.scale_of_last_mm;
//    strm << ", vmax: " << obj.vmax;
//    strm << ", x: " << obj.x;
//    strm << ", y: " << obj.y;
//    strm << ", z: " << obj.z;
//    strm << ", vx: " << obj.vx;
//    strm << ", vy: " << obj.vy;
//    strm << ", vz: " << obj.vz;
//    strm << ", jx: " << obj.jx;
//    strm << ", jy: " << obj.jy;
//    strm << ", jz: " << obj.jz;
//    strm << ", spin: " << obj.spin;
//    // unsigned            breadth_first_id; // Breadth-first ordering of halos within a tree.
//    // unsigned            depth_first_id; // Depth-first ordering of halos withing a tree.
//    // long long  tree_root_id; // ID of the halo at the last timestep in the tree.
//    // long long  orig_halo_id; // Original halo ID from the halo finder.
//    strm << ", snap_num: " << obj.snap_num << "}";
//    // unsigned            next_coprog_depthfirst_id; // Depth-first ID of next coprogenitor.
//    // unsigned            last_coprog_depthfirst_id; // Depth-first ID of last progenitor.
//    // float              rs_klypin; // Scale radius determined using Vmax and Mvir (see Rockstar paper).
//    // float              mvir_all; // Mass enclosed within the specified overdensity, including unbound
//    //                               // particles (Msun/h).
// }

bool
id_map_has( const vector<std::pair<long long,unsigned>>& id_map,
	    long long id )
{
   return std::binary_search( id_map.begin(), id_map.end(), id, hpc::less_first<std::pair<long long,unsigned>>() );
}

unsigned
id_map_get( const vector<std::pair<long long,unsigned>>& id_map,
	    long long id )
{
   auto it = std::lower_bound( id_map.begin(), id_map.end(), id, hpc::less_first<std::pair<long long,unsigned>>() );
   ASSERT( it != id_map.end() );
   return it->second;
}

void
convert_ids( vector<bolshoi_halo_type>& forest,
	     const vector<std::pair<long long,unsigned>>& id_map )
{
   for( auto& halo : forest )
   {
      halo.id = id_map_get( id_map, halo.id );
      if( halo.desc_id != -1 )
	 halo.desc_id = id_map_get( id_map, halo.desc_id );
      if( halo.pid != -1 )
	 halo.pid = id_map_get( id_map, halo.pid );
      // if( halo.upid != -1 )
      // 	 halo.upid = id_map_get( id_map, halo.upid );
      // if( halo.desc_pid != -1 )
      // 	 halo.desc_pid = id_map_get( id_map, halo.desc_pid );
   }
}

long long
top_pid( const vector<bolshoi_halo_type>& forest,
	 long long pid )
{
   long long tmp = pid;
   while( tmp != -1 )
   {
      pid = tmp;
      tmp = forest[pid].pid;
   }
   return pid;
}

void
bolshoi_to_sage( const vector<bolshoi_halo_type>& forest,
		 const csr<unsigned>& progens,
		 const map<long long,unsigned>& fof_group_map,
		 const csr<unsigned>& fof_groups,
		 vector<sage_halo_type>& sage_tree )
{
   LOGDLN( "Converting structures.", setindent( 2 ) );

   // Copy basic info.
   LOGD( "Copying basic information... " );
   for( unsigned ii = 0; ii < forest.size(); ++ii )
   {
      const bolshoi_halo_type& bh = forest[ii];
      sage_halo_type& sh = sage_tree[ii];

      sh.descendant              = bh.desc_id;
      sh.first_progenitor        = -1;
      sh.next_progenitor         = -1;
      sh.first_halo_in_fof_group = -1;
      sh.next_halo_in_fof_group  = -1;
      sh.Len                     = bh.mvir/particle_mass;
      sh.M_Mean200               = 0.0;
      sh.Mvir                    = bh.mvir/1e10;
      sh.M_TopHat                = 0.0;
      sh.Pos[0]                  = bh.x;
      sh.Pos[1]                  = bh.y;
      sh.Pos[2]                  = bh.z;
      sh.Vel[0]                  = bh.vx;
      sh.Vel[1]                  = bh.vy;
      sh.Vel[2]                  = bh.vz;
      sh.VelDisp                 = bh.vrms;
      sh.Vmax                    = bh.vmax;

      // Spin needs some processing.
      float vvir = 1e-3*sqrt( newton_gravitation*bh.mvir*solar_mass/(bh.rvir*m_per_kpc) );
      float spin_amp = bh.spin/sqrt( bh.jx*bh.jx + bh.jy*bh.jy + bh.jz*bh.jz );
      if( spin_amp == spin_amp )
      {
	 sh.Spin[0]              = 1e-3*bh.jx*spin_amp*sqrt( 2.0 )*bh.rvir*vvir;
	 sh.Spin[1]              = 1e-3*bh.jy*spin_amp*sqrt( 2.0 )*bh.rvir*vvir;
	 sh.Spin[2]              = 1e-3*bh.jz*spin_amp*sqrt( 2.0 )*bh.rvir*vvir;
      }
      else
      {
	 sh.Spin[0]              = 0;
	 sh.Spin[1]              = 0;
	 sh.Spin[2]              = 0;
      }

      sh.MostBoundID             = 0;
      sh.SnapNum                 = bh.snap_num;
      sh.FileNr                  = 0;
      sh.SubhaloIndex            = 0;
      sh.SubHalfMass             = 0.0;

      // Sanity checks.
      ASSERT( sh.Spin[0] == sh.Spin[0] );
      ASSERT( sh.Spin[1] == sh.Spin[1] );
      ASSERT( sh.Spin[2] == sh.Spin[2] );
      ASSERT( sh.Mvir == sh.Mvir );
      ASSERT( sh.Pos[0] == sh.Pos[0] );
      ASSERT( sh.Pos[1] == sh.Pos[1] );
      ASSERT( sh.Pos[2] == sh.Pos[2] );
      ASSERT( sh.Vel[0] == sh.Vel[0] );
      ASSERT( sh.Vel[1] == sh.Vel[1] );
      ASSERT( sh.Vel[2] == sh.Vel[2] );
      ASSERT( sh.VelDisp == sh.VelDisp );
      ASSERT( sh.Vmax == sh.Vmax );
   }
   LOGDLN( "done." );

   // Process progenitor information.
   LOGD( "Processing progenitor information... " );
   for( unsigned ii = 0; ii < progens.num_rows(); ++ii )
   {
      const auto& row = progens[ii];
      if( row.size() )
      {
	 ASSERT( sage_tree[ii].first_progenitor == -1 );
	 sage_tree[ii].first_progenitor = row[0];
	 for( unsigned jj = 1; jj < row.size(); ++jj )
	 {
	    ASSERT( sage_tree[row[jj - 1]].next_progenitor == -1 );
	    sage_tree[row[jj - 1]].next_progenitor = row[jj];
	 }
      }
   }
   LOGDLN( "done." );

   // Process FOF group information.
   LOGD( "Processing FOF group information... " );
   for( const auto& fof_pair : fof_group_map )
   {
      long long pid = fof_pair.first;
      unsigned ii = fof_pair.second;
      const auto& row = fof_groups[ii];

      // First halo in FOF group should be the most massive, i.e. the
      // first halo in any group should point to itself.
      ASSERT( forest[row[0]].id == pid &&
	      (forest[row[0]].pid == -1 || forest[row[0]].pid == forest[row[0]].id ) );

      for( unsigned jj = 0; jj < row.size(); ++jj )
      {
	 ASSERT( sage_tree[row[jj]].first_halo_in_fof_group == -1 );
	 sage_tree[row[jj]].first_halo_in_fof_group = pid;
      }
      for( unsigned jj = 1; jj < row.size(); ++jj )
      {
	 ASSERT( sage_tree[row[jj - 1]].next_halo_in_fof_group == -1 );
	 sage_tree[row[jj - 1]].next_halo_in_fof_group = row[jj];
      }
   }
   LOGDLN( "done." );

#ifndef NDEBUG
   LOGDLN( "Checking sanity of conversion.", setindent( 2 ) );

   // Trees must be the same size.
   ASSERT( sage_tree.size() == forest.size() );

   // Check progenitors against dscendents.
   LOGDLN( "Checking progenitors.", setindent( 2 ) );
   for( unsigned ii = 0; ii < sage_tree.size(); ++ii )
   {
      int prog = sage_tree[ii].first_progenitor;
      unsigned num_prog = 0;
      while( prog != -1 )
      {
	 LOGDLN( "Progenitor with local ID ", prog, setindent( 2 ) );

	 int desc = sage_tree[prog].descendant;

	 // Must be self-consistent.
	 ASSERT( desc == ii );
	 LOGDLN( "Is self-consistent." );

	 // Must match Bolshoi.
	 ASSERT( forest[prog].desc_id != -1 );
	 ASSERT( forest[prog].desc_id == desc );
	 LOGDLN( "Matches Bolshoi." );

	 // Advance.
	 prog = sage_tree[prog].next_progenitor;
	 ++num_prog;

	 LOGD( setindent( -2 ) );
      }

      // Counts must match.
      ASSERT( num_prog == forest[ii].num_prog );
   }
   LOGDLN( "Done.", setindent( -2 ) );

   // Check FOF groups.
   LOGDLN( "Checking FOF groups.", setindent( 2 ) );
   for( unsigned ii = 0; ii < sage_tree.size(); ++ii )
   {
      int fof = sage_tree[ii].first_halo_in_fof_group;
      LOGDLN( "FOF with local halo ID ", fof );

      // In SAGE every halo belongs to a FOF group.
      ASSERT( fof >= 0 );

      // Single FOF groups will have -1 in pid or pid == id.
      bool single = (sage_tree[fof].next_halo_in_fof_group == -1);
      if( single )
      {
	 LOGDLN( "Single FOF group." );
      	 ASSERT( forest[ii].pid == -1 || forest[ii].pid == forest[ii].id );
      }

      // Otherwise pid must match fof or ID must match fof.
      else
      {
	 LOGDLN( "Multi FOF group." );
	 if( forest[ii].pid != -1 )
	    ASSERT( forest[ii].pid == fof );
	 else
	    ASSERT( forest[ii].id == fof );
      }
   }
   LOGDLN( "Done.", setindent( -2 ) );

   LOGDLN( "Done.", setindent( -2 ) );
#endif

   LOGDLN( "Done.", setindent( -2 ) );
}

void
build_progenitors( const vector<bolshoi_halo_type>& forest,
		   csr<unsigned>& progens )
{
   LOG_ENTER();
   LOGDLN( "Building progenitors.", setindent( 2 ) );

   progens.clear();

   // First construct the counts and setup the arrays.
   progens.num_rows( forest.size() );
   {
      auto cnts = progens.counts();
      unsigned ii = 0;
      for( const auto& bh : forest )
	 cnts[ii++] = bh.num_prog;
      progens.setup_array( true );
   }

   // Fill each array.
   vector<unsigned> cnts( forest.size() );
   std::fill( cnts.begin(), cnts.end(), 0 );
   for( const auto& bh : forest )
   {
      if( bh.desc_id != -1 )
      {
	 unsigned id = bh.desc_id;
	 progens( id, cnts[id]++ ) = bh.id;
      }
   }

   // Order the progenitors for each halo with most massive first.
   mass_compare_type cmp( forest );
   for( unsigned ii = 0; ii < progens.num_rows(); ++ii )
   {
      auto row = progens[ii];
      std::sort( row.begin(), row.end(), cmp );
   }

#ifndef NDEBUG
   // Check each set of progenitors.
   LOGDLN( "Checking consistency.", setindent( 2 ) );
   for( unsigned ii = 0; ii < progens.num_rows(); ++ii )
   {
      auto row = progens[ii];
      LOGD( "Halo ", ii, " with ", row.size(), " progenitors... " );

      // No duplicates.
      set<unsigned> dup;
      dup.insert( row.begin(), row.end() );
      ASSERT( dup.size() == row.size() );

      // Descendants must match.
      for( auto prog : row )
	 ASSERT( forest[prog].desc_id == ii );

      // Must be in most massive order.
      for( unsigned jj = 1; jj < row.size(); ++jj )
	 ASSERT( forest[row[jj]].mvir <= forest[row[jj - 1]].mvir );

      LOGDLN( "okay." );
   }
   LOGDLN( "Done.", setindent( -2 ) );
#endif

   LOGDLN( "Done.", setindent( -2 ) );
   LOG_EXIT();
}

void
build_fof_groups( const vector<bolshoi_halo_type>& forest,
		  map<long long,unsigned>& fof_group_map,
		  csr<unsigned>& fof_groups )
{
   LOG_ENTER();
   LOGDLN( "Building FOF groups.", setindent( 2 ) );

   fof_group_map.clear();
   fof_groups.clear();

   // I need a mapping from the UPID identifying the most massive
   // halo in a FOF group (the parent halo) to the local FOF group
   // index.
   for( const auto& halo : forest )
   {
      long long pid = (halo.pid != -1) ? halo.pid : halo.id;
      auto it = fof_group_map.insert( pid );
      if( it.second )
	 it.first->second = fof_group_map.size() - 1;
   }
   LOGDLN( "Found ", fof_group_map.size(), " FOF groups." );

   // First construct the counts and setup the arrays.
   fof_groups.num_rows( fof_group_map.size() );
   {
      auto cnts = fof_groups.counts();
      std::fill( cnts.begin(), cnts.end(), 0 );
      for( const auto& bh : forest )
      {
	 long long pid = (bh.pid != -1) ? bh.pid : bh.id;
	 unsigned id = fof_group_map.get( pid );
	 cnts[id]++;
      }
      LOGDLN( "FOF group counts: ", cnts );
      fof_groups.setup_array( true );
   }

   // Fill each array.
   vector<unsigned> cnts( fof_group_map.size() );
   std::fill( cnts.begin(), cnts.end(), 0 );
   for( const auto& bh : forest )
   {
      long long pid = (bh.pid != -1) ? bh.pid : bh.id;
      unsigned id = fof_group_map.get( pid );
      fof_groups( id, cnts[id]++ ) = bh.id;
   }

   // Order the FOF groups for each halo with most massive first.
   mass_compare_type cmp( forest );
   for( unsigned ii = 0; ii < fof_groups.num_rows(); ++ii )
   {
      auto row = fof_groups[ii];
      std::sort( row.begin(), row.end(), cmp );
   }

   // There are cases where the most massive FOF is not
   // the host halo. We need to find these cases and swap
   // them.
   for( auto fof_pair : fof_group_map )
   {
      unsigned fof_id = fof_pair.first;
      unsigned ii = fof_pair.second;
      auto row = fof_groups[ii];

      if( forest[row[0]].id != fof_id )
      {
	 unsigned jj = 1;
	 for( ; jj < row.size(); ++jj )
	 {
	    if( forest[row[jj]].id == fof_id )
	       break;
	 }
	 ASSERT( jj < row.size() );
	 LOGDLN( "Swapping misplaced host halo FOF index." );
	 std::swap( row[0], row[jj] );
      }
   }

#ifndef NDEBUG
   // Check each set of FOF halos.
   LOGDLN( "Checking consistency.", setindent( 2 ) );
   for( auto fof_pair : fof_group_map )
   {
      unsigned fof_id = fof_pair.first;
      unsigned ii = fof_pair.second;
      auto row = fof_groups[ii];
      LOGDLN( "FOF group ", ii, " with ", row.size(), " halos.", setindent( 2 ) );
      LOGDLN( "upid: ", fof_id );
      LOGD( "FOF IDs: [", forest[row[0]].id );
      for( unsigned jj = 1; jj < row.size(); ++jj )
	 LOGD( ", ", forest[row[jj]].id );
      LOGDLN( "]" );
      LOGD( "FOF masses: [", forest[row[0]].mvir );
      for( unsigned jj = 1; jj < row.size(); ++jj )
	 LOGD( ", ", forest[row[jj]].mvir );
      LOGDLN( "]" );

      // No duplicates.
      set<unsigned> dup;
      dup.insert( row.begin(), row.end() );
      ASSERT( dup.size() == row.size() );

      // FOF entries must match up.
      for( auto halo : row )
      {
	 if( forest[halo].pid != -1 )
	    ASSERT( forest[halo].pid == fof_id );
	 else
	    ASSERT( forest[halo].id == fof_id );
      }

      // // Must be in most massive order. Unfortunately this won't be
      // // true now that I'm reordering to fix Bolshoi problems.
      // for( unsigned jj = 1; jj < row.size(); ++jj )
      // 	 ASSERT( forest[row[jj]].mvir <= forest[row[jj - 1]].mvir );

      // The first halo in the FOF group must have the same ID as
      // the FOF group has been given.
      ASSERT( forest[row[0]].id == fof_id );

      LOGDLN( "Done.", setindent( -2 ) );
   }
   LOGDLN( "Done.", setindent( -2 ) );
#endif

   LOGDLN( "Done.", setindent( -2 ) );
   LOG_EXIT();
}

///
/// Finish processing a single Bolshoi forest.
///
void
process_forest( vector<bolshoi_halo_type>& forest,
		vector<std::pair<long long,unsigned>>& id_map,
		exporter& exp )
{
   LOGDLN( "Post processing forest.", setindent( 2 ) );

   // First need to prepare the ID map by sorting it.
   // I do this to use a binary search.
   std::sort( id_map.begin(), id_map.end(), hpc::less_first<std::pair<long long,unsigned>>() );

#ifndef NDEBUG
   // Check that the ID map and the tree are consistent.
   for( auto ids : id_map )
      ASSERT( ids.first == forest[ids.second].id );

   for( const auto& halo : forest )
   {
      // Check that each referenced ID is in this forest.
      ASSERT( id_map_has( id_map, halo.id ) );
      if( halo.desc_id != -1 )
	 ASSERT( id_map_has( id_map, halo.desc_id ) );
      if( halo.pid != -1 )
	 ASSERT( id_map_has( id_map, halo.pid ) );
      // if( halo.upid != -1 )
      // 	 ASSERT( id_map_has( id_map, halo.upid ) );
      // if( halo.desc_pid != -1 )
      // 	 ASSERT( id_map_has( id_map, halo.desc_pid ) );

      // Check that the x,y,z position is in range.
      ASSERT( halo.x >= 0.0 && halo.x <= 250.0 );
      ASSERT( halo.y >= 0.0 && halo.y <= 250.0 );
      ASSERT( halo.z >= 0.0 && halo.z <= 250.0 );
   }

   // Walk each halo's host halo chain and make sure all masses are larger
   // the further up the chain we go. Also check that host and subhalos
   // exist in the same snapshot.
   for( const auto& halo : forest )
   {
      long long pid = halo.pid;
      float mass = halo.mvir;
      unsigned snap_num = halo.snap_num;
      while( pid != -1 )
      {
   	 unsigned id = id_map_get( id_map, pid );

// 	 // Need to confirm that each host halo has greater mass than
// 	 // all subhalos, as these subhalos are included in the mass
// 	 // of parent halos.
//    	 // ASSERT( forest[id].mvir >= mass );
// 	 if( forest[id].mvir < mass )
// 	 {
// 	    LOGDLN( "Found a host halo with lower mass than subhalo." );
// #pragma omp critical
// 	    {
// 	       std::ofstream bmf( "bad_mass.dat", std::ios::out | std::ios::app );
// 	       ASSERT( bmf );
// 	       bmf << pid << "  " << forest[id].mvir << "  " << halo.id << "  " << mass << "\n";
// 	    }
// 	 }

	 // Snapshot numbers must match parents.
	 ASSERT( forest[id].snap_num == snap_num );

	 // Move to parent.
   	 pid = forest[id].pid;
      }
   }
#endif

   // First, convert all the IDs so we can free the ID map, which can be huge.
   convert_ids( forest, id_map );
   id_map.deallocate();

   // Convert all the PID values to their actual top most containing halo.
   // I need to do this because it seems PID does not actually point to
   // this, and I need it to.
   for( auto& halo : forest )
      halo.pid = top_pid( forest, halo.pid );

   // Build the progenitors.
   csr<unsigned> progens;
   build_progenitors( forest, progens );

   // Build the FOF groups.
   map<long long,unsigned> fof_group_map;
   csr<unsigned> fof_groups;
   build_fof_groups( forest, fof_group_map, fof_groups );

   // Convert to sage format.
   vector<sage_halo_type> sage_tree( forest.size() );
   bolshoi_to_sage( forest, progens, fof_group_map, fof_groups, sage_tree );

   // Save current sage tree.
   exp.export_forest( sage_tree );

   LOGDLN( "Done.", setindent( -2 ) );
}

///
///
///
int
main( int argc,
      char* argv[] )
{
   mpi::initialise( argc, argv );
   LOG_CONSOLE();
   // LOG_PUSH( new logging::file( "log.2", logging::info ) );
   // LOG_PUSH( new logging::stdout( logging::info ) );
   // LOG_PUSH( new logging::omp::file( "conv.", logging::debug ) );

   // Get arguments.
   fs::path path;
   if( argc >= 2 )
      path = argv[1];
   else
      path = ".";

   // Keep track of how many Halos each thread has seen
   // for updating purposes.
   unsigned long long halos_seen[num_threads];
   unsigned long long forests_seen[num_threads];
   unsigned long long total_halos_seen, total_forests_seen;
   posix::time_type since_update = posix::timer();
   std::fill( halos_seen, halos_seen + num_threads, 0 );
   std::fill( forests_seen, forests_seen + num_threads, 0 );

   // Load in tables needed to locate FOF groups.
   LOGILN( "Reading auxilliary tables.", setindent( 2 ) );
   vector<size_t> forest_tree_displs, forest_tree_counts;
   vector<unsigned> forest_halo_cnts;
   vector<size_t> location_file_offsets;
   vector<unsigned> location_file_idxs;
   vector<string> location_files;
   {
      h5::file aux_file( "forests.h5", H5F_ACC_RDONLY );
      aux_file.reada<size_t>( "forest_trees_displs", forest_tree_displs );
      aux_file.reada<size_t>( "forest_trees_counts", forest_tree_counts );
      aux_file.reada<unsigned>( "forest_halo_counts", forest_halo_cnts );
      aux_file.reada<size_t>( "location_file_offsets", location_file_offsets );
      aux_file.reada<unsigned>( "location_file_indices", location_file_idxs );
      aux_file.reada( "location_files", location_files );
   }
   LOGILN( "Done.", setindent( -2 ) );

#ifndef NDEBUG
   // Check sizes.
   ASSERT( forest_tree_displs.size() == forest_halo_cnts.size() );
   ASSERT( forest_tree_counts.size() == forest_halo_cnts.size() );
#endif

   // Cache the number of forests.
   unsigned num_forests = forest_halo_cnts.size();
   LOGILN( "Have ", num_forests, " forests to process." );

   // Setup the number of threads to use.
#ifdef _OPENMP
   omp_set_num_threads( num_threads );
   LOGILN( "Running with ", num_threads, " threads." );
#endif

   // Calculate for each thread the start and finish offsets. I will
   // do this by balancing the number of halos (approximately) each
   // will have to process. First count the total number of halos
   // we have.
   LOGILN( "Distributing.", setindent( 2 ) );
   long long net_halos = 0;
   for( unsigned ii = 0; ii < num_forests; ++ ii )
      net_halos += forest_halo_cnts[ii];
   long long halos_per_thread = net_halos/num_threads;
   LOGILN( net_halos, " total halos, giving ~", halos_per_thread, " per thread." );

   // Now add forests to each thead until they reach the balance point.
   unsigned start_forest[num_threads], finish_forest[num_threads];
   long long num_thread_halos[num_threads], num_thread_forests[num_threads];
   for( unsigned ii = 0; ii < num_threads; ++ii )
   {
      // First thread begins at zero, other threads
      // begin where the prior thread left off.
      if( ii == 0 )
	 start_forest[ii] = 0;
      else
	 start_forest[ii] = finish_forest[ii - 1];

      // Build up forests on each thread based on how many halos
      // they contribute.
      num_thread_forests[ii] = 0;
      num_thread_halos[ii] = 0;
      finish_forest[ii] = start_forest[ii];
      while( finish_forest[ii] != num_forests && num_thread_halos[ii] < halos_per_thread )
      {
	 num_thread_halos[ii] += forest_halo_cnts[finish_forest[ii]++];
	 ++num_thread_forests[ii];
      }
      LOGILN( "Thread ", ii, " has ", num_thread_forests[ii], " forests and ", num_thread_halos[ii], " halos."  );
      LOGILN( "And will process forests from ", start_forest[ii], " to ", finish_forest[ii], "." );
   }
   LOGILN( "Done.", setindent( -2 ) );

   // Form a team of threads here.
#pragma omp parallel
   {
      // Cache our thread ID.
      int tid = OMP_TID;

      // Create the exporter.
      exporter exp( start_forest[tid], finish_forest[tid], forest_halo_cnts );

      // Loop over our section of the forests.
      for( unsigned ii = start_forest[tid]; ii < finish_forest[tid]; ++ii )
      {
	 LOGDLN( "Processing forest ", ii, ".", setindent( 2 ) );

	 // Setup tree displacement and counts.
	 unsigned tree_displ = forest_tree_displs[ii];
	 unsigned num_trees = forest_tree_counts[ii];
	 LOGDLN( num_trees, " trees in forest." );
	 LOGDLN( "Looking at tree displacement ", tree_displ );

	 // How many halos in this forest?
	 unsigned num_halos = forest_halo_cnts[ii];
	 LOGDLN( num_halos, " halos in forest." );

	 // Use a mapping to map from global Bolshoi tree IDs to
	 // local tree IDs for SAGE.
	 vector<std::pair<long long,unsigned>> id_map( num_halos );

	 // Store current forest in a vector of halos.
	 vector<bolshoi_halo_type> forest( num_halos );

	 // Track how many halos we load in this forest.
	 unsigned cur_halo = 0;

	 // Select the range of trees we need to process and
	 // iterate over each tree.
	 for( unsigned jj = 0; jj < num_trees; ++jj )
	 {
	    LOGDLN( "Processing tree ", jj, ".", setindent( 2 ) );

	    // Map to the location information.
	    auto file_offs = location_file_offsets[tree_displ + jj];
	    LOGDLN( "File offset ", file_offs, "." );
	    LOGDLN( "File index ", location_file_idxs[tree_displ + jj], "." );
	    auto filename = location_files[location_file_idxs[tree_displ + jj]];

	    // Open the file in which the tree lives then scan to
	    // the appropriate offset for the tree.
	    LOGDLN( "Opening file \"", filename, "\"" );
	    std::ifstream file( (path/filename).c_str(), std::ios::in );
	    ASSERT( file );
	    LOGDLN( "Seeking to ", file_offs, "." );
	    file.seekg( file_offs );
	    ASSERT( file );

	    // We aren't given information about how many halos are in
	    // each tree, so we just have to loop until we hit the end
	    // of the file.
	    while( file.peek() != '#' && !file.eof()  )
	    {
	       LOGDLN( "Processing halo ", cur_halo, "." );

	       // Put this halo on the back of the current tree.
	       bolshoi_halo_type& bh = forest[cur_halo];

	       // Read each field in turn.
	       float scale, desc_scale, sam_mvir, rs, scale_of_last_mm;
	       long long upid, desc_pid, tree_root_id, orig_halo_id;
	       unsigned breadth_first_id, depth_first_id;
	       int phantom, mmp;
	       file >> scale >> bh.id >> desc_scale >> bh.desc_id >> bh.num_prog >> bh.pid
	    	    >> upid >> desc_pid >> phantom >> sam_mvir >> bh.mvir
	    	    >> bh.rvir >> rs >> bh.vrms >> mmp >> scale_of_last_mm >> bh.vmax
	    	    >> bh.x >> bh.y >> bh.z >> bh.vx >> bh.vy >> bh.vz >> bh.jx >> bh.jy
	    	    >> bh.jz >> bh.spin >> breadth_first_id >> depth_first_id
	    	    >> tree_root_id >> orig_halo_id >> bh.snap_num;
	       finish_line( file );
	       ASSERT( !file.fail() );

	       // Insert the ID into the ID map.
	       id_map[cur_halo].first = bh.id;
	       id_map[cur_halo].second = cur_halo;

	       // Advance.
	       ++cur_halo;
	    }

	    LOGDLN( "Done.", setindent( -2 ) );
	 }
	 LOGDLN( "Loaded ", cur_halo, " halos." );

	 // Must have filled the entire forest.
	 ASSERT( cur_halo == num_halos );
// 	 if( cur_halo != num_halos )
// 	 {
// #pragma omp critical( dummy )
// 	    LOGILN( "Wrong counts: ", cur_halo, ", ", num_halos );
// 	 }

	 // Process the forest.
	 process_forest( forest, id_map, exp );

	 // Update the number of halos seen.
	 halos_seen[tid] += cur_halo;
	 ++forests_seen[tid];
	 LOGDLN( "Seen ", halos_seen[tid], " halos." );

	 // Check if we need to update the user.
#pragma omp critical( update )
	 if( posix::seconds( posix::timer() - since_update ) > update_every )
	 {
	    // Sum the net results.
	    total_halos_seen = halos_seen[0];
	    total_forests_seen = forests_seen[0];
	    for( unsigned jj = 1; jj < num_threads; ++jj )
	    {
	       total_halos_seen += halos_seen[jj];
	       total_forests_seen += forests_seen[jj];
	    }

	    // Print some info.
	    LOGILN( "Total halos seen:    ", total_halos_seen );
	    LOGILN( "Total forests seen:  ", total_forests_seen );
	    LOGILN( "Percentage complete: ", 100.0*(float)total_forests_seen/(float)num_forests );

	    // Reset the update clock.
	    since_update = posix::timer();
	 }

	 LOGDLN( "Done.", setindent( -2 ) );
      }

#pragma omp critical( done )
      LOGDLN( "Thread ", tid, " done." );
   }

   LOGILN( "Complete." );
   mpi::finalise();
   return EXIT_SUCCESS;
}
