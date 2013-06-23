#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <libhpc/libhpc.hh>
#include "types.hh"
#include "utils.hh"

using namespace hpc;

///
/// Tunable parameters.
///
const unsigned num_threads = 10;
const float    update_every = 10.0;

///
///
///
int
main( int argc,
      char* argv[] )
{
   mpi::initialise( argc, argv );
   LOG_PUSH( new logging::file( "log", logging::info ) );
   // LOG_PUSH( new logging::stdout( logging::info ) );
   // LOG_PUSH( new logging::omp::file( "conv.", logging::debug ) );

   // Keep track of how many Halos each thread has seen
   // for updating purposes.
   unsigned long long halos_seen[num_threads];
   unsigned long long forests_seen[num_threads];
   unsigned long long total_halos_seen, total_forests_seen;
   unix::time_type since_update = unix::timer();
   std::fill( halos_seen, halos_seen + num_threads, 0 );
   std::fill( forests_seen, forests_seen + num_threads, 0 );

   // Load in tables needed to locate FOF groups.
   LOGILN( "Reading auxilliary tables.", setindent( 2 ) );
   csr<long long> forest_trees;
   vector<long long> forest_ids;
   vector<unsigned> forest_halo_cnts;
   vector<size_t> location_file_offsets;
   vector<unsigned> location_file_idxs;
   vector<string> location_files;
   {
      h5::file aux_file( "bolshoi.h5", H5F_ACC_RDONLY );
      aux_file.reada<long long>( "/info/forest_trees", forest_trees );
      aux_file.reada<long long>( "/info/forest_ids", forest_ids );
      aux_file.reada<unsigned>( "/info/forest_halo_counts", forest_halo_cnts );
      aux_file.reada<size_t>( "/info/location_file_offsets", location_file_offsets );
      aux_file.reada<unsigned>( "/info/location_file_indices", location_file_idxs );
      aux_file.reada( "/info/location_files", location_files );
   }
   LOGILN( "Done.", setindent( -2 ) );

#ifndef NDEBUG
   // Check sizes.
   ASSERT( forest_tree.displs().size() == forest_halo_cnts.size() + 1 );
   ASSERT( forest_tree.array().size() == location_file_offsets.size() );
   ASSERT( forest_tree.array().size() == location_file_idxs.size() );
#endif

   // Cache the number of forests.
   unsigned num_forests = forest_halo_cnts.size();
   LOGILN( "Have ", num_forests, " forests to process." );

   // Setup the number of threads to use.
   omp_set_num_threads( num_threads );
   LOGILN( "Running with ", num_threads, " threads." );

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

   // Open the H5 file.
   h5::file h5_file( "forests.h5", H5F_FILE_TRUNC );
   h5::group h5_group( h5_file, "forests" );

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

	 // Cache forest ID.
	 long long forest_id = forest_ids[ii];

	 // Setup tree displacement and counts.
	 unsigned tree_displ = forest_tree.displs()[ii];
	 unsigned num_trees = forest_tree.displs()[ii + 1] - tree_displ;
	 LOGDLN( num_trees, " trees in forest." );
	 LOGDLN( "Looking at tree displacement ", tree_displ );

	 // How many halos in this forest?
	 unsigned num_halos = forest_halo_cnts[ii];
	 LOGDLN( num_halos, " halos in forest." );

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
	    std::ifstream file( filename, std::ios::in );
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

	       // Advance.
	       ++cur_halo;
	    }

	    LOGDLN( "Done.", setindent( -2 ) );
	 }
	 LOGDLN( "Loaded ", cur_halo, " halos." );

	 // Must have filled the entire forest.
	 ASSERT( cur_halo == num_halos );

	 // Write out the forest.
	 {
	    h5::dataspace dspace;
	    dspace.create( num_halos );
	    h5::dataset dset( h5_group, to_string( forest_id ), file_type, dspace );
	    dset.write( forest.data(), mem_type );
	 }

	 // Update the number of halos seen.
	 halos_seen[tid] += cur_halo;
	 ++forests_seen[tid];
	 LOGDLN( "Seen ", halos_seen[tid], " halos." );

	 // Check if we need to update the user.
#pragma omp critical( update )
	 if( unix::seconds( unix::timer() - since_update ) > update_every )
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
	    since_update = unix::timer();
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
