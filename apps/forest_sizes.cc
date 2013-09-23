#include <iostream>
#include <boost/filesystem.hpp>
#include <libhpc/libhpc.hh>
#include "substar/utils.hh"

using namespace hpc;
namespace fs = boost::filesystem;

extern const unsigned num_threads;
extern unsigned file_offsets[];

///
/// Tunable parameters.
///
const double update_every = 10; // seconds

///
///
///
int
main( int argc,
      char* argv[] )
{
   LOG_CONSOLE();

   // Get arguments.
   fs::path path;
   if( argc >= 2 )
      path = argv[1];
   else
      path = ".";

   // Use this to track when to update.
   time_type since_update = timer();

   // Reset the output file.
   {
      std::ofstream cnt_file( "forest_sizes.dat", std::ios::out );
   }

   // Load the forests. The file contains rows of
   // forest ID and tree ID, creating a multimap.
   LOGILN( "Loading forests.", setindent( 2 ) );
   multimap<long long,long long> forests_map;
   load_forests_map( path, forests_map );
   LOGILN( "Done.", setindent( -2 ) );

   // Load the forest locations.
   LOGILN( "Loading forest locations.", setindent( 2 ) );
   map<long long,std::pair<unsigned,size_t>> locations_map;
   map<unsigned,string> file_map;
   load_locations_map( path, locations_map, file_map );
   LOGILN( "Done.", setindent( -2 ) );

#ifndef NDEBUG
   // Locations and forests must be the same size (number of trees).
   ASSERT( forests_map.size() == locations_map.size() );

   // Check that every value in the multimap corresponds
   // to a key in the locations map.
   for( const auto& pair : forests_map )
      ASSERT( locations_map.has( pair.second ) );
#endif

   // Calculate a set of forest IDs.
   LOGILN( "Calculating forest IDs.", setindent( 2 ) );
   set<long long> forest_ids;
   for( const auto& pair : forests_map )
      forest_ids.insert( pair.first );
   LOGILN( "Done.", setindent( -2 ) );

   // Cache the number of forests.
   unsigned num_forests = forest_ids.size();
   LOGILN( "Have ", num_forests, " to process." );

   // Use a shared iterator to indicate which
   // forest we are up to.
   auto forest_it = forest_ids.cbegin();

   // Keep track of how many forests we've completed.
   unsigned forests_done = 0;

   // Start our workers.
   LOGILN( "Launching ", num_threads, " threads." );
#ifdef _OPENMP
   omp_set_num_threads( num_threads );
#endif
#pragma omp parallel
   {
      // Enter the loop.
      while( forest_it != forest_ids.end() )
      {
	 // Cache the forest we will work on next.
	 set<long long>::const_iterator my_forest_it;
#pragma omp critical( get_forest )
	 my_forest_it = forest_it++;
	 long long my_forest = *my_forest_it;
	 LOGDLN( "Working on forest ", my_forest );

	 // Store a local count.
	 unsigned num_halos = 0;

	 // Iterate over the trees in this forest.
	 auto tree_rng = forests_map.equal_range( my_forest );
	 while( tree_rng.first != tree_rng.second )
	 {
	    // Cache the tree index.
	    auto tree = (*tree_rng.first).second;

	    // Map to the location information.
	    auto loc = locations_map.get( tree );

	    // Open the file in which the tree lives then scan to
	    // the appropriate offset for the tree.
	    LOGDLN( "Opening file \"", file_map.get( loc.first ), "\"." );
	    std::ifstream file( (path/file_map.get( loc.first )).c_str(), std::ios::in );
	    ASSERT( file );
	    LOGDLN( "Seeking to ", loc.second, "." );
	    file.seekg( loc.second );
	    ASSERT( file );

	    // We aren't given information about how many halos are in
	    // each tree, so we just have to loop until we hit the end
	    // of the file.
#ifndef NDEBUG
	    bool first = true;
#endif
	    while( file.peek() != '#' && !file.eof()  )
	    {
#ifndef NDEBUG
	       // The first line should have an ID that matches the
	       // tree root ID we're after.
	       if( first )
	       {
		  float scale;
		  long long id;
		  file >> scale >> id;
		  ASSERT( !file.fail() );
		  ASSERT( id == tree );
		  first = false;
	       }
#endif

	       // Read over the line.
	       finish_line( file );
	       ASSERT( !file.fail() );

	       // Add a halo.
	       ++num_halos;
	    }

	    // Advance the trpee iterator.
	    ++tree_rng.first;
	 }
	 LOGDLN( "Found ", num_halos, " halos." );

	 // Update the number of forests we've done.
#pragma omp critical( increment_done )
	 ++forests_done;

	 // Write out to file the halo count for this forest.
#pragma omp critical( write_out )
	 {
	    std::ofstream cnt_file( "forest_sizes.dat", std::ios::out | std::ios::app );
	    cnt_file << my_forest << "  " << num_halos << "\n";
	 }

	 // Have each thread dump information every so often.
#pragma omp critical( update )
	 if( seconds( timer() - since_update ) > update_every )
	 {
	    LOGILN( forests_done, " of ", num_forests, " completed." );

	    // Reset the update clock.
	    since_update = timer();
	 }
      }

      LOGDLN( "Shutting down thread." );
   }

   LOGILN( "Finished." );
   return EXIT_SUCCESS;
}
