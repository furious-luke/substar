#include <iostream>
#include <boost/filesystem.hpp>
#include <libhpc/libhpc.hh>
#include "utils.hh"

using namespace hpc;
namespace fs = boost::filesystem;

///
///
///
int
main( int argc,
      char* argv[] )
{
   mpi::initialise( argc, argv );
   LOG_CONSOLE();

   // Get arguments.
   fs::path path;
   if( argc >= 2 )
      path = argv[1];
   else
      path = ".";

   // Load the forests. The file contains rows of
   // forest ID and tree ID, creating a multimap.
   LOGILN( "Loading forests.", setindent( 2 ) );
   multimap<long long,long long> forests_map;
   load_forests_map( path, forests_map );
   LOGILN( "Done.", setindent( -2 ) );

   // Load the precomputed counts of halos to each
   // forest, called forest sizes.
   LOGILN( "Loading forest sizes.", setindent( 2 ) );
   map<long long,unsigned> forest_sizes_map;
   load_forest_sizes_map( ".", forest_sizes_map );
   LOGILN( "Done.", setindent( -2 ) );

   // Load locations.
   LOGILN( "Loading locations.", setindent( 2 ) );
   map<long long,std::pair<unsigned,size_t>> locations_map;
   map<unsigned,string> file_map_map;
   load_locations_map( path, locations_map, file_map_map );
   LOGILN( "Done.", setindent( -2 ) );

#ifndef NDEBUG
   // Locations and forests must be the same size (number of trees).
   ASSERT( forests_map.size() == locations_map.size() );

   // Check that every value in the multimap corresponds
   // to a key in the locations map.
   for( const auto& pair : forests_map )
      ASSERT( locations_map.has( pair.second ) );
#endif

   // Cache the number of forests. Can't get this from
   // the multimap as that will be the number of
   // trees, not forests.
   unsigned num_forests = forest_sizes_map.size();

   // Allocate for the forest sizes array.
   vector<unsigned> forest_sizes( num_forests );

   // Allocate for the forests CSR and construct both
   // the CSR and the forest sizes array.
   LOGILN( "Computing flattened arrays.", setindent( 2 ) );
   csr<long long> forests;
   forests.num_rows( num_forests );
   {
      vector<size_t>::view cnts = forests.counts();
      unsigned idx = 0;
      for( const auto& pair : forest_sizes_map )
      {
	 long long fid = pair.first;
	 cnts[idx] = forests_map.count( fid );
	 forest_sizes[idx] = pair.second;
	 ++idx;
      }
      forests.setup_array( true );
   }
   {
      vector<long long>& array = forests.mod_array();

      // The array needs to be the same size as the forest
      // multimap.
      ASSERT( array.size() == forests_map.size() );

      unsigned idx = 0;
      for( const auto& pair : forests_map )
	 array[idx++] = pair.second;
   }
   LOGILN( "Done.", setindent( -2 ) );

   // Clear unneeded maps now.
   forests_map.clear();
   forest_sizes_map.clear();

   // Now convert the locations.
   vector<size_t> location_offsets( locations_map.size() );
   vector<unsigned> location_files( locations_map.size() );
   vector<string> file_map( file_map_map.size() );
   for( unsigned ii = 0; ii < forests.array().size(); ++ii )
   {
      auto loc = locations_map.get( forests.array()[ii] );
      location_files[ii] = loc.first;
      location_offsets[ii] = loc.second;
      if( file_map[loc.first].empty() )
	 file_map[loc.first] = file_map_map.get( loc.first );

#ifndef NDEBUG
      // Be sure every file reference matches.
      else
	 ASSERT( file_map[loc.first] == file_map_map.get( loc.first ) );
#endif
   }

   // Clear more maps.
   locations_map.clear();
   file_map_map.clear();

   // Dump results as HDF5.
   LOGILN( "Saving results.", setindent( 2 ) );
   h5::file out( "forests.h5", H5F_ACC_TRUNC );
   out.write<long long>( "forest_trees", forests );
   out.write<unsigned>( "forest_halo_counts", forest_sizes );
   out.write<size_t>( "location_file_offsets", location_offsets );
   out.write<unsigned>( "location_file_indices", location_files );
   out.write<string>( "location_files", file_map );
   LOGILN( "Done.", setindent( -2 ) );

   mpi::finalise();
   return EXIT_SUCCESS;
}
