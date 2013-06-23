#include <fstream>
#include "utils.hh"

using namespace hpc;
namespace fs = boost::filesystem;

///
///
///
void
load_forests_map( const fs::path& path,
		  multimap<long long,long long>& forests )
{
   LOG_ENTER();

   // Open the gzip stream.
   std::ifstream file( (path/"forests.list").c_str(), std::ios_base::in );
   ASSERT( file );

   // There will be at least one commented line at the start.
   skip_comments( file );

   // Then each line will have a tree ID and a forest ID.
   forests.clear();
   while( file.peek() != '\n' && !file.eof() )
   {
      long long tree, forest;
      file >> tree >> forest;
      finish_line( file );
      ASSERT( !file.fail() );
      forests.insert( forest, tree );
   }

   LOG_EXIT();
}

///
///
///
void
load_locations_map( const fs::path& path,
		    map<long long,std::pair<unsigned,size_t>>& locations,
		    map<unsigned,string>& file_map )
{
   LOG_ENTER();

   // Clear structures.
   locations.clear();
   file_map.clear();

   // Open the gzip stream.
   std::ifstream file( (path/"locations.dat").c_str(), std::ios_base::in );
   ASSERT( file );

   // There will be at least one commented line at the start.
   skip_comments( file );

   // Then each line will have a tree ID and a forest ID.
   while( file.peek() != '\n' && !file.eof() )
   {
      long long tree, file_id;
      size_t offset;
      string filename;
      file >> tree >> file_id >> offset >> filename;
      finish_line( file );
      ASSERT( !file.fail() );
      locations.insert( tree, std::make_pair( file_id, offset ) );

#ifndef NDEBUG
      // If the file map has the id, the filename must match.
      if( file_map.has( file_id ) )
	 ASSERT( file_map.get( file_id ) == filename );
#endif

      file_map.insert( file_id, filename );
   }

   LOG_EXIT();
}

///
///
///
void
load_forest_sizes_map( const fs::path& path,
		       map<long long,unsigned>& forest_sizes )
{
   LOG_ENTER();

   forest_sizes.clear();

   // Open the gzip stream.
   std::ifstream file( (path/"forest_sizes.dat").c_str(), std::ios_base::in );
   ASSERT( file );

   // Then each line will have a tree ID and a forest ID.
   while( file.peek() != '\n' && !file.eof() )
   {
      long long forest;
      unsigned cnt;
      file >> forest >> cnt;
      finish_line( file );
      ASSERT( !file.fail() );
      forest_sizes.insert( forest, cnt );
   }

   LOG_EXIT();
}
