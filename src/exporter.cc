#include "exporter.hh"

using namespace hpc;

extern const unsigned num_threads;
extern unsigned file_offsets[];

const unsigned exporter::default_halos_per_file = 150000;

exporter::exporter( unsigned first_forest,
		    unsigned last_forest,
		    const vector<unsigned>& forest_halo_cnts,
		    unsigned halos_per_file )
  : _hpf( halos_per_file ),
    _first( first_forest ),
    _last( last_forest ),
    _halo_cnts( forest_halo_cnts ),
    _num_halos( 0 ),
    _next_forests( 0 ),
    _next_halos( 0 )
{
   _calc_file_offset();
}

exporter::~exporter()
{
}

void
exporter::export_forest( const vector<sage_halo_type>& forest )
{
   LOGDLN( "Exporting forest.", setindent( 2 ) );

   // If _next_forests is zero we need to calculate the next file.
   if( !_next_forests )
   {
      LOGDLN( "Calculating size of next file, index ", _cur_file, ".", setindent( 2 ) );
      _calc_next_count( _first );
      _open_file();
      _write_header();
      ++_cur_file;
      LOGDLN( "done." );
   }

   // Write forest.
   ASSERT( _first != _last );
   LOGT( "Writing forest... " );
   _file.write( (const char*)forest.data(), sizeof(sage_halo_type)*forest.size() );
   ASSERT( _file );
   LOGTLN( "done." );

   // Update my counts.
   --_next_forests;
   LOGTLN( "Updated the remaining forests in this file to ", _next_forests, "." );
   _num_halos += forest.size();
   LOGTLN( _num_halos, " halos written so far for this thread." );
   ++_first;
}

unsigned
exporter::_calc_next_count( unsigned first )
{
   for( _next_forests = 0, _next_halos = 0;
	_next_halos < _hpf && first != _last;
	++_next_forests, ++first )
   {
      _next_halos += _halo_cnts[first];
   }
   LOGTLN( "Next file will have ", _next_halos, " in ", _next_forests, "." );
   return first;
}

void
exporter::_open_file()
{
   // Close the file if it's open.
   if( _file.is_open() )
   {
      LOGTLN( "Closing file." );
      _file.close();
   }

   // Overwrite the file.
   string fn = _cur_filename();
   LOGTLN( "Opening next file for export: ", fn );
   _file.open( _cur_filename(), std::ios::out | std::ios::binary | std::ios::trunc );
   ASSERT( _file );
}

void
exporter::_write_header()
{
   _file.write( (const char*)&_next_forests, sizeof(unsigned) );
   _file.write( (const char*)&_next_halos, sizeof(unsigned) );
   auto it = _first;
   for( unsigned ii = 0; ii < _next_forests; ++ii )
   {
      ASSERT( it != _last );
      unsigned num_halos = _halo_cnts[it++];
      _file.write( (const char*)&num_halos, sizeof(unsigned) );
   }
}

string
exporter::_cur_filename()
{
   return string( "subfind." ) + to_string( _cur_file );
}

void
exporter::_calc_file_offset()
{
   _calc_num_local_files();
   if( OMP_TID == 0 )
      file_offsets[OMP_TID] = 0;
   if( OMP_TID < num_threads - 1 )
      file_offsets[OMP_TID + 1] = _num_local_files;
#pragma omp barrier
#pragma omp master
   for( unsigned ii = 1; ii < num_threads; ++ii )
      file_offsets[ii] += file_offsets[ii - 1];
#pragma omp barrier
   _cur_file = file_offsets[OMP_TID];

   // Need to reset some things.
   _next_forests = 0;
}

void
exporter::_calc_num_local_files()
{
   _num_local_files = 0;
   auto it = _first;
   long long num_halos = 0;
   while( it != _last )
   {
      it = _calc_next_count( it );
      ++_num_local_files;
      num_halos += _next_halos;
   }
}
