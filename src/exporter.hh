#ifndef bolshoi_exporter_hh
#define bolshoi_exporter_hh

#include <fstream>
#include <libhpc/libhpc.hh>
#include "types.hh"

class exporter
{
public:

   static const unsigned default_halos_per_file;

public:

   exporter( unsigned first_forest,
	     unsigned last_forest,
	     const hpc::vector<unsigned>& forest_halo_cnts,
	     unsigned halos_per_file=default_halos_per_file );

   ~exporter();

   void
   export_forest( const hpc::vector<sage_halo_type>& forest );

protected:

   unsigned
   _calc_next_count( unsigned first );

   void
   _open_file();

   void
   _write_header();

   hpc::string
   _cur_filename();

   void
   _calc_file_offset();

   void
   _calc_num_local_files();

protected:

   std::ofstream _file;
   unsigned _num_local_files;
   unsigned _cur_file;
   unsigned _hpf;
   unsigned _num_halos;
   unsigned _last_forests, _num_forests;
   unsigned _next_forests, _next_halos;
   unsigned _first, _last;
   const hpc::vector<unsigned>& _halo_cnts;
};

#endif
