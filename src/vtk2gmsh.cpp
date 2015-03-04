/*  Copyright (C) 2010 Imperial College London and others.
 *
 *  Please see the AUTHORS file in the main source directory for a
 *  full list of copyright holders.
 *
 *  Gerard Gorman
 *  Applied Modelling and Computation Group
 *  Department of Earth Science and Engineering
 *  Imperial College London
 *
 *  g.gorman@imperial.ac.uk
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above
 *  copyright notice, this list of conditions and the following
 *  disclaimer in the documentation and/or other materials provided
 *  with the distribution.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 *  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 *  INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 *  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 *  ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
 *  TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
 *  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 *  SUCH DAMAGE.
 */

#include "writers.h"
#include "mesh_conversion.h"

#include <getopt.h>

void usage(char *cmd){
  std::cerr<<
    "Converts the VTK mesh output from Cleaver2 into a GMSH file. In order to "
    "avoid having disconnected domains we sweep across adjacent elements "
    "between two parallel sides of the domain and only keep mesh elements "
    "that are visited.\n"
    
    "Usage: "<<cmd<<" [options ...] [VTK mesh file]\n"
	     <<"\nOptions:\n"
           <<" -h, --help\n\tHelp! Prints this message.\n"
           <<" -n, --nhdr\n\nSpecify the NHDR file so that the meta-data can be read.\n"
           <<" -v, --verbose\n\tVerbose output.\n"
	   <<" -x, --x\n\tApply sweep align the x-axis (i.e. between the Y-Z parallel planes). This is the default.\n"
	   <<" -y, --y\n\tApply sweep align the y-axis (i.e. between the X-Z parallel planes).\n"
	   <<" -z, --z\n\tApply sweep align the z-axis (i.e. between the X-Y parallel planes).\n";
  return;
}

int parse_arguments(int argc, char **argv,
                    std::string &filename,
		    bool &verbose,
		    std::string &nhdr_filename,
		    int &axis){

  // Set defaults
  verbose = false;
  axis = 0;
  
  if(argc==1){
    usage(argv[0]);
    exit(0);
  }

  struct option longOptions[] = {
    {"help", 0, 0, 'h'},
    {"nhdr", optional_argument, 0, 'n'},
    {"verbose", 0, 0, 'v'},
    {"x",  0, 0, 'x'},
    {"y",  0, 0, 'y'},
    {"z",  0, 0, 'z'},
    {0, 0, 0, 0}
  };

  int optionIndex = 0;
  int verbosity = 0;
  int c;

  const char *shortopts = "hn:vxyz";

  // Set opterr to nonzero to make getopt print error messages
  opterr=1;
  while (true){
    c = getopt_long(argc, argv, shortopts, longOptions, &optionIndex);
    
    if (c == -1) break;
    
    switch (c){
    case 'h':
      usage(argv[0]);
      break;
    case 'n':
      nhdr_filename = std::string(optarg);
      break;
    case 'v':
      verbose = true;
      break;
    case 'x':
      axis = 0;
      break;
    case 'y':
      axis = 1;
      break;
    case 'z':
      axis = 2;
      break;
    case '?':
      // missing argument only returns ':' if the option string starts with ':'
      // but this seems to stop the printing of error messages by getopt?
      std::cerr<<"ERROR: unknown option or missing argument\n";
      usage(argv[0]);
      exit(-1);
    case ':':
      std::cerr<<"ERROR: missing argument\n";
      usage(argv[0]);
      exit(-1);
    default:
      // unexpected:
      std::cerr<<"ERROR: getopt returned unrecognized character code\n";
      exit(-1);
    }
  }

  filename = std::string(argv[argc-1]);
  std::ifstream infile;
  infile.open(filename.c_str());
  if (!infile.good()){
    std::cerr<<"ERROR: Cannot read file: "<<filename<<std::endl;
    usage(argv[0]);
    exit(0);
  }

}

int main(int argc, char **argv){
  std::string filename, nhdr_filename;
  bool verbose;
  int axis = 0;
  parse_arguments(argc, argv, filename, verbose, nhdr_filename, axis);

  std::string basename = filename.substr(0, filename.size()-4);
  
  if(verbose)
    std::cout<<"INFO: Reading "<<filename<<std::endl;

  std::vector<double> xyz;
  std::vector<int> tets;
  read_vtk_mesh_file(filename, nhdr_filename, xyz, tets);
  
  if(verbose)
    std::cout<<"INFO: Finished reading "<<filename<<std::endl;
  
  // Generate facets and trim disconnnected parts of the domain.
  std::vector<int> facets, facet_ids;
  if(verbose){
    std::cout<<"INFO: Create the active domain."<<std::endl;
    write_vtk_file(basename+"_original", xyz, tets, facets, facet_ids);
  }  

  create_domain(axis, xyz, tets, facets, facet_ids);
  
  if(tets.empty()){
    std::cerr<<"ERROR: There is no active region in the mesh. ";
    if(verbose)
      std::cerr<<"Check ";
    else
      std::cerr<<"Rerun the command with the -v option and check ";
    std::cerr<<"the file "<<basename+"_original.vtu to confirm there is indeed no connected region going between the two Y-Z planes.";

    return -1;
  }

  if(verbose) 
    std::cout<<"INFO: Active domain created."<<std::endl;
  
  if(verbose){
    std::cout<<"INFO: Writing out mesh."<<std::endl;
    write_vtk_file(basename, xyz, tets, facets, facet_ids);
  }

  write_gmsh_file(basename, xyz, tets, facets, facet_ids);
  if(verbose)
    std::cout<<"INFO: Finished."<<std::endl;

  return 0;
}
