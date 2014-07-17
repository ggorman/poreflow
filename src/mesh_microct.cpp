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

#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>

#include "CTImage.h"

void usage(char *cmd){
  std::cout<<"Usage: "<<cmd<<" CT-image [options]\n"
           <<"\nOptions:\n"
           <<" -h, --help\n\tHelp! Prints this message.\n"
           <<" -v, --verbose\n\tVerbose output.\n"
           <<" -s width, --slab width\n\tExtract a square block of size 'width' from the data.\n";
  return;
}

int parse_arguments(int argc, char **argv,
                    std::string &filename, bool &verbose, int &slab_width){

  // Set defaults
  verbose = false;
  slab_width = -1;

  if(argc==1){
    usage(argv[0]);
    exit(0);
  }

  struct option longOptions[] = {
    {"help",    0,                 0, 'h'},
    {"verbose", 0,                 0, 'v'},
    {"slab",    optional_argument, 0, 's'},
    {0, 0, 0, 0}
  };

  int optionIndex = 0;
  int verbosity = 0;
  int c;
  const char *shortopts = "hvs:";

  // Set opterr to nonzero to make getopt print error messages
  opterr=1;
  while (true){
    c = getopt_long(argc, argv, shortopts, longOptions, &optionIndex);
    
    if (c == -1) break;
    
    switch (c){
    case 'h':
      usage(argv[0]);
      break;
    case 'v':
      verbose = true;
      break;
    case 's':
      slab_width = atoi(optarg);
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

  return 0;
}

int main(int argc, char **argv){
  if(argc==1){
    usage(argv[0]);
    exit(-1);
  }
    
  std::string filename;
  bool verbose;
  int slab_width;

  parse_arguments(argc, argv, filename, verbose, slab_width);

  CTImage image;
  if(verbose)
    image.verbose_on();
  
  image.read_raw_ese_image(filename.c_str(), slab_width);
  
  if(verbose)
    std::cout<<"INFO: Generate mesh using CGAL.\n";

  image.mesh();
    
  if(verbose)
    std::cout<<"INFO: Trim disconnected regions.\n";

  image.trim_channels(1, 2);
    
  if(verbose){
    std::cout<<"INFO: Write out VTK file.\n";

    image.write_vtu();
  }

  if(verbose)
    std::cout<<"INFO: Write out GMSH file.\n";

  image.write_gmsh();

  return 0;
}
