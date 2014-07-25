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
  std::cout<<"Usage: "<<cmd<<" [options]\n"
           <<"\nOptions:\n"
           <<" -h, --help\n\tHelp! Prints this message.\n"
           <<" -v, --verbose\n\tVerbose output.\n"
	   <<" -c format, --convert format\n\tConvert image to another format. Options are vox, nrrd, inr.\n"
           <<" -s width, --slab width\n\tImage width.\n"
           <<" -t width, --throat width\n\tWidth of throat.\n"
           <<" -o filename, --output filename\n\tName of outfile.\n";
  return;
}

int parse_arguments(int argc, char **argv,
                    std::string &filename, bool &verbose, std::string &convert, int &slab_width, int &throat_width){

  // Set defaults
  filename = std::string("hourglass.vox");
  verbose = false;
  convert = std::string("vox");
  slab_width = 100;
  throat_width = 10;

  if(argc==1){
    usage(argv[0]);
    exit(0);
  }

  struct option longOptions[] = {
    {"help",    0,                 0, 'h'},
    {"verbose", 0,                 0, 'v'},
    {"convert", optional_argument, 0, 'c'},
    {"slab",    optional_argument, 0, 's'},
    {"throat",  optional_argument, 0, 't'},
    {"output",  optional_argument, 0, 'o'},
    {0, 0, 0, 0}
  };

  int optionIndex = 0;
  int verbosity = 0;
  int c;
  const char *shortopts = "hvc:s:t:o:";

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
    case 'c':
      convert = std::string(optarg);
      break;
    case 's':
      slab_width = atoi(optarg);
      break;    
    case 't':
      throat_width = atoi(optarg);
      break;
    case 'o':
      filename = std::string(optarg);
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

  return 0;
}

int main(int argc, char **argv){
  if(argc==1){
    usage(argv[0]);
    exit(-1);
  }
    
  std::string filename, convert;
  bool verbose, generate_mesh;
  int slab_width, throat_width;

  parse_arguments(argc, argv, filename, verbose, convert, slab_width, throat_width);

  CTImage image;
  if(verbose)
    image.verbose_on();
  
  image.create_hourglass(slab_width, throat_width);
  
  if(convert==std::string("vox")){
    if(verbose)
      std::cout<<"INFO: Write VOX file\n";

    image.write_vox(filename.c_str());
  }else if(convert==std::string("nrrd")){
    if(verbose)
      std::cout<<"INFO: Write NRRD file\n";

    image.write_nrrd(filename.c_str());
  }else if(convert==std::string("inr")){
    if(verbose)
      std::cout<<"INFO: Write INR file\n";

    image.write_inr(filename.c_str());
  }


  return 0;
}
