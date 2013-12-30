#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <limits>

#include <cassert>
#include <cstdlib>

#include <getopt.h>

#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>

void usage(char *cmd){
  std::cerr<<"\nUsage: "<<cmd<<" [options ...] [Tarantula mesh file]\n"
           <<"\nOptions:\n"
           <<" -h, --help\n\tHelp! Prints this message.\n"
           <<" -v, --verbose\n\tVerbose output.\n"
           <<" -b ID0,ID1, --boundary ID0, ID1\n\tBoundary ID's corresponding to the inflow and outflow.\n";
  return;
}

int parse_arguments(int argc, char **argv,
                    std::string &filename, bool &verbose, int *IDs){

  // Set defaults
  verbose = false;
  IDs[0]=-1;
  IDs[0]=-1;

  if(argc==1){
    usage(argv[0]);
    exit(0);
  }

  struct option longOptions[] = {
    {"help", 0, 0, 'h'},
    {"verbose", 0, 0, 'v'},
    {"boundary", optional_argument, 0, 'b'},
    
    {0, 0, 0, 0}
  };

  int optionIndex = 0;
  int verbosity = 0;
  int c;
  const char *shortopts = "hvb:";

  std::string boundary_str;
  size_t pos;
  std::stringstream id0, id1;

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
    case 'b':
      boundary_str = std::string(optarg);
      pos = boundary_str.find(",", 0);
      if(pos == std::string::npos){
        std::cerr<<"ERROR: boundary not specified correctly.\n";
        usage(argv[0]);
        exit(0);
      }
      id0<<boundary_str.substr(0, pos);
      if(!(id0>>IDs[0])) // Give the value to Result using the characters in the string
        IDs[0] = -1;
      id1<<boundary_str.substr(pos+1);
      if(!(id1>>IDs[1])) // Give the value to Result using the characters in the string
        IDs[1] = -1;

      // parse optarg
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

double volume(const double *x0, const double *x1, const double *x2, const double *x3){

    double x01 = (x0[0] - x1[0]);
    double x02 = (x0[0] - x2[0]);
    double x03 = (x0[0] - x3[0]);

    double y01 = (x0[1] - x1[1]);
    double y02 = (x0[1] - x2[1]);
    double y03 = (x0[1] - x3[1]);

    double z01 = (x0[2] - x1[2]);
    double z02 = (x0[2] - x2[2]);
    double z03 = (x0[2] - x3[2]);

    return (-x03*(z02*y01 - z01*y02) + x02*(z03*y01 - z01*y03) - x01*(z03*y02 - z02*y03))/6;
  }

int read_tarantula_mesh_file(std::string filename,
                             std::vector<double> &xyz,
                             std::vector<int> &tets, 
                             std::vector<int> &facets,
                             std::vector<int> &facet_ids){
  std::ifstream infile;
  infile.open(filename.c_str());
  
  // Read header
  std::string throwaway;
  std::getline(infile, throwaway); // line 1
  std::getline(infile, throwaway); // line 2
  int NNodes;
  infile>>NNodes;
  
  // Read vertices
  xyz.resize(NNodes*3);
  for(int i=0;i<NNodes;i++){
    infile>>xyz[i*3];
    infile>>xyz[i*3+1]; 
    infile>>xyz[i*3+2];
  }

  // throwaway trash
  std::getline(infile, throwaway); 
  std::getline(infile, throwaway); 
  std::getline(infile, throwaway);

  // Read elements
  int NTetra, nloc;
  infile>>NTetra;
  tets.resize(NTetra*4);
  for(int i=0;i<NTetra;i++){
    infile>>nloc;
    assert(nloc==4);
    infile>>tets[i*4];
    infile>>tets[i*4+1];
    infile>>tets[i*4+2];
    infile>>tets[i*4+3];
  }

  // Read facets
  while(infile.good()){
    // Stream through file until we find facet data.
    std::getline(infile, throwaway); 
    if(throwaway.substr(0, 4)!="Face")
      continue;

    // Get facet ID
    std::stringstream id_ss(throwaway.substr(4));
    int id;
    id_ss>>id;

    std::getline(infile, throwaway); 
    int nfacets;
    infile>>nfacets;
    nfacets/=2;

    for(int i=0;i<nfacets;i++){
      int eid, index;
      infile>>eid;
      infile>>index;
      assert(index<4);

      facet_ids.push_back(id);
      if(index==0){
        facets.push_back(tets[eid*4+0]); facets.push_back(tets[eid*4+1]); facets.push_back(tets[eid*4+2]);
      }else if(index==1){
        facets.push_back(tets[eid*4+1]); facets.push_back(tets[eid*4+0]); facets.push_back(tets[eid*4+3]);
      }else if(index==2){
        facets.push_back(tets[eid*4+2]); facets.push_back(tets[eid*4+1]); facets.push_back(tets[eid*4+3]);
      }else if(index==3){
        facets.push_back(tets[eid*4+0]); facets.push_back(tets[eid*4+2]); facets.push_back(tets[eid*4+3]);
      }
    }
  }

  infile.close();
}

int write_vtk_file(std::string filename,
                   std::vector<double> &xyz,
                   std::vector<int> &tets, 
                   std::vector<int> &facets,
                   std::vector<int> &facet_ids){
  
  // Write out points
  int NNodes = xyz.size()/3;
  vtkPoints *pts = vtkPoints::New();
  pts->SetNumberOfPoints(NNodes);
  for(int i=0;i<NNodes;i++)
    pts->SetPoint(i, &(xyz[i*3]));
  
  // Initalise the vtk mesh
  vtkUnstructuredGrid *ug_tets = vtkUnstructuredGrid::New();
  ug_tets->SetPoints(pts);

  int NTetra = tets.size()/4;
  for(int i=0;i<NTetra;i++){
    vtkIdList *idlist = vtkIdList::New();
    for(int j=0;j<4;j++)
      idlist->InsertNextId(tets[i*4+j]);
    ug_tets->InsertNextCell(10, idlist);
    idlist->Delete();
  }
  
  vtkXMLUnstructuredGridWriter *tet_writer = vtkXMLUnstructuredGridWriter::New();
  tet_writer->SetFileName(std::string(filename+".vtu").c_str());
  tet_writer->SetInput(ug_tets);
  tet_writer->Write();
  
  ug_tets->Delete();
  tet_writer->Delete();

  // Write out facets
  vtkUnstructuredGrid *ug_facets = vtkUnstructuredGrid::New();
  ug_facets->SetPoints(pts);
  int NFacets = facet_ids.size();
  for(int i=0;i<NFacets;i++){
    vtkIdList *idlist = vtkIdList::New();
    for(int j=0;j<3;j++){
      idlist->InsertNextId(facets[i*3+j]);
    }
    ug_facets->InsertNextCell(5, idlist);
    idlist->Delete();
  }

  vtkIntArray *vtk_facet_ids = vtkIntArray::New();
  vtk_facet_ids->SetNumberOfTuples(NFacets);
  vtk_facet_ids->SetNumberOfComponents(1);
  vtk_facet_ids->SetName("Facet IDs");
  for(int i=0;i<NFacets;i++){
    vtk_facet_ids->SetValue(i, facet_ids[i]);
  }
  ug_facets->GetCellData()->AddArray(vtk_facet_ids);
  vtk_facet_ids->Delete();
  
  vtkXMLUnstructuredGridWriter *tri_writer = vtkXMLUnstructuredGridWriter::New();
  tri_writer->SetFileName(std::string(filename+"_facets.vtu").c_str());
  tri_writer->SetInput(ug_facets);
  tri_writer->Write();

  pts->Delete();
  ug_facets->Delete();
  tri_writer->Delete();

  return 0;
}

int trim_channels(const int *id,
                  std::vector<double> &xyz,
                  std::vector<int> &tets, 
                  std::vector<int> &facets,
                  std::vector<int> &facet_ids){

  // Create node-element adjancy list - delete invested elements as we go.
  std::vector< std::set<int> > NEList(xyz.size()/3);
  int NTetra = tets.size()/4;
  int count_positive=0, count_negative=0;
  for(int i=0;i<NTetra;i++){
    if(tets[i*4]==-1)
      continue;
    
    double v = volume(&xyz[3*tets[i*4]], &xyz[3*tets[i*4+1]], &xyz[3*tets[i*4+2]], &xyz[3*tets[i*4+3]]);
    if(v<0){
      tets[i*4] = -1;
      count_negative++;
      continue;
    }else{
      count_positive++; 
    }
  
    for(int j=0;j<4;j++){
      NEList[tets[i*4+j]].insert(i);
    }
  }
  std::cout<<"Count of positive and negative volumes = "<<count_positive<<", "<<count_negative<<std::endl;

  // Create element-element adjancy list
  std::vector<int> EEList(NTetra*4, -1);
  for(int i=0;i<NTetra;i++){
    if(tets[i*4]==-1)
      continue;
  
    for(int j=0;j<4;j++){  
      std::set<int> edge_neighbours;
      set_intersection(NEList[tets[i*4+(j+1)%4]].begin(), NEList[tets[i*4+(j+1)%4]].end(),
                       NEList[tets[i*4+(j+2)%4]].begin(), NEList[tets[i*4+(j+2)%4]].end(),
                       inserter(edge_neighbours, edge_neighbours.begin()));
      
      std::set<int> neighbours;
      set_intersection(NEList[tets[i*4+(j+3)%4]].begin(), NEList[tets[i*4+(j+3)%4]].end(),
                       edge_neighbours.begin(), edge_neighbours.end(),
                       inserter(neighbours, neighbours.begin()));

      if(neighbours.size()==2){
        if(*neighbours.begin()==i)
          EEList[i*4+j] = *neighbours.rbegin();
        else
          EEList[i*4+j] = *neighbours.begin();
      }
#ifndef NDEBUG
      else{
	assert(neighbours.size()==1);
	assert(*neighbours.begin()==i);
      }
#endif
    }
  }

  // Create full facet ID list. Also, create the initial fronts for
  // the active region detection.
  std::map< std::set<int>, int> facet_id_lut;
  int NFacets = facet_ids.size();
  for(int i=0;i<NFacets;i++){
    std::set<int> facet;
    for(int j=0;j<3;j++){
      facet.insert(facets[i*3+j]);
    }
    assert(facet.size()==3);
    assert(facet_id_lut.find(facet)==facet_id_lut.end());
    facet_id_lut[facet] = facet_ids[i];
  }

  std::set<int> front0, front1;
  std::map< std::set<int>, int> facet_element_lut;
  for(int i=0;i<NTetra;i++){
    if(tets[i*4]==-1)
      continue;
    
    for(int j=0;j<4;j++){
      if(EEList[i*4+j]==-1){
        std::set<int> facet;
        for(int k=1;k<4;k++)
          facet.insert(tets[i*4+(j+k)%4]);
        assert(facet.size()==3);
	assert(facet_element_lut.find(facet)==facet_element_lut.end());
	facet_element_lut[facet] = i;

        std::map< std::set<int>, int>::iterator facet_id_pair = facet_id_lut.find(facet);
        if(facet_id_pair==facet_id_lut.end()){
          facet_ids.push_back(7);
	  
	  if(j==0){
            facets.push_back(tets[i*4+1]); facets.push_back(tets[i*4+3]); facets.push_back(tets[i*4+2]); 
	  }else if(j==1){
            facets.push_back(tets[i*4]); facets.push_back(tets[i*4+3]); facets.push_back(tets[i*4+2]);
	  }else if(j==2){
            facets.push_back(tets[i*4]); facets.push_back(tets[i*4+3]); facets.push_back(tets[i*4+1]);
	  }else if(j==3){
            facets.push_back(tets[i*4]); facets.push_back(tets[i*4+2]); facets.push_back(tets[i*4+1]);
	  }
        }else{
          if(facet_id_pair->second==id[0])
            front0.insert(i);
          else if(facet_id_pair->second==id[1])
            front1.insert(i);
        }
      }
    }
  }
  
  // Advance front0
  std::vector<int> label(NTetra, 0);
  while(!front0.empty()){
    // Get the next unprocessed element in the set.
    int seed = *front0.begin();
    front0.erase(front0.begin());
    if(label[seed]==1)
      continue;
    label[seed] = 1;
    
    for(int i=0;i<4;i++){
      int eid = EEList[seed*4+i];
      if(eid!=-1 && label[eid]!=1){
        front0.insert(eid);
      }
    }
  }

  // Advance back sweep using front1.
  while(!front1.empty()){
    // Get the next unprocessed element in the set.
    int seed = *front1.begin();
    front1.erase(front1.begin());
    if(label[seed]!=1) // ie was either never of interest or has been processed in the backsweep.
      continue;
    label[seed] = 2;
    
    for(int i=0;i<4;i++){
      int eid = EEList[seed*4+i];
      if(eid!=-1 && label[eid]==1){
        front1.insert(eid);
      }
    }
  }

  // Find active vertex set and create renumbering.
  std::map<int, int> renumbering;
  for(int i=0;i<NTetra;i++){
    if(tets[i*4]==-1)
      continue;
    
    if(label[i]==2){
      for(int j=0;j<4;j++)
        renumbering.insert(std::pair<int, int>(tets[i*4+j], -1));
    }
  }

  // Create new compressed mesh.
  std::vector<double> xyz_new;
  std::vector<int> tets_new;
  std::vector<int> facets_new;
  std::vector<int> facet_ids_new;
  int cnt=0;
  for(std::map<int, int>::iterator it=renumbering.begin();it!=renumbering.end();++it){
    it->second = cnt++;
    
    xyz_new.push_back(xyz[(it->first)*3]);
    xyz_new.push_back(xyz[(it->first)*3+1]);
    xyz_new.push_back(xyz[(it->first)*3+2]);
  }
  for(int i=0;i<NTetra;i++){
    if(tets[i*4]==-1)
      continue;
    
    if(label[i]==2){
      for(int j=0;j<4;j++){
        tets_new.push_back(renumbering[tets[i*4+j]]);
      }
    }
  }
  NFacets = facet_ids.size();
  for(int i=0;i<NFacets;i++){
    std::set<int> facet;
    for(int j=0;j<3;j++)
      facet.insert(facets[i*3+j]);

    // Check if this is an orphaned facet from previous purge.
    if(facet_element_lut.find(facet)==facet_element_lut.end())
      continue;

    // Check if this is a newly orphaned facet.
    if(label[facet_element_lut[facet]]!=2){
      continue;
    }
    
    for(int j=0;j<3;j++){
      std::map<int, int>::iterator it=renumbering.find(facets[i*3+j]);
      assert(it!=renumbering.end());
      facets_new.push_back(it->second);
    }
    
    facet_ids_new.push_back(facet_ids[i]);
  }
  xyz.swap(xyz_new);
  tets.swap(tets_new);
  facets.swap(facets_new);
  facet_ids.swap(facet_ids_new);
}

int purge_locked(size_t NNodes, std::vector<int> &tets){
  // Create node-element adjancy list.
  std::vector< std::set<int> > NEList(NNodes);
  size_t NTetra = tets.size()/4;
  for(size_t i=0;i<NTetra;i++){
    for(size_t j=0;j<4;j++){
      NEList[tets[i*4+j]].insert(i);
    }
  }
  
  // Create element-element adjancy list.
  std::vector<int> EEList(NTetra*4, -1);
  for(size_t i=0;i<NTetra;i++){
    for(size_t j=0;j<4;j++){
      
      std::set<int> edge_neighbours;
      set_intersection(NEList[tets[i*4+(j+1)%4]].begin(), NEList[tets[i*4+(j+1)%4]].end(),
                       NEList[tets[i*4+(j+2)%4]].begin(), NEList[tets[i*4+(j+2)%4]].end(),
                       inserter(edge_neighbours, edge_neighbours.begin()));
      
      std::set<int> neighbours;
      set_intersection(NEList[tets[i*4+(j+3)%4]].begin(), NEList[tets[i*4+(j+3)%4]].end(),
                       edge_neighbours.begin(), edge_neighbours.end(),
                       inserter(neighbours, neighbours.begin()));

      if(neighbours.size()==2){
        if(*neighbours.begin()==i)
          EEList[i*4+j] = *neighbours.rbegin();
        else
          EEList[i*4+j] = *neighbours.begin();
      }
#ifndef NDEBUG
      else{
	assert(neighbours.size()==1);
	assert(*neighbours.begin()==i);
      }
#endif
    }
  }

  // Create list of boundary nodes.
  std::set<int> boundary_nodes;
  for(size_t i=0;i<NTetra;i++){
    for(size_t j=0;j<4;j++){
      if(EEList[i*4+j]==-1){
        for(size_t k=1;k<4;k++)
          boundary_nodes.insert(tets[i*4+(j+k)%4]);
      }
    }
  }

  // Mask out elements that are locked (all nodes on the boundary).
  for(size_t i=0;i<NTetra;i++){
    bool locked = true;
    for(size_t j=0;j<4;j++){
      if(boundary_nodes.find(tets[i*4+j])==boundary_nodes.end()){
	locked = false;
	break;
      }
    }
    if(locked){
      tets[i*4] = -1;
    }
  }
}

int write_triangle_file(std::string basename,
                        std::vector<double> &xyz,
                        std::vector<int> &tets, 
                        std::vector<int> &facets,
                        std::vector<int> &facet_ids){
  std::string filename_node = basename+".node";
  std::string filename_face = basename+".face";
  std::string filename_ele = basename+".ele";
  
  int NNodes = xyz.size()/3;
  int NTetra = tets.size()/4;
  int NFacets = facet_ids.size();
  assert(NFacets==facets.size()/3);

  ofstream nodefile;
  nodefile.open(std::string(basename+".node").c_str());
  nodefile<<NNodes<<" "<<3<<" "<<0<<" "<<0<<std::endl;
  nodefile<<std::setprecision(std::numeric_limits<double>::digits10+1);
  
  for(int i=0;i<NNodes;i++){
    nodefile<<i+1<<" "<<xyz[i*3]<<" "<<xyz[i*3+1]<<" "<<xyz[i*3+2]<<std::endl;
  }

  ofstream elefile;
  elefile.open(std::string(basename+".ele").c_str());
  elefile<<NTetra<<" "<<4<<" "<<1<<std::endl;

  for(int i=0;i<NTetra;i++){
    elefile<<i+1<<" "<<tets[i*4]+1<<" "<<tets[i*4+1]+1<<" "<<tets[i*4+2]+1<<" "<<tets[i*4+3]+1<<" 1"<<std::endl;
  }

  ofstream facefile;
  facefile.open(std::string(basename+".face").c_str());
  facefile<<NFacets<<" "<<1<<std::endl;
  for(int i=0;i<NFacets;i++){
    facefile<<i+1<<" "<<facets[i*3]+1<<" "<<facets[i*3+1]+1<<" "<<facets[i*3+2]+1<<" "<<facet_ids[i]<<std::endl;
  }
  
  return 0;
}

int write_gmsh_file(std::string basename,
		    std::vector<double> &xyz,
		    std::vector<int> &tets, 
		    std::vector<int> &facets,
		    std::vector<int> &facet_ids){
  
  int NNodes = xyz.size()/3;
  int NTetra = tets.size()/4;
  int NFacets = facet_ids.size();
  assert(NFacets==facets.size()/3);

  ofstream file;
  file.open(std::string(basename+".msh").c_str());
  file<<"$MeshFormat"<<std::endl
      <<"2.2 0 8"<<std::endl
      <<"$EndMeshFormat"<<std::endl
      <<"$Nodes"<<std::endl
      <<NNodes<<std::endl;
  file<<std::setprecision(std::numeric_limits<double>::digits10+1);
  for(size_t i=0;i<NNodes;i++){
    file<<i+1<<" "<<xyz[i*3]<<" "<<xyz[i*3+1]<<" "<<xyz[i*3+2]<<std::endl;
  }
  file<<"$EndNodes"<<std::endl
      <<"$Elements"<<std::endl
      <<NTetra+NFacets<<std::endl;
  for(size_t i=0;i<NTetra;i++){
    file<<i+1<<" 4 1 1 "<<tets[i*4]+1<<" "<<tets[i*4+1]+1<<" "<<tets[i*4+2]+1<<" "<<tets[i*4+3]+1<<std::endl;
  }
  for(size_t i=0;i<NFacets;i++){
    file<<i+NTetra+1<<" 2 1 "<<facet_ids[i]<<" "<<facets[i*3]+1<<" "<<facets[i*3+1]+1<<" "<<facets[i*3+2]+1<<std::endl;
  }
  file<<"$EndElements"<<std::endl;
  file.close();

  return 0;
}

int main(int argc, char **argv){
  std::string filename;
  bool verbose;
  int IDs[2];
  parse_arguments(argc, argv, filename, verbose, IDs);

  std::vector<double> xyz;
  std::vector<int> tets, facets, facet_ids;
  std::string basename = filename.substr(0, filename.size()-4);

  if(verbose) std::cout<<"INFO: Reading "<<filename<<std::endl;
  read_tarantula_mesh_file(filename, xyz, tets, facets, facet_ids);
  if(verbose) std::cout<<"INFO: Finished reading "<<filename<<std::endl;

  write_vtk_file(basename+"_original", xyz, tets, facets, facet_ids);

  // Characterise boundaries
  std::map<int, std::vector<double> > mean_boundary_coordinate;
  std::map<int, int> boundary_count;
  for(size_t i=0;i<facet_ids.size();i++){
    std::vector<double> coord(3, 0.0);
    for(size_t j=0;j<3;j++){
      for(size_t k=0;k<3;k++){
        coord[k]+=xyz[facets[i*3+j]*3+k];
      }
    }
    if(mean_boundary_coordinate.count(facet_ids[i])){
      for(size_t k=0;k<3;k++){
        mean_boundary_coordinate[facet_ids[i]][k]+=coord[k];
      }
      boundary_count[facet_ids[i]]++;
    }else{
      mean_boundary_coordinate[facet_ids[i]]=coord;
      boundary_count[facet_ids[i]] = 1;
    }
  }
  std::vector<double> bbox(6);
  std::vector<int> sorted_ids(6);
  {
    std::map<int, std::vector<double> >::iterator it=mean_boundary_coordinate.begin();
    for(size_t i=0;i<3;i++){
      it->second[i]/=(3*boundary_count[it->first]);

      bbox[i*2] = it->second[i];
      bbox[i*2+1] = it->second[i];
  
      sorted_ids[i*2] = it->first;
      sorted_ids[i*2+1] = it->first;
    }
    ++it;
    for(;it!=mean_boundary_coordinate.end();++it){
      for(size_t i=0;i<3;i++){
        it->second[i]/=(3*boundary_count[it->first]);

        if(it->second[i]<bbox[i*2]){
          bbox[i*2] = it->second[i];
          sorted_ids[i*2] = it->first;
        }

        if(it->second[i]>bbox[i*2+1]){
          bbox[i*2+1] = it->second[i];
          sorted_ids[i*2+1] = it->first;
        }
      }
    }
  }

  // 
  std::vector<double> orig_xyz = xyz;
  std::vector<int> orig_tets = tets, orig_facets = facets, orig_facet_ids = facet_ids;

  for(size_t i=0;i<3;i++){
    if(i!=0){
      xyz = orig_xyz;
      tets = orig_tets;
      facets = orig_facets;
      facet_ids = orig_facet_ids;
    }

    if(verbose) std::cout<<"INFO: Trimming inactive regions."<<filename<<std::endl;
    trim_channels(&(sorted_ids[i*2]), xyz, tets, facets, facet_ids);
    if(verbose) std::cout<<"INFO: Finished trimming."<<filename<<std::endl;

    std::ostringstream filename;
    filename<<basename<<"_axis_"<<i<<"_ids_"<<sorted_ids[i*2]<<"_"<<sorted_ids[i*2+1];
    write_vtk_file(filename.str(), xyz, tets, facets, facet_ids);
    write_gmsh_file(filename.str(), xyz, tets, facets, facet_ids);
  }

  return 0;
}
