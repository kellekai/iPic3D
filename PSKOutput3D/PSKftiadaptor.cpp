
#include "PSKftiadaptor.h"

using namespace PSK;

//FTIOutputAdaptor fti_agent;

// TODO hack
int FTIOutputAdaptor::ns = 2;

// definition of pointer array for protected variables
std::vector<void*> FTIOutputAdaptor::buffers;
int FTIOutputAdaptor::var_id;
std::vector<std::pair<FTIT_H5Group, std::string>> FTIOutputAdaptor::groups;

/* 
 * 
 * abort if zero length
 * 
 * add leading "/" if doesn't start with "/"
 * 
 */
std::string FTIOutputAdaptor::purify_object_name(const std::string & objname) {
  if (objname.length() == 0)
    throw PSK::OutputException("Zero length tag name", "FTIOutputAdaptor::purify_object_name()");

  return objname[0] != '/' ? "/" + objname : objname;

}

FTIT_H5Group* FTIOutputAdaptor::get_group_ptr(const std::string & tag) {
    for(std::vector< std::pair<FTIT_H5Group,std::string> >::iterator it=groups.begin(); it!=groups.end(); ++it) {
        if( tag.compare((*it).second) == 0 ) {
            return &((*it).first);
        }
    }
    return NULL;
}

void FTIOutputAdaptor::init() {

  FTI_Init("config.fti", MPI_COMM_WORLD);
   
  var_id = 0;
#  ifdef FTI_HDF5
  init_groups();
#  endif // FTI_HDF5

}

void FTIOutputAdaptor::init_groups() {
  
  int i;
  FTIT_H5Group H5Group;
  std::pair<FTIT_H5Group, std::string> H5Pair;
  
  FTI_InitGroup( &H5Group, "topology", NULL );
  H5Pair.first = H5Group; 
  H5Pair.second = "/topology";
  groups.push_back(H5Pair);

  FTI_InitGroup( &H5Group, "moments", NULL );
  H5Pair.first = H5Group; 
  H5Pair.second = "/moments";
  groups.push_back(H5Pair);

  FTI_InitGroup( &H5Group, "rho", get_group_ptr("/moments") );
  H5Pair.first = H5Group; 
  H5Pair.second = "/moments/rho";
  groups.push_back(H5Pair);

  for(i=0; i<ns; i++){
    FTI_InitGroup( &H5Group, std::string("species_" + to_string(i)).c_str(), get_group_ptr("/moments") );
    H5Pair.first = H5Group; 
    H5Pair.second = "/moments/species_" + to_string(i);
    groups.push_back(H5Pair);
    FTI_InitGroup( &H5Group, "rho", get_group_ptr("/moments/species_" + to_string(i)) );
    H5Pair.first = H5Group; 
    H5Pair.second = "/moments/species_" + to_string(i) + "/rho";
    groups.push_back(H5Pair);
  }
  
  FTI_InitGroup( &H5Group, "particles", NULL );
  H5Pair.first = H5Group; 
  H5Pair.second = "/particles";
  groups.push_back(H5Pair);
  for(i=0; i<ns; i++) {
    FTI_InitGroup( &H5Group, std::string("species_" + to_string(i)).c_str(), get_group_ptr("/particles"));
    H5Pair.first = H5Group; 
    H5Pair.second = "/particles/species_" + to_string(i);
    groups.push_back(H5Pair);
    FTI_InitGroup( &H5Group, "q", get_group_ptr("/particles/species_" + to_string(i)));
    H5Pair.first = H5Group; 
    H5Pair.second = "/particles/species_" + to_string(i) + "/q";
    groups.push_back(H5Pair);
    FTI_InitGroup( &H5Group, "u", get_group_ptr("/particles/species_" + to_string(i)));
    H5Pair.first = H5Group; 
    H5Pair.second = "/particles/species_" + to_string(i) + "/u";
    groups.push_back(H5Pair);
    FTI_InitGroup( &H5Group, "v", get_group_ptr("/particles/species_" + to_string(i)));
    H5Pair.first = H5Group; 
    H5Pair.second = "/particles/species_" + to_string(i) + "/v";
    groups.push_back(H5Pair);
    FTI_InitGroup( &H5Group, "w", get_group_ptr("/particles/species_" + to_string(i)));
    H5Pair.first = H5Group; 
    H5Pair.second = "/particles/species_" + to_string(i) + "/w";
    groups.push_back(H5Pair);
    FTI_InitGroup( &H5Group, "x", get_group_ptr("/particles/species_" + to_string(i)));
    H5Pair.first = H5Group; 
    H5Pair.second = "/particles/species_" + to_string(i) + "/x";
    groups.push_back(H5Pair);
    FTI_InitGroup( &H5Group, "y", get_group_ptr("/particles/species_" + to_string(i)));
    H5Pair.first = H5Group; 
    H5Pair.second = "/particles/species_" + to_string(i) + "/y";
    groups.push_back(H5Pair);
    FTI_InitGroup( &H5Group, "z", get_group_ptr("/particles/species_" + to_string(i)));
    H5Pair.first = H5Group; 
    H5Pair.second = "/particles/species_" + to_string(i) + "/z";
    groups.push_back(H5Pair);
  }

  FTI_InitGroup( &H5Group, "fields", NULL );
  H5Pair.first = H5Group; 
  H5Pair.second = "/fields";
  groups.push_back(H5Pair);
  
  FTI_InitGroup( &H5Group, "Bx", get_group_ptr("/fields") );
  H5Pair.first = H5Group; 
  H5Pair.second = "/fields/Bx";
  groups.push_back(H5Pair);
  
  FTI_InitGroup( &H5Group, "By", get_group_ptr("/fields") );
  H5Pair.first = H5Group; 
  H5Pair.second = "/fields/By";
  groups.push_back(H5Pair);
  
  FTI_InitGroup( &H5Group, "Bz", get_group_ptr("/fields") );
  H5Pair.first = H5Group; 
  H5Pair.second = "/fields/Bz";
  groups.push_back(H5Pair);

  FTI_InitGroup( &H5Group, "Ex", get_group_ptr("/fields") );
  H5Pair.first = H5Group; 
  H5Pair.second = "/fields/Ex";
  groups.push_back(H5Pair);
  
  FTI_InitGroup( &H5Group, "Ey", get_group_ptr("/fields") );
  H5Pair.first = H5Group; 
  H5Pair.second = "/fields/Ey";
  groups.push_back(H5Pair);
  
  FTI_InitGroup( &H5Group, "Ez", get_group_ptr("/fields") );
  H5Pair.first = H5Group; 
  H5Pair.second = "/fields/Ez";
  groups.push_back(H5Pair);
  

}

void FTIOutputAdaptor::write(const std::string & tag, int i_value) {
  
    std::string ptag = purify_object_name(tag);
    
    // get group name
    std::string group_name = ptag.substr(0, ptag.rfind("/"));
    
    // get name of data set
    std::string dataset_name = ptag.substr(ptag.rfind("/")+1);
   
    // allocate buffer and add do FTI buffer list
    void * buffer = malloc(sizeof(int));
    memcpy(buffer, &i_value, sizeof(int));
    buffers.push_back(buffer);

    FTI_Protect(var_id, buffer, 1, FTI_INTG);
    
    FTI_DefineDataset( var_id, 0, NULL, dataset_name.c_str(), get_group_ptr(group_name));
    var_id++;

}

void FTIOutputAdaptor::write(const std::string & tag, const Dimens dimens, const int *i_array){

    // add '/' as first character if missing 
    std::string ptag = purify_object_name(tag);
    
    // get group name
    std::string group_name = ptag.substr(0, tag.rfind("/"));
    
    // get name of data set
    std::string dataset_name = ptag.substr(tag.rfind("/")+1);
   
    // allocate buffer and add do FTI buffer list
    void * buffer = malloc(sizeof(int)*dimens.nels());
    memcpy(buffer, i_array, sizeof(int)*dimens.nels());
    buffers.push_back(buffer);

    FTI_Protect(var_id, buffer, dimens.nels(), FTI_INTG);
    
    int *dimensions = new int[dimens.size()];
    int i=0;
    for(; i<dimens.size(); i++) {
        dimensions[i] = dimens[i];
    }

    FTI_DefineDataset( var_id, dimens.size(), dimensions, dataset_name.c_str(), get_group_ptr(group_name));
    // TODO check if can be freed here !! delete[]dimensions;
    var_id++;

}

void FTIOutputAdaptor::write(const std::string & tag, const Dimens dimens, double *d_array) {

    // add '/' as first character if missing 
    std::string ptag = purify_object_name(tag);
    
    // get group name
    std::string group_name = ptag.substr(0, tag.rfind("/"));
    
    // get name of data set
    std::string dataset_name = ptag.substr(tag.rfind("/")+1);
    
    FTI_Protect(var_id, d_array, dimens.nels(), FTI_DBLE);
    
    int *dimensions = new int[dimens.size()];
    int i=0;
    for(; i<dimens.size(); i++) {
        dimensions[i] = dimens[i];
    }

    FTI_DefineDataset( var_id, dimens.size(), dimensions, dataset_name.c_str(), get_group_ptr(group_name));
    // TODO check if can be freed here !! delete[]dimensions;
    var_id++;

}

void FTIOutputAdaptor::write(const std::string & tag, const Dimens dimens, double ***d_array){

    int nels = dimens.nels();
    double *d_array_p = new double[nels];
    const int di = dimens[0];
    const int dj = dimens[1];
    const int dk = dimens[2];
    const int djk = dk * dj;
    for (int i = 0; i < di; ++i)
      for (int j = 0; j < dj; ++j)
        for (int k = 0; k < dk; ++k) {

          if (dk != 1)
            d_array_p[i * djk + j * dk + k] = d_array[i + 1][j + 1][k + 1]; // I am not writing ghost cells
          else if (dj != 1)
            d_array_p[i * djk + j * dk] = d_array[i + 1][j + 1][0];
          else {
            d_array_p[i * djk + j * dk] = d_array[i + 1][0][0];

          }
        }
    
    // add '/' as first character if missing 
    std::string ptag = purify_object_name(tag);
    
    // get group name
    std::string group_name = ptag.substr(0, tag.rfind("/"));
    
    // get name of data set
    std::string dataset_name = ptag.substr(tag.rfind("/")+1);

    // add to buffer of protected variables
    buffers.push_back((void*)d_array_p);
    
    int *dimensions = new int[dimens.size()];
    for(int i=0; i<dimens.size(); i++) {
        dimensions[i] = dimens[i];
    }

    FTI_Protect( var_id, (void*)d_array_p, nels, FTI_DBLE );  
    FTI_DefineDataset( var_id, dimens.size(), dimensions, dataset_name.c_str(), get_group_ptr(group_name));
    var_id++;

}

void FTIOutputAdaptor::write(const std::string & tag, const Dimens dimens, const int ns, double ****d_array) {
    int nels = dimens.nels();
    double *d_array_p = new double[nels];
    const int di = dimens[0];
    const int dj = dimens[1];


    const int dk = dimens[2];
    const int djk = dk * dj;
    for (int i = 0; i < di; ++i)
      for (int j = 0; j < dj; ++j)
        for (int k = 0; k < dk; ++k) {
          if (dk != 1)
            d_array_p[i * djk + j * dk + k] = d_array[ns][i + 1][j + 1][k + 1]; // I am not writing ghost cells
          else if (dj != 1)
            d_array_p[i * djk + j * dk] = d_array[ns][i + 1][j + 1][0];
          else
            d_array_p[i * djk + j * dk] = d_array[ns][i + 1][0][0];
        }
    
    // add '/' as first character if missing 
    std::string ptag = purify_object_name(tag);
    
    // get group name
    std::string group_name = ptag.substr(0, tag.rfind("/"));
    
    // get name of data set
    std::string dataset_name = ptag.substr(tag.rfind("/")+1);

    // add to buffer of protected variables
    buffers.push_back((void*)d_array_p);
    
    int *dimensions = new int[dimens.size()];
    for(int i=0; i<dimens.size(); i++) {
        dimensions[i] = dimens[i];
    }

    FTI_Protect( var_id, (void*)d_array_p, nels, FTI_DBLE );  
    FTI_DefineDataset( var_id, dimens.size(), dimensions, dataset_name.c_str(), get_group_ptr(group_name));
    var_id++;
}
















    void FTIOutputAdaptor::write(const std::string & tag, long i_value){}//
    void FTIOutputAdaptor::write(const std::string & tag, const Dimens dimens, const long *i_array){}

    void FTIOutputAdaptor::write(const std::string & tag, const Dimens dimens, const std::vector < int >&i_array){}

    void FTIOutputAdaptor::write(const std::string & tag, const Dimens dimens, const std::vector < long >&i_array){}

    void FTIOutputAdaptor::write(const std::string & objname, const Dimens dimens, const int ***i_array){}


    // w FTIOutputAdaptor::ite float functions
    void FTIOutputAdaptor::write(const std::string & objname, float f){}
    void FTIOutputAdaptor::write(const std::string & objname, const Dimens dimens, const float *f_array){}
    void FTIOutputAdaptor::write(const std::string & objname, const Dimens dimens, const std::vector < float >&f_array){}
    void FTIOutputAdaptor::write(const std::string & objname, const Dimens dimens, const float ***f_array){}

    // w FTIOutputAdaptor::ite double functions
    void FTIOutputAdaptor::write(const std::string & objname, double d){}
    //void FTIOutputAdaptor::write(const std::string & objname, const Dimens dimens, const double *d_array){}
    void FTIOutputAdaptor::write(const std::string & objname, const Dimens dimens, const std::vector < double >&d_array){}
    //void FTIOutputAdaptor::write(const std::string & objname, const Dimens dimens, double ***d_array){}
    //void FTIOutputAdaptor::write(const std::string & objname, const Dimens dimens, const int i, double ****d_array){}

    void FTIOutputAdaptor::write(const std::string & objname, const Dimens dimens, double **d_array){}

    //void FTIOutputAdaptor::write(const std::string & objname, const Dimens dimens, const int i, double ***d_array){}


















//
//void FTIOutputAdaptor::write(const std::string & tag, const Dimens dimens, const std::vector < int >&i_array) {
//  try {
//    int n = dimens.nels();
//    int *i_array_p = new int[n];
//    for (int i = 0; i < n; ++i)
//      i_array_p[i] = i_array[i];
//    write(tag, dimens, i_array_p);
//    delete[]i_array_p;
//  } catch(PSK::Exception & e) {
//    e.push("In FTIOutputAdaptor::write(vector<int> array)");
//    throw e;
//  }
//}
//
//void FTIOutputAdaptor::write(const std::string & tag, const Dimens dimens, const std::vector < long >&i_array) {
//  try {
//    int n = dimens.nels();
//    long *i_array_p = new long[n];
//    for (int i = 0; i < n; ++i)
//      i_array_p[i] = i_array[i];
//    write(tag, dimens, i_array_p);
//    delete[]i_array_p;
//  } catch(PSK::Exception & e) {
//    e.push("In FTIOutputAdaptor::write(vector<long> array)");
//    throw e;
//  }
//}
//
//
//void FTIOutputAdaptor::write(const std::string & tag, const Dimens dimens, const std::vector < double >&d_array) {
//  try {
//    int n = dimens.nels();
//    double *d_array_p = new double[n];
//    for (int i = 0; i < n; ++i)
//      d_array_p[i] = d_array[i];
//    write(tag, dimens, d_array_p);
//    delete[]d_array_p;
//  } catch(PSK::Exception & e) {
//    e.push("In FTIOutputAdaptor::write(vector<double> array)");
//    throw e;
//  }
//}
//
//void FTIOutputAdaptor::Init_Checkpoint(const std::string & tag, const Dimens dimens, void *d_array) {
//  try {
//    if (dimens.size() == 0) {
//      PSK::OutputException e("Zero Dimens size", "HDF5OutputAdaptor::write(double* array)");
//      throw e;
//    }
//   
//    buffers.push_back();
//    std::string ptag = purify_object_name(tag);
//
//    std::vector < hid_t > hid_array;
//    std::string dataset_name;
//
//  } catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(double* array)");
//    throw e;
//  }
//}
//void FTIOutputAdaptor::write(const std::string & objname, const Dimens dimens, double ***d_array) {
//  if (dimens.size() != 3) {
//    PSK::OutputException e("Dimens size not 3 for object " + objname, "FTIOutputAdaptor::FTI_Init_Checkpoint(double*** array)");
//    throw e;
//  }
//
//  try {
//
//    //FTI_Protect
//  }
//  catch(PSK::Exception & e) {
//    e.push("In FTIOutputAdaptor::write(double*** array)");
//    throw e;
//  }
//}
//
//
