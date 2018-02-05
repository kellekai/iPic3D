
#include "PSKftiadaptor.h"

using namespace PSK;

//FTIOutputAdaptor fti_agent;

// TODO hack
int ns = 2;

// definition of pointer array for protected variables
vector<void*> FTIOutputAdaptor::buffers;
int FTIOutputAdaptor::var_id;
vector<std::pair<FTIT_H5Group*, std::string>> FTIOutputAdaptor::groups;

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
    for(std::vector< std::pair<FTIT_H5Group*,std::string> >::iterator it=groups.begin(); it!=groups.end(); ++it) {
        if( tag.compare((*it).second) ) {
            return (*it).first;
        }
    }
}

void FTIOutputAdaptor::init() {

  var_id = 0;
#  ifdef FTI_HDF5
  init_groups();
#  endif // FTI_HDF5

}

void FTIOutputAdaptor::init_groups() {

  int i;
  FTIT_H5Group *group_ptr = (FTIT_H5Group*) malloc(sizeof(FTIT_H5Group));
  FTI_InitGroup( group_ptr, "topology", NULL );
  groups.push_back(std::pair<FTIT_H5Group*, std::string>(group_ptr, "/topology"));
 
  group_ptr = (FTIT_H5Group*) malloc(sizeof(FTIT_H5Group)); 
  FTI_InitGroup( group_ptr, "moments", NULL );
  groups.push_back(std::pair<FTIT_H5Group*, std::string>(group_ptr, "/moments"));
  
  group_ptr = (FTIT_H5Group*) malloc(sizeof(FTIT_H5Group)); 
  FTI_InitGroup( group_ptr, "rho", get_group_ptr("/moments") );
  groups.push_back(std::pair<FTIT_H5Group*, std::string>(group_ptr, "/rho"));
  for(i=0; i<ns; i++){
    group_ptr = (FTIT_H5Group*) malloc(sizeof(FTIT_H5Group)); 
    FTI_InitGroup( group_ptr, std::string("species_" + to_string(i)).c_str(), get_group_ptr("/rho") );
    groups.push_back(std::pair<FTIT_H5Group*, std::string>(group_ptr, "/rho/species_" + to_string(i)));
    group_ptr = (FTIT_H5Group*) malloc(sizeof(FTIT_H5Group)); 
    FTI_InitGroup( group_ptr, "rho", get_group_ptr("/rho/species_" + to_string(i)) );
    groups.push_back(std::pair<FTIT_H5Group*, std::string>(group_ptr, "/rho/species_" + to_string(i)+"/rho"));
  }
  
  group_ptr = (FTIT_H5Group*) malloc(sizeof(FTIT_H5Group)); 
  FTI_InitGroup( group_ptr, "particles", NULL );
  groups.push_back(std::pair<FTIT_H5Group*, std::string>(group_ptr, "/particles"));
  for(i=0; i<ns; i++) {
    group_ptr = (FTIT_H5Group*) malloc(sizeof(FTIT_H5Group)); 
    FTI_InitGroup( group_ptr, std::string("species_" + to_string(i)).c_str(), get_group_ptr("/particles"));
    groups.push_back(std::pair<FTIT_H5Group*, std::string>(group_ptr, "/particles/species_" + to_string(i)));
  }

  group_ptr = (FTIT_H5Group*) malloc(sizeof(FTIT_H5Group)); 
  FTI_InitGroup( group_ptr, "fields", NULL );
  groups.push_back(std::pair<FTIT_H5Group*, std::string>(group_ptr, "/fields"));
  
  group_ptr = (FTIT_H5Group*) malloc(sizeof(FTIT_H5Group)); 
  FTI_InitGroup( group_ptr, "Bx", get_group_ptr("/fields") );
  groups.push_back(std::pair<FTIT_H5Group*, std::string>(group_ptr, "/fields/Bx"));
  
  group_ptr = (FTIT_H5Group*) malloc(sizeof(FTIT_H5Group)); 
  FTI_InitGroup( group_ptr, "By", get_group_ptr("/fields") );
  groups.push_back(std::pair<FTIT_H5Group*, std::string>(group_ptr, "/fields/By"));
  
  group_ptr = (FTIT_H5Group*) malloc(sizeof(FTIT_H5Group)); 
  FTI_InitGroup( group_ptr, "Bz", get_group_ptr("/fields") );
  groups.push_back(std::pair<FTIT_H5Group*, std::string>(group_ptr, "/fields/Bz"));

  group_ptr = (FTIT_H5Group*) malloc(sizeof(FTIT_H5Group)); 
  FTI_InitGroup( group_ptr, "Ex", get_group_ptr("/fields") );
  groups.push_back(std::pair<FTIT_H5Group*, std::string>(group_ptr, "/fields/Ex"));
  
  group_ptr = (FTIT_H5Group*) malloc(sizeof(FTIT_H5Group)); 
  FTI_InitGroup( group_ptr, "Ey", get_group_ptr("/fields") );
  groups.push_back(std::pair<FTIT_H5Group*, std::string>(group_ptr, "/fields/Ey"));
  
  group_ptr = (FTIT_H5Group*) malloc(sizeof(FTIT_H5Group)); 
  FTI_InitGroup( group_ptr, "Ez", get_group_ptr("/fields") );
  groups.push_back(std::pair<FTIT_H5Group*, std::string>(group_ptr, "/fields/Ez"));

}

void FTIOutputAdaptor::write(const std::string & tag, int i_value) {
  
    std::string ptag = purify_object_name(tag);
    
    // get group name
    std::string group_name = ptag.substr(0, ptag.rfind("/"));
    
    // get name of data set
    std::string dataset_name = ptag.substr(ptag.rfind("/")+1);
    
    FTI_Protect(var_id, &i_value, 1, FTI_INTG);
    FTI_DefineDataset( var_id, 0, NULL, dataset_name.c_str(), get_group_ptr(group_name));
    var_id++;

}

void FTIOutputAdaptor::write(const std::string & tag, const Dimens dimens, const int *i_array){

    // get group name
    std::string group_name = tag.substr(0, tag.rfind("/"));
    
    // get name of data set
    std::string dataset_name = tag.substr(tag.rfind("/")+1);
    
    FTI_Protect(var_id, &i_array, dimens.nels(), FTI_INTG);
    int *dimensions = new int[dimens.size()];
    int i=0;
    for(; i<dimens.size(); i++) {
        dimensions[i] = dimens[i];
    }
    FTI_DefineDataset( var_id, dimens.size(), dimensions, dataset_name.c_str(), get_group_ptr(group_name));
    // TODO check if can be freed here !! delete[]dimensions;
    var_id++;

}
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
//void FTIOutputAdaptor::write(const std::string & tag, const Dimens dimens, const double *d_array) {
//  try {
//    if (dimens.size() == 0) {
//      PSK::OutputException e("Zero Dimens size", "FTIOutputAdaptor::write(double* array)");
//      throw e;
//    }
//
//    std::string ptag = purify_object_name(tag);
//
//    std::vector < hid_t > hid_array;
//    std::string dataset_name;
//
//    get_dataset_context(tag, hid_array, dataset_name);
//
//    hsize_t *hdf5dims = new hsize_t[dimens.size()];
//    for (int i = 0; i < dimens.size(); ++i)
//      hdf5dims[i] = dimens[i];
//
//    herr_t hdf5err = H5LTfind_dataset(hid_array[hid_array.size() - 1],
//                                      dataset_name.c_str());
//    if (hdf5err < 1) {
//      herr_t hdf5err = H5LTmake_dataset_double(hid_array[hid_array.size() - 1],
//                                               dataset_name.c_str(),
//                                               dimens.size(), hdf5dims, d_array);
//    }
//    else {
//      hid_t dataset_id = H5Dopen2(hid_array[hid_array.size() - 1], dataset_name.c_str(), H5P_DEFAULT); // HDF 1.8
//      hdf5err = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, d_array);
//      hdf5err = H5Dclose(dataset_id);
//    }
//
//
//    if (hdf5err < 0) {
//      PSK::OutputException e("make_dataset fails for " + tag, "FTIOutputAdaptor::write(double* array)");
//      throw e;
//    }
//
//    // close groups, if any, but don't try to close the file id at [0]
//    for (int i = hid_array.size() - 1; i > 0; --i)
//      hdf5err = H5Gclose(hid_array[i]);
//
//    delete [] hdf5dims;
//  } catch(PSK::Exception & e) {
//    e.push("In FTIOutputAdaptor::write(double* array)");
//    throw e;
//  }
//}
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
//    int nels = dimens.nels();
//    double *d_array_p = new double[nels];
//    const int di = dimens[0];
//    const int dj = dimens[1];
//    const int dk = dimens[2];
//    const int djk = dk * dj;
//    for (int i = 0; i < di; ++i)
//      for (int j = 0; j < dj; ++j)
//        for (int k = 0; k < dk; ++k) {
//
//          if (dk != 1)
//            d_array_p[i * djk + j * dk + k] = d_array[i + 1][j + 1][k + 1]; // I am not writing ghost cells
//          else if (dj != 1)
//            d_array_p[i * djk + j * dk] = d_array[i + 1][j + 1][0];
//          else {
//            d_array_p[i * djk + j * dk] = d_array[i + 1][0][0];
//
//          }
//        }
//
//    // add to buffer of protected variables
//    buffers.push_back((void*)d_array_p);
//    
//    // get name compontents of tag
//    std::string ptag = purify_object_name(objname);
//    std::vector < std::string > name_components;
//    split_name(ptag, name_components);
//
//    FTI_Protect(var_id, 
//    FTI_DefineDataset( var_id, dimens.size(), &dimens[0], ss.str().c_str(), &BxGroup); 
//    var_id++;
//
//    //FTI_Protect
//  }
//  catch(PSK::Exception & e) {
//    e.push("In FTIOutputAdaptor::write(double*** array)");
//    throw e;
//  }
//}
//
//void FTIOutputAdaptor::write(const std::string & objname, const Dimens dimens, const int ns, double ****d_array) {
//  if (dimens.size() != 3) {
//    PSK::OutputException e("Dimens size not 3 for object " + objname, "FTIOutputAdaptor::write(double**** array)");
//    throw e;
//  }
//
//  try {
//    int nels = dimens.nels();
//    double *d_array_p = new double[nels];
//    const int di = dimens[0];
//    const int dj = dimens[1];
//
//
//    const int dk = dimens[2];
//    const int djk = dk * dj;
//    for (int i = 0; i < di; ++i)
//      for (int j = 0; j < dj; ++j)
//        for (int k = 0; k < dk; ++k) {
//          if (dk != 1)
//            d_array_p[i * djk + j * dk + k] = d_array[ns][i + 1][j + 1][k + 1]; // I am not writing ghost cells
//          else if (dj != 1)
//            d_array_p[i * djk + j * dk] = d_array[ns][i + 1][j + 1][0];
//          else
//            d_array_p[i * djk + j * dk] = d_array[ns][i + 1][0][0];
//        }
//    write(objname, dimens, d_array_p);
//    delete[]d_array_p;
//  }
//  catch(PSK::Exception & e) {
//    e.push("In FTIOutputAdaptor::write(double**** array)");
//    throw e;
//  }
//}
//
