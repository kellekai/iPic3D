// File: PSKftiadaptor.h

#ifndef _PSK_FTI_ADAPTOR_H_
#define _PSK_FTI_ADAPTOR_H_

#include "PSKOutput.h"

#include <algorithm>

#include "fti.h"

namespace PSK {

  class FTIOutputAdaptor:public OutputAdaptor {

    //std::string _fti_file_name;

    static vector<void*> buffers;
    static int var_id;
    static vector<std::pair<FTIT_H5Group*, std::string>> groups;

    static std::string purify_object_name(const std::string & objname);
    //static void split_name(const std::string & name, std::vector < std::string > &elements);
    //static void free_buffer();
    static FTIT_H5Group* get_group_ptr(const std::string & tag); 
    static void init_groups();
    //void get_dataset_context(const std::string & name, std::vector < hid_t > &hid_array, std::string & dataset_name);

  public:
      static void init();
      FTIOutputAdaptor(void) {;
    }
    //void open(const std::string & name);
    //void open_append(const std::string & name);
    //void close(void);

    void write(const std::string & tag, int i_value);
    void write(const std::string & tag, long i_value);//
    void write(const std::string & tag, const Dimens dimens, const int *i_array);
    void write(const std::string & tag, const Dimens dimens, const long *i_array);

    void write(const std::string & tag, const Dimens dimens, const std::vector < int >&i_array);

    void write(const std::string & tag, const Dimens dimens, const std::vector < long >&i_array);

    void write(const std::string & objname, const Dimens dimens, const int ***i_array);


    // write float functions
    void write(const std::string & objname, float f);
    void write(const std::string & objname, const Dimens dimens, const float *f_array);
    void write(const std::string & objname, const Dimens dimens, const std::vector < float >&f_array);
    void write(const std::string & objname, const Dimens dimens, const float ***f_array);

    // write double functions
    void write(const std::string & objname, double d);
    void write(const std::string & objname, const Dimens dimens, const double *d_array);
    void write(const std::string & objname, const Dimens dimens, const std::vector < double >&d_array);
    void write(const std::string & objname, const Dimens dimens, double ***d_array);
    void write(const std::string & objname, const Dimens dimens, const int i, double ****d_array);

    void write(const std::string & objname, const Dimens dimens, double **d_array);

    void write(const std::string & objname, const Dimens dimens, const int i, double ***d_array);

  };

}                               // namespace

#endif

