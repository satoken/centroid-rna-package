%module CentroidFold
%{
#undef PACKAGE_VERSION
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_NAME
#undef PACKAGE_BUGREPORT
#undef VERSION
#include "../config.h"
#include "../src/centroid_fold_interface.h"
%}
%include stl.i
%template() std::vector<std::string>;
%template() std::vector<float>;
%template() std::pair<std::string,float>;
//%apply std::string& OUTPUT { std::string& paren };
%ignore CentroidFold::decode_structure(float, std::string&) const;
%include "../src/centroid_fold_interface.h"
