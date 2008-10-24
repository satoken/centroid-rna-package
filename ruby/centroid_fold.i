%module CentroidFold
%{
#undef PACKAGE_VERSION
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_NAME
#undef PACKAGE_BUGREPORT
#undef VERSION
#include "../config.h"
#include "../src/centroid_fold.h"
%}
%include stl.i
%template() std::vector<std::string>;
%template() std::vector<float>;
%template() std::pair<std::string,float>;
//%apply std::string& OUTPUT { std::string& paren };
%ignore CentroidFold::decode_structure(float, std::string&) const;
%ignore CentroidFold::print(std::ostream&, const std::string&, const std::string&, const std::vector<double>&) const;
%ignore CentroidFold::print_posterior(std::ostream&, const std::string&, float) const;
%include "../src/centroid_fold.h"
