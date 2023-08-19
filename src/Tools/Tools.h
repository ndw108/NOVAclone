#ifndef TOOLS_H
#define TOOLS_H

#include "Field/Field.h"
#include <memory>

namespace tools 
{

void CFL( const Field<vector>& U, const Mesh& mesh );


void tdma( std::vector<scalar>&, std::vector<scalar>&, std::vector<scalar>&, std::vector<scalar>&, std::vector<scalar>& );

}


#endif
