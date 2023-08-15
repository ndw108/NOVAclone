#ifndef TOOLS_H
#define TOOLS_H

#include "Field/Field.h"

namespace tools 
{

void CFL( const Field<vector>& U, const Mesh& mesh );

}


#endif
