#ifndef SCALAR_H
#define SCALAR_H

#include "pTraits/pTraits.h"

#ifdef HAVE_MPI
#include "boost/mpi.hpp"
#endif

typedef double scalar;

class vector;

template <>
class pTraits<scalar>
{
public:
    static const int nComp;
    static const std::string name;
};

template <>
class outerProductType<scalar, vector>
{
public:
    typedef vector type;
};



inline scalar& component( scalar& d, const int dir )
{
    return d;
}


#ifdef HAVE_MPI
BOOST_IS_MPI_DATATYPE( scalar );
#endif

#endif
