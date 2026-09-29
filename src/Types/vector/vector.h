#ifndef VECTOR_H
#define VECTOR_H

#include <ostream>
#include "pTraits/pTraits.h"
#include "Types/scalar/scalar.h"
#include "Types/symmTensor/symmTensor.h"

#ifdef HAVE_MPI
#include "boost/mpi.hpp"
#endif

class vector
{
#ifdef HAVE_MPI
    friend class boost::serialization::access ;

    template<class Archive> void serialize( Archive &ar, const unsigned int version )
    {
        ar & v_ ;
    }
#endif

    scalar v_[3];

public:

    //constructors
    vector();
    vector( scalar );
    vector( scalar, scalar, scalar );

    //access
    inline scalar x() const
    {
        return v_[0];
    }

    inline scalar y() const
    {
        return v_[1];
    }

    inline scalar z() const
    {
        return v_[2];
    }

    inline scalar& x()
    {
        return v_[0];
    }

    inline scalar& y()
    {
        return v_[1];
    }

    inline scalar& z()
    {
        return v_[2];
    }

    inline scalar& operator[](const int dir)
    {
        return v_[dir];
    }

    //arithmatic 
    void operator+=( const vector );
    vector& operator-=( const vector );
    vector& operator*=( const vector );
    vector& operator/=( const vector );

    friend std::ostream& operator<<(std::ostream&, const vector&);
};

#ifdef HAVE_MPI
BOOST_IS_MPI_DATATYPE( vector );
#endif

template <>
class pTraits<vector>
{
public:
    static const int nComp;
    static const std::string name;
};

template <>
class outerProductType<vector, scalar>
{
public:
    typedef vector type;
};

template <>
class innerProductType<vector, symmTensor>
{
public:
    typedef vector type;
};

template <>
class innerProductType<vector, vector>
{
public:
    typedef scalar type;
};

inline scalar& component( vector& v, const int dir )
{
    return v[dir];
}

vector operator*
( 
    scalar,
    vector
);

vector operator*
( 
    vector,
    scalar
);

vector operator/
( 
    vector,
    scalar
);

vector operator-
( 
    vector
);

vector operator+
( 
    const vector,
    const vector
);

vector operator-
( 
    const vector,
    const vector
);

scalar operator&
( 
    const vector,
    const vector
);

vector operator&
(
    const vector,
    const symmTensor
);

vector operator^
(
    const vector,
    const vector
);

#endif
