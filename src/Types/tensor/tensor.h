#ifndef TENSOR_H
#define TENSOR_H

#include <ostream>
#include "Types/scalar/scalar.h"
#include "pTraits/pTraits.h"
#include "Types/symmTensor/symmTensor.h"
#include "Types/vector/vector.h"

#ifdef HAVE_MPI
#include "boost/mpi.hpp"
#endif

class tensor
{
#ifdef HAVE_MPI
    friend class boost::serialization::access ;

    template<class Archive> void serialize( Archive &ar, const unsigned int version )
    {   
        ar & v_ ;
    }
#endif

    scalar v_[9];

public:

    //constructors
    tensor();
    tensor( scalar );
    tensor
    ( 
        scalar, scalar, scalar, 
        scalar, scalar, scalar,
        scalar, scalar, scalar
    );

    tensor( const tensor& );
    tensor( const symmTensor& );

    //access
    scalar& xx()
    {
        return v_[0];
    }

    scalar& xy()
    {
        return v_[1];
    }

    scalar& xz()
    {
        return v_[2];
    }

    scalar& yx()
    {
        return v_[3];
    }

    scalar& yy()
    {
        return v_[4];
    }

    scalar& yz()
    {
        return v_[5];
    }

    scalar& zx()
    {
        return v_[6];
    }

    scalar& zy()
    {
        return v_[7];
    }


    scalar& zz()
    {
        return v_[8];
    }

    scalar xx() const
    {
        return v_[0];
    }

    scalar xy() const
    {
        return v_[1];
    }

    scalar xz() const
    {
        return v_[2];
    }

    scalar yx() const
    {
        return v_[3];
    }

    scalar yy() const
    {
        return v_[4];
    }

    scalar yz() const
    {
        return v_[5];
    }

    scalar zx() const
    {
        return v_[6];
    }

    scalar zy() const
    {
        return v_[7];
    }


    scalar zz() const
    {
        return v_[8];
    }



    scalar& operator[](const int dir)
    {
        return v_[dir];
    }

    scalar operator[](const int dir) const
    {
        return v_[dir];
    }


    //arithmatic
    tensor& operator+( const tensor& );
    tensor& operator+=( const tensor& );
    tensor& operator-( const tensor& );
    tensor& operator-=( const tensor& );
    tensor& operator*=( const tensor& );
    tensor& operator/=( const tensor& );

    friend std::ostream& operator<<(std::ostream&, const tensor&);
};

#ifdef HAVE_MPI
BOOST_IS_MPI_DATATYPE( tensor );
#endif

template <>
class pTraits<tensor>
{
public:
    static const int nComp;
    static const tensor I;
    static const std::string name;
};


inline scalar& component( tensor& v, const int dir )
{
    return v[dir];
}

tensor operator*
( 
    scalar,
    const tensor&
);

tensor operator*
( 
    const tensor&,
    scalar
);



tensor operator-
( 
    const tensor&
);


scalar tr
( 
    const tensor& 
);

tensor dev
( 
    const tensor& 
);

tensor dev2
( 
    const tensor& 
);

symmTensor twoSymm
(
    const tensor&
);

symmTensor symm
(
    const tensor&
);


tensor transpose
(
    const tensor&
);

vector operator&
(
    const vector& v,
    const tensor& t
);

scalar operator&&
(
    tensor&,
    tensor&
);

#endif
