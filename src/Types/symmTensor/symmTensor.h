#ifndef SYMMTENSOR_H
#define SYMMTENSOR_H

#include <ostream>
#include "pTraits/pTraits.h"
#include "Types/scalar/scalar.h"

class symmTensor
{
#ifdef HAVE_MPI
    friend class boost::serialization::access ;

    template<class Archive> void serialize( Archive &ar, const unsigned int version )
    {
        ar & v_ ;
    }
#endif

    scalar v_[6];

public:

    //constructors
    symmTensor();
    symmTensor( scalar );
    symmTensor
    ( 
        scalar, scalar, scalar, 
                scalar, scalar,
                        scalar
    );
    symmTensor( const symmTensor& );

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

    scalar& yy()
    {
        return v_[3];
    }

    scalar& yz()
    {
        return v_[4];
    }

    scalar& zz()
    {
        return v_[5];
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
        return v_[1];
    }

    scalar yy() const
    {
        return v_[3];
    }

    scalar yz() const
    {
        return v_[4];
    }

    scalar zx() const
    {
        return v_[2];
    }

    scalar zy() const
    {
        return v_[4];
    }

    scalar zz() const
    {
        return v_[5];
    }



    scalar& operator[](const int dir)
    {
        return v_[dir];
    }

    //arithmatic
    symmTensor& operator+( const symmTensor& );
    symmTensor& operator+=( const symmTensor& );
    symmTensor& operator-( const symmTensor& );
    symmTensor& operator-=( const symmTensor& );
    symmTensor& operator*=( const symmTensor& );
    symmTensor& operator/=( const symmTensor& );

    friend std::ostream& operator<<(std::ostream&, const symmTensor&);
};

#ifdef HAVE_MPI
BOOST_IS_MPI_DATATYPE( symmTensor );
#endif



template <>
class pTraits<symmTensor>
{
public:
    static const int nComp;
    static const symmTensor I;
    static const std::string name;
};


inline scalar& component( symmTensor& v, const int dir )
{
    return v[dir];
}

symmTensor operator*
( 
    scalar,
    const symmTensor&
);

symmTensor operator*
( 
    const symmTensor&,
    const scalar
);

symmTensor operator-
( 
    const symmTensor&
);


scalar tr
( 
    const symmTensor& 
);

symmTensor dev
( 
    const symmTensor& 
);

symmTensor dev2
( 
    const symmTensor& 
);



#endif
