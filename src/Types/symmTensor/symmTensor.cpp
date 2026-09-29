#include "Types/symmTensor/symmTensor.h"

const int pTraits<symmTensor>::nComp = 6;
const symmTensor pTraits<symmTensor>::I = symmTensor
(
    1, 0, 0,
       1, 0,
          1
);
const std::string pTraits<symmTensor>::name = std::string( "symmTensor" );

//constructors
symmTensor::symmTensor()
{
    v_[0] = v_[1] = v_[2] = v_[3] = v_[4] = v_[5] = 0.0;
}

symmTensor::symmTensor( scalar val )
{
    v_[0] = v_[1] = v_[2] = v_[3] = v_[4] = v_[5] = val;
}

symmTensor::symmTensor( scalar xx, scalar xy, scalar xz, scalar yy, scalar yz, scalar zz )
{
    v_[0] = xx;
    v_[1] = xy;
    v_[2] = xz;
    v_[3] = yy;
    v_[4] = yz;
    v_[5] = zz;
}

symmTensor::symmTensor( const symmTensor& st )
:
symmTensor::symmTensor( st.xx(), st.xy(), st.xz(), st.yy(), st.yz(), st.zz() )
{
}

//====
symmTensor& symmTensor::operator+=( const symmTensor& r )
{
    this->xx() += r.xx();
    this->xy() += r.xy();
    this->xz() += r.xz();
    this->yy() += r.yy();
    this->yz() += r.yz();
    this->zz() += r.zz();

    return *this;
}

symmTensor& symmTensor::operator+( const symmTensor& r )
{
    this->xx() += r.xx();
    this->xy() += r.xy();
    this->xz() += r.xz();
    this->yy() += r.yy();
    this->yz() += r.yz();
    this->zz() += r.zz();

    return *this;
}

symmTensor& symmTensor::operator-=( const symmTensor& r )
{
    this->xx() -= r.xx();
    this->xy() -= r.xy();
    this->xz() -= r.xz();
    this->yy() -= r.yy();
    this->yz() -= r.yz();
    this->zz() -= r.zz();

    return *this;
}

symmTensor& symmTensor::operator-( const symmTensor& r )
{
    this->xx() -= r.xx();
    this->xy() -= r.xy();
    this->xz() -= r.xz();
    this->yy() -= r.yy();
    this->yz() -= r.yz();
    this->zz() -= r.zz();

    return *this;
}

symmTensor& symmTensor::operator*=( const symmTensor& r )
{
    this->xx() *= r.xx();
    this->xy() *= r.xy();
    this->xz() *= r.xz();
    this->yy() *= r.yy();
    this->yz() *= r.yz();
    this->zz() *= r.zz();

    return *this;
}

symmTensor& symmTensor::operator/=( const symmTensor& r )
{
    this->xx() /= r.xx();
    this->xy() /= r.xy();
    this->xz() /= r.xz();
    this->yy() /= r.yy();
    this->yz() /= r.yz();
    this->zz() /= r.zz();

    return *this;
}

std::ostream& operator<<(std::ostream& os, const symmTensor& v)
{
    os << "(" << v.v_[0] << " " << v.v_[1] << " " << v.v_[2] << " " << v.v_[3] << " " << v.v_[4] << " " << v.v_[5] << ")";
    return os;
}


//global operations
symmTensor operator*
(
    scalar s,
    const symmTensor& v
)
{
    return symmTensor
    ( 
        s*v.xx(), s*v.xy(), s*v.xz(),
                  s*v.yy(), s*v.yz(), 
                            s*v.zz() 
    );
}

symmTensor operator*
( 
    const symmTensor& st,
    const scalar s
)
{
    return s*st;
}

symmTensor operator-
(
    const symmTensor& v
)
{
    return -1*v; 
}


scalar tr
(
    const symmTensor& st
)
{
    return st.xx() + st.yy() + st.zz();
}


symmTensor dev
(
    const symmTensor& st
)
{
    return symmTensor( st ) -0.33333*tr(st)*pTraits<symmTensor>::I;
}

symmTensor dev2
(
    const symmTensor& st
)
{
    return symmTensor( st ) -0.66666*tr(st)*pTraits<symmTensor>::I;
}
