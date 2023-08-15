#include "Types/tensor/tensor.h"

const int pTraits<tensor>::nComp = 9;
const tensor pTraits<tensor>::I = tensor
(
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
);

const std::string pTraits<tensor>::name = std::string( "tensor" );

//constructors
tensor::tensor()
{
    v_[0] = v_[1] = v_[2] = v_[3] = v_[4] = v_[5] = v_[6] = v_[7] = v_[8] = 0.0;
}

tensor::tensor( scalar val )
{
    v_[0] = v_[1] = v_[2] = v_[3] = v_[4] = v_[5] = v_[6] = v_[7] = v_[8] = val;
}

tensor::tensor( scalar xx, scalar xy, scalar xz, scalar yx, scalar yy, scalar yz, scalar zx, scalar zy, scalar zz )
{
    v_[0] = xx;
    v_[1] = xy;
    v_[2] = xz;
    v_[3] = yx;
    v_[4] = yy;
    v_[5] = yz;
    v_[6] = zx;
    v_[7] = zy;
    v_[8] = zz; 
}

tensor::tensor( const tensor& t )
:
tensor::tensor( t.xx(), t.xy(), t.xz(), t.yx(), t.yy(), t.yz(), t.zx(), t.zy(), t.zz() )
{
}

tensor::tensor( const symmTensor& t )
:
tensor::tensor( t.xx(), t.xy(), t.xz(), t.yx(), t.yy(), t.yz(), t.zx(), t.zy(), t.zz() )
{
}

//====
tensor& tensor::operator+=( const tensor& r )
{
    for( int c=0; c<9; c++ )
    {
        this->operator[](c) += r[c];
    }

    return *this;
}

tensor& tensor::operator+( const tensor& r )
{
    for( int c=0; c<9; c++ )
    {
        this->operator[](c) += r[c];
    }

    return *this;
}

tensor& tensor::operator-=( const tensor& r )
{
    for( int c=0; c<9; c++ )
    {
        this->operator[](c) -= r[c];
    }

    return *this;
}

tensor& tensor::operator-( const tensor& r )
{
    for( int c=0; c<9; c++ )
    {
        this->operator[](c) -= r[c];
    }

    return *this;
}

tensor& tensor::operator*=( const tensor& r )
{
    for( int c=0; c<9; c++ )
    {
        this->operator[](c) *= r[c];
    }

    return *this;
}

tensor& tensor::operator/=( const tensor& r )
{
    for( int c=0; c<9; c++ )
    {
        this->operator[](c) /= r[c];
    }
    
    return *this;
}

std::ostream& operator<<(std::ostream& os, const tensor& v)
{
    os << "(";

    for( int c=0; c<9; c++ )
    {
        os << v.v_[c] << " ";    
    }
    os << ")";

    return os;
}


//global operations
tensor operator*
(
    scalar s,
    const tensor& v
)
{
    return tensor
    ( 
        s*v.xx(), s*v.xy(), s*v.xz(),
        s*v.yx(), s*v.yy(), s*v.yz(), 
        s*v.zx(), s*v.zy(), s*v.zz() 
    );
}

tensor operator*
(
    const tensor& v,
    scalar s
)
{
    return v*s;
}

tensor operator-
(
    const tensor& v
)
{
    return -1*v; 
}


scalar tr
(
    const tensor& t
)
{
    return t.xx() + t.yy() + t.zz();
}


tensor dev
(
    const tensor& t
)
{
    return tensor( t ) -0.33333*tr(t)*pTraits<tensor>::I;
}

tensor dev2
(
    const tensor& t
)
{
    return tensor( t ) -0.666666*tr(t)*pTraits<tensor>::I;
}

symmTensor twoSymm
(
    const tensor& t
)
{
    return symmTensor
    (
        2*t.xx(), t.xy() + t.yx(), t.xz() + t.zx(),
                  2*t.yy()       , t.yz() + t.zy(),
                                   2*t.zz()
    );
}

symmTensor symm
(
    const tensor& t
)
{
    return symmTensor
    (
        t.xx(), 0.5*(t.xy() + t.yx()), 0.5*(t.xz() + t.zx()),
                t.yy()               , 0.5*(t.yz() + t.zy()),
                                       t.zz()
    );
}

tensor transpose
(
    const tensor& t
)
{
    return tensor
    (
        t.xx(), t.yx(), t.zx(),
        t.xy(), t.yy(), t.zy(),
        t.xz(), t.yz(), t.zz()
    );
}

vector operator&
(
    const vector& v,
    const tensor& t
)
{
    vector res( 0,0,0 );

    res.x() = v.x()*t.xx() + v.y()*t.yx() + v.z()*t.zx();
    res.y() = v.x()*t.xy() + v.y()*t.yy() + v.z()*t.zy();
    res.z() = v.x()*t.xz() + v.y()*t.yz() + v.z()*t.zz();

    return res;
}


scalar operator&&
(
    tensor& t1, 
    tensor& t2
)
{
    scalar res=0;
    for( int i=0; i<9; i++ )
    {
        for( int j=0; j<9; j++ )
        {
            res+= component( t1, i ) * component( t2, j );
        }
    }   

    return res;
}
