#include "Types/vector/vector.h"

const int pTraits<vector>::nComp = 3;
const std::string pTraits<vector>::name = std::string( "vector" );

//constructors
vector::vector()
{
    v_[0] = v_[1] = v_[2] = 0.0;
}

vector::vector( scalar val )
{
    v_[0] = v_[1] = v_[2] = val;
}

vector::vector( scalar x, scalar y, scalar z )
{
    v_[0] = x;
    v_[1] = y;
    v_[2] = z;
}

//====
void vector::operator+=( vector r )
{
    this->x() += r.x();
    this->y() += r.y();
    this->z() += r.z();
}

vector& vector::operator-=( vector r )
{
    this->x() -= r.x();
    this->y() -= r.y();
    this->z() -= r.z();

    return *this;
}

vector& vector::operator*=( vector r )
{
    this->x() *= r.x();
    this->y() *= r.y();
    this->z() *= r.z();

    return *this;
}

vector& vector::operator/=( vector r )
{
    this->x() /= r.x();
    this->y() /= r.y();
    this->z() /= r.z();

    return *this;
}

std::ostream& operator<<(std::ostream& os, const vector& v)
{
    os << "(" << v.v_[0] << " " << v.v_[1] << " " << v.v_[2] << ")";
    return os;
}


//global operations
vector operator*
(
    scalar s,
    vector v
)
{
    return vector( s*v.x(), s*v.y(), s*v.z() );
}

vector operator*
(
    vector v,
    scalar s
)
{
    return vector( s*v.x(), s*v.y(), s*v.z() );
}

vector operator/
(
    vector v,
    scalar s
)
{
    return vector( v.x()/s, v.y()/s, v.z()/s );
}

vector operator-
(
    vector v
)
{
    return vector( -v.x(), -v.y(), -v.z() );
}

vector operator+
(
    const vector v1,
    const vector v2
)
{
    return vector( v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z() );
}

vector operator-
(
    const vector v1,
    const vector v2
)
{
    return vector( v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z() );
}

scalar operator&
(
    const vector v1,
    const vector v2
)
{
    return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}

vector operator&
(
    const vector v,
    const symmTensor t
)
{
    return vector
    (
        v.x() * t.xx() + v.y() * t.xy() + v.z() * t.xz(),
        v.x() * t.yx() + v.y() * t.yy() + v.z() * t.yz(),
        v.x() * t.zx() + v.y() * t.zy() + v.z() * t.zz()
    );
}


vector operator^
(
    const vector v1,
    const vector v2
)
{
    return vector
    (
        v1.y()*v2.z() - v1.z()*v2.y(),
        v1.z()*v2.x() - v1.x()*v2.z(),
        v1.x()*v2.y() - v1.y()*v2.x()
    );
}
