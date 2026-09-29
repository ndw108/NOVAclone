#include "shapes/box.h"

bool Box::intersect( Box& b )
{
    return
    (
        b.max().x() > this->min().x() &&
        b.min().x() < this->max().x() &&
        b.max().y() > this->min().y() &&
        b.min().y() < this->max().y() &&
        b.max().z() > this->min().z() &&
        b.min().z() < this->max().z()
    );
}


Box operator &&
(
    Box& b1,
    Box& b2
)
{
    if( !b1.intersect( b2 ) )
    {
        return Box( vector(0,0,0), vector( 0,0,0 ) );
    }

    vector min
    ( 
        std::max( b1.min().x(), b2.min().x() ), 
        std::max( b1.min().y(), b2.min().y() ), 
        std::max( b1.min().z(), b2.min().z() )
    );

    vector max
    ( 
        std::min( b1.max().x(), b2.max().x() ), 
        std::min( b1.max().y(), b2.max().y() ), 
        std::min( b1.max().z(), b2.max().z() )
    );


    return Box( min, max );
}
