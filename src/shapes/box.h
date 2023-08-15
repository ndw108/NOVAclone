#ifndef BOX_H
#define BOX_H

#include <sstream>
#include "shapes/shape.h"

class Box : public shape
{
    vector min_;
    vector max_;

    public:
    Box( vector min, vector max ) : min_(min), max_(max) {}

    Box( std::string line )
    {
        std::stringstream strProcessor(line);
        std::string buf;
        strProcessor>>buf;
        assert( buf == std::string( "Box:" ) );

        strProcessor>>buf>>min_.x()>>min_.y()>>min_.z()>>buf>>buf>>max_.x()>>max_.y()>>max_.z();  
    }    

    bool inside( vector xp )
    {
        return 
        (
            (xp.x() > min_.x() && xp.x() <= max_.x()) && 
            (xp.y() > min_.y() && xp.y() <= max_.y()) && 
            (xp.z() > min_.z() && xp.z() <= max_.z())
        );
    }

    vector& min()
    {
        return min_;
    }

    vector& max()
    {
        return max_;
    }

    bool intersect( Box& );
};

Box operator && (Box&,Box&);

#endif

