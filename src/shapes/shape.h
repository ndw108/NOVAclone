#ifndef SHAPE_H
#define SHAPE_H

#include "Types/vector/vector.h"
#include <cassert>

class shape
{
    public:
    virtual bool inside(vector) =0;
    virtual vector& min() =0;
    virtual vector& max() =0;
};

#endif
