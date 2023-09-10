#ifndef TOOLS_H
#define TOOLS_H

#include "Field/Field.h"
#include <memory>
#include <numbers>

namespace tools 
{

extern const scalar pi;
extern const scalar eps;

void CFL( const Field<vector>& U, const Mesh& mesh );


void tdma( std::vector<scalar>&, std::vector<scalar>&, std::vector<scalar>&, std::vector<scalar>&, std::vector<scalar>& );

template <class T>
T maximum( std::shared_ptr<Field<T> > f )
{
    T max=0.0;
    for( int i=settings::m()/2; i<f->ni()-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<f->nj()-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<f->nk()-settings::m()/2; k++ )
            {
                max = std::max(std::abs(f->operator()(i, j, k)),max);
            }
        }
    }

    reduce(max, maxOp<T>() ); 
    return max;
}
}

#endif
