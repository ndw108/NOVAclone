/*  This file is part of NOVA.

    NOVA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    NOVA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with NOVA.  If not, see <http://www.gnu.org/licenses/>.

    Author: Alex Skillen. alex.skillen@manchester.ac.uk

*/
#include <math.h>
namespace im
{

template <class T>
std::shared_ptr<fdMatrix<T> > laplacian
( 
    const scalar mu,
    std::shared_ptr<Field<T> > fld 
)
{
    std::shared_ptr<fdMatrix<T> > matrix = std::make_shared<fdMatrix<T> >(fld);
    
    for( int i=settings::m()/2; i<fld->ni()-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<fld->nj()-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<fld->nk()-settings::m()/2; k++ )
            {
                matrix->al(0)(i, j, k) = mu / pow( fld->mesh().dy(), 2 ); //N
                matrix->al(1)(i, j, k) = mu / pow( fld->mesh().dy(), 2 ); //S
                matrix->al(2)(i, j, k) = mu / pow( fld->mesh().dx(), 2 ); //E
                matrix->al(3)(i, j, k) = mu / pow( fld->mesh().dx(), 2 ); //W
                matrix->al(4)(i, j, k) = mu / pow( fld->mesh().dz(), 2 ); //T
                matrix->al(5)(i, j, k) = mu / pow( fld->mesh().dz(), 2 ); //B

                matrix->ap()(i, j, k) = 
                -2.0 * mu / pow( fld->mesh().dx(), 2 )
                -2.0 * mu / pow( fld->mesh().dy(), 2 )
                -2.0 * mu / pow( fld->mesh().dz(), 2 ); //P 
            }
        }
    } 

    //apply deferred correction
    std::shared_ptr<Field<T> > hoa = ex::laplacian( mu, fld );
    matrix->deferredCorrection( hoa );

    return matrix;
}


template <class T>
std::shared_ptr<fdMatrix<T> > laplacian
( 
    std::shared_ptr<Field<T> > f
)
{
    return im::laplacian( scalar(1.0), f );
}


}
