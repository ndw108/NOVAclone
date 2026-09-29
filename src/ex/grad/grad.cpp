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

#include "ex/grad/grad.h"

namespace ex
{

std::shared_ptr<Field<vector> > grad
(
    const Field<scalar>& fld
)
{
    reuseTmp<vector> ptrt( fld.mesh() );
    std::shared_ptr<Field<vector> > ptr( ptrt() );

    for( int l=0; l<settings::m()/2; l++ )
    {   
        for( int i=settings::m()/2; i<fld.ni()-settings::m()/2; i++ )
        {   
            for( int j=settings::m()/2; j<fld.nj()-settings::m()/2; j++ )
            {   
                vector* gptr = &ptr->operator()(i,j,0); 

                for( int k=settings::m()/2; k<fld.nk()-settings::m()/2; k++ )
                {
                    gptr[k] +=  
                    vector
                    (
                        settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dx() * fld(i+1+l, j, k)
                       -settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dx() * fld(i-1-l, j, k),

                        settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dy() * fld(i, j+1+l, k)
                       -settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dy() * fld(i, j-1-l, k),

                        settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dz() * fld(i, j, k+1+l) / fld.mesh().dzdxi(k)
                       -settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dz() * fld(i, j, k-1-l) / fld.mesh().dzdxi(k)
                    );
                }
            }
        }
    }

    return ptr;

}


std::shared_ptr<Field<tensor> > grad
(
    const Field<vector>& fld
)
{
    reuseTmp<tensor> ptrt( fld.mesh() );
    std::shared_ptr<Field<tensor> > ptr( ptrt() );

    #pragma omp parallel for collapse(3)
    for( int i=settings::m()/2; i<fld.ni()-settings::m()/2; i++ )
    {   
        for( int j=settings::m()/2; j<fld.nj()-settings::m()/2; j++ )
        {   
            for( int k=settings::m()/2; k<fld.nk()-settings::m()/2; k++ )
            {   
                ptr->operator()(i, j, k) = 0.0;

                for( int l=0; l<settings::m()/2; l++ )
                {   
                    scalar dx = settings::coef()[l] / (l+1) / fld.mesh().dx();
                    scalar dy = settings::coef()[l] / (l+1) / fld.mesh().dy();
                    scalar dz = settings::coef()[l] / (l+1) / fld.mesh().dz() / fld.mesh().dzdxi(k);

                    ptr->operator()(i, j, k) +=  
                    tensor
                    (
                        dx*(fld(i+l, j, k).x()-fld(i-1-l, j, k).x()), dy*(fld(i, j+l, k).x()-fld(i, j-1-l, k).x()), dz*(fld(i, j, k+l).x()-fld(i, j, k-1-l).x()),
                        dx*(fld(i+l, j, k).y()-fld(i-1-l, j, k).y()), dy*(fld(i, j+l, k).y()-fld(i, j-1-l, k).y()), dz*(fld(i, j, k+l).y()-fld(i, j, k-1-l).y()),
                        dx*(fld(i+l, j, k).z()-fld(i-1-l, j, k).z()), dy*(fld(i, j+l, k).z()-fld(i, j-1-l, k).z()), dz*(fld(i, j, k+l).z()-fld(i, j, k-1-l).z())
                    );
                }
            }
        }
    }

    return ptr;

}

std::shared_ptr<Field<vector> > grad
( 
    const std::shared_ptr<Field<scalar> > fld
)
{
    return ex::grad( (*fld) );
}

std::shared_ptr<Field<tensor> > grad
( 
    const std::shared_ptr<Field<vector> > fld
)
{
    return ex::grad( (*fld) );
}


}
