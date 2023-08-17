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
namespace ex
{

template <class T>
std::shared_ptr<Field<T> > laplacian
( 
    const Field<T>& fld 
)
{
    reuseTmp<T> lapt( fld.mesh() );
    std::shared_ptr<Field<T> > lap( lapt() );

    #pragma omp parallel for collapse(3)
    for( int i=settings::m()/2; i<fld.ni()-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<fld.nj()-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<fld.nk()-settings::m()/2; k++ )
            {
                for( int l=0; l<settings::m()/2; l++ )
                {
                    lap->operator()(i, j, k) += 
                        settings::coef()[l] / pow( (l+1) * fld.mesh().dy(), 2 ) * fld(i, j+1+l, k) //N
                       +settings::coef()[l] / pow( (l+1) * fld.mesh().dy(), 2 ) * fld(i, j-1-l, k) //S
                       +settings::coef()[l] / pow( (l+1) * fld.mesh().dx(), 2 ) * fld(i+1+l, j, k) //E
                       +settings::coef()[l] / pow( (l+1) * fld.mesh().dx(), 2 ) * fld(i-1-l, j, k) //W
                       +settings::coef()[l] / pow( (l+1) * fld.mesh().dz(), 2 ) * fld(i, j, k+1+l) //T
                       +settings::coef()[l] / pow( (l+1) * fld.mesh().dz(), 2 ) * fld(i, j, k-1-l) //B 
                       -(
                            2.0*settings::coef()[l] / pow( (l+1) * fld.mesh().dx(), 2 )
                           +2.0*settings::coef()[l] / pow( (l+1) * fld.mesh().dy(), 2 )
                           +2.0*settings::coef()[l] / pow( (l+1) * fld.mesh().dz(), 2 )
                        ) * fld(i, j, k); //P 
                }
            }
        }
    }
    
    return lap;
}




template <class T> 
std::shared_ptr<Field<T> > laplacian
(
    const Field<scalar>& nu,
    const Field<T>& fld
)
{
    reuseTmp<T> lapt( fld.mesh() );
    std::shared_ptr<Field<T> > lap( lapt() );


    for( int l=0; l<settings::m()/2; l++ )
    {   
        #pragma omp parallel for collapse(2)
        for( int i=settings::m()/2; i<fld.ni()-settings::m()/2; i++ )
        {   
            for( int j=settings::m()/2; j<fld.nj()-settings::m()/2; j++ )
            {  
                T* lptr = &lap->operator()(i,j,0); 
               
                for( int k=settings::m()/2; k<fld.nk()-settings::m()/2; k++ )
                {
                    lptr[k] += 
                    (
                        nu(i, j+l, k) * settings::coef()[l] / pow( (l+1) * fld.mesh().dy(), 2 ) * fld(i, j+1+l, k) //N
                       +nu(i, j-1-l, k) * settings::coef()[l] / pow( (l+1) * fld.mesh().dy(), 2 ) * fld(i, j-1-l, k) //S
                       +nu(i+l, j, k) * settings::coef()[l] / pow( (l+1) * fld.mesh().dx(), 2 ) * fld(i+1+l, j, k) //E
                       +nu(i-1-l, j, k) * settings::coef()[l] / pow( (l+1) * fld.mesh().dx(), 2 ) * fld(i-1-l, j, k) //W
                       +nu(i, j, k+l) * settings::coef()[l] / pow( (l+1) * fld.mesh().dz(), 2 ) * fld(i, j, k+1+l) / std::pow( fld.mesh().dzdxi(k), 2 )
                       +nu(i, j, k-1-l) * settings::coef()[l] / pow( (l+1) * fld.mesh().dz(), 2 ) * fld(i, j, k-1-l) / std::pow( fld.mesh().dzdxi(k), 2 )
                       +(
                           -(nu(i+l, j, k) + nu(i-1-l, j, k) ) * settings::coef()[l] / pow( (l+1) * fld.mesh().dx(), 2 )
                           -(nu(i, j+l, k) + nu(i, j-1-l, k) ) * settings::coef()[l] / pow( (l+1) * fld.mesh().dy(), 2 )
                           -(nu(i, j, k+l) + nu(i, j, k-1-l) ) * settings::coef()[l] / pow( (l+1) * fld.mesh().dz(), 2 ) / std::pow( fld.mesh().dzdxi(k), 2 )
                        ) * fld(i, j, k) //P
                     );
                }
            }
        }
    }


    for( int l=0; l<settings::m()/2; l++ )
    {
        for( int i=settings::m()/2; i<fld.ni()-settings::m()/2; i++ )
        {
            for( int j=settings::m()/2; j<fld.nj()-settings::m()/2; j++ )
            {
                T* lptr = &lap->operator()(i,j,0);

                for( int k=settings::m()/2; k<fld.nk()-settings::m()/2; k++ )
                {
                    lptr[k] -=
                    (
                        settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dz() * fld(i, j, k+1+l)
                       -settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dz() * fld(i, j, k-1-l)
                    ) * nu(i, j, k ) * fld.mesh().d2zdxi2(k) / std::pow( fld.mesh().dzdxi(k), 3 );
                }
            }
        }
    }

    return lap;
}


template <class T> 
std::shared_ptr<Field<T> > laplacian
(
    const scalar& nu,
    const Field<T>& fld
)
{
    reuseTmp<T> lapt( fld.mesh() );
    std::shared_ptr<Field<T> > lap( lapt() );

    for( int l=0; l<settings::m()/2; l++ )
    {
        #pragma omp parallel for collapse (2)
        for( int i=settings::m()/2; i<fld.ni()-settings::m()/2; i++ )
        {
            for( int j=settings::m()/2; j<fld.nj()-settings::m()/2; j++ )
            {
                T* lptr = &lap->operator()(i,j,0);
                
                for( int k=settings::m()/2; k<fld.nk()-settings::m()/2; k++ )
                {
                    lptr[k] += 
                    (
                        nu * settings::coef()[l] / pow( (l+1) * fld.mesh().dy(), 2 ) * fld(i, j+1+l, k) //N
                       +nu * settings::coef()[l] / pow( (l+1) * fld.mesh().dy(), 2 ) * fld(i, j-1-l, k) //S
                       +nu * settings::coef()[l] / pow( (l+1) * fld.mesh().dx(), 2 ) * fld(i+1+l, j, k) //E
                       +nu * settings::coef()[l] / pow( (l+1) * fld.mesh().dx(), 2 ) * fld(i-1-l, j, k) //W
                       +nu * settings::coef()[l] / pow( (l+1) * fld.mesh().dz(), 2 ) * fld(i, j, k+1+l) / std::pow( fld.mesh().dzdxi(k), 2 )
                       +nu * settings::coef()[l] / pow( (l+1) * fld.mesh().dz(), 2 ) * fld(i, j, k-1-l) / std::pow( fld.mesh().dzdxi(k), 2 )
                       +(
                           -(nu + nu ) * settings::coef()[l] / pow( (l+1) * fld.mesh().dx(), 2 )
                           -(nu + nu ) * settings::coef()[l] / pow( (l+1) * fld.mesh().dy(), 2 )
                           -(nu + nu ) * settings::coef()[l] / pow( (l+1) * fld.mesh().dz(), 2 ) / std::pow( fld.mesh().dzdxi(k), 2 )
                        ) * fld(i, j, k) //P
                    );
                }
            }
        }
    }

    for( int l=0; l<settings::m()/2; l++ )
    {
        for( int i=settings::m()/2; i<fld.ni()-settings::m()/2; i++ )
        {
            for( int j=settings::m()/2; j<fld.nj()-settings::m()/2; j++ )
            {
                T* lptr = &lap->operator()(i,j,0);

                for( int k=settings::m()/2; k<fld.nk()-settings::m()/2; k++ )
                {
                    lptr[k] -=
                    (
                        settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dz() * fld(i, j, k+1+l)
                       -settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dz() * fld(i, j, k-1-l)
                    ) * nu * fld.mesh().d2zdxi2(k) / std::pow( fld.mesh().dzdxi(k), 3 );
                }
            }
        }
    }

    return lap;
}

template<class T>
std::shared_ptr<Field<T> > laplacian
(
    const scalar& s,
    const std::shared_ptr<Field<T> > f
)
{
    return ex::laplacian( s, (*f) );
}

}

