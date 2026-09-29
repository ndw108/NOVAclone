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
namespace ex
{

std::shared_ptr<Field<scalar> > div
( 
    const Field<vector>& fld 
)
{
    reuseTmp<scalar> divfldt( fld.mesh() );
    std::shared_ptr<Field<scalar> > divfld( divfldt() );

    std::unique_ptr<scalar[]> Up(new scalar[fld.nk()]);
    std::unique_ptr<scalar[]> Um(new scalar[fld.nk()]);
    std::unique_ptr<scalar[]> Vp(new scalar[fld.nk()]);
    std::unique_ptr<scalar[]> Vm(new scalar[fld.nk()]);
    std::unique_ptr<scalar[]> Wp(new scalar[fld.nk()]);
    std::unique_ptr<scalar[]> Wm(new scalar[fld.nk()]);

    for( int l=0; l<settings::m()/2; l++ )
    {
        for( int i=settings::m()/2; i<fld.ni()-settings::m()/2; i++ )
        {
            for( int j=settings::m()/2; j<fld.nj()-settings::m()/2; j++ )
            {
                           
                scalar* dptr = &divfld->operator()(i,j,0); 

                for( int k=settings::m()/2; k<fld.nk()-settings::m()/2; k++ )
                { 
                    Vp[k] = fld(i, j+1+l, k).y(); 
                    Vm[k] = fld(i, j-1-l, k).y();
                    Up[k] = fld(i+1+l, j, k).x();
                    Um[k] = fld(i-1-l, j, k).x();
                    Wp[k] = fld(i, j, k+1+l).z();
                    Wm[k] = fld(i, j, k-1-l).z();
                }

                for( int k=settings::m()/2; k<fld.nk()-settings::m()/2; k++ )
                {
                    dptr[k] += 
                    (
                        settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dy() * Vp[k] //N
                       -settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dy() * Vm[k] //S
                       +settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dx() * Up[k] //E
                    );
                }

                for( int k=settings::m()/2; k<fld.nk()-settings::m()/2; k++ )
                {
                    dptr[k] += 
                    (
                       -settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dx() * Um[k] //W
                       +settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dz() * Wp[k] / fld.mesh().dzdxi(k) //T
                       -settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dz() * Wm[k] / fld.mesh().dzdxi(k) //B
                    );
                }
            }
        }
    } 

    return divfld;

}


template<class T>
std::shared_ptr<Field<T> > div
( 
    const Field<vector>& fld1, 
    const Field<T>& fld2
)
{
    reuseTmp<T> divfldt( fld2.mesh() );
    std::shared_ptr<Field<T> > divfld( divfldt() );

    std::unique_ptr<scalar[]> Up(new scalar[fld2.nk()]);
    std::unique_ptr<scalar[]> Um(new scalar[fld2.nk()]);
    std::unique_ptr<scalar[]> Vp(new scalar[fld2.nk()]);
    std::unique_ptr<scalar[]> Vm(new scalar[fld2.nk()]);
    std::unique_ptr<scalar[]> Wp(new scalar[fld2.nk()]);
    std::unique_ptr<scalar[]> Wm(new scalar[fld2.nk()]);

    for( int l=0; l<settings::m()/2; l++ )
    {
        for( int i=settings::m()/2; i<fld1.ni()-settings::m()/2; i++ )
        {
            for( int j=settings::m()/2; j<fld1.nj()-settings::m()/2; j++ )
            {
            
                T* dptr = &divfld->operator()(i,j,0);

                for( int k=settings::m()/2; k<fld1.nk()-settings::m()/2; k++ )
                { 
                    Vp[k] = fld1(i, j+1+l, k).y(); 
                    Vm[k] = fld1(i, j-1-l, k).y();
                    Up[k] = fld1(i+1+l, j, k).x();
                    Um[k] = fld1(i-1-l, j, k).x();
                    Wp[k] = fld1(i, j, k+1+l).z();
                    Wm[k] = fld1(i, j, k-1-l).z();
                }
 
                for( int k=settings::m()/2; k<fld1.nk()-settings::m()/2; k++ )
                { 
                    dptr[k] +=  
                    (
                        0.5*settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dy() * fld2(i, j+1+l, k) * Vp[k] //N 
                       -0.5*settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dy() * fld2(i, j-1-l, k) * Vm[k] //S 
                       +0.5*settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dx() * fld2(i+1+l, j, k) * Up[k] //E
                    );
                }
                for( int k=settings::m()/2; k<fld1.nk()-settings::m()/2; k++ )
                { 
                    dptr[k] +=  
                    (

                       -0.5*settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dx() * fld2(i-1-l, j, k) * Um[k] //W 
                       +0.5*settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dz() * fld2(i, j, k+1+l) * Wp[k] / fld1.mesh().dzdxi(k) //T 
                       -0.5*settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dz() * fld2(i, j, k-1-l) * Wm[k] / fld1.mesh().dzdxi(k)//B 
                    );
                }
    

                for( int k=settings::m()/2; k<fld1.nk()-settings::m()/2; k++ )
                { 
                    dptr[k] +=  
                    (
                        0.5*fld1(i, j, k).y() * settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dy() * fld2(i, j+1+l, k) //N
                       -0.5*fld1(i, j, k).y() * settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dy() * fld2(i, j-1-l, k) //S
                       +0.5*fld1(i, j, k).x() * settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dx() * fld2(i+1+l, j, k) //E
                       -0.5*fld1(i, j, k).x() * settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dx() * fld2(i-1-l, j, k) //W
                       +0.5*fld1(i, j, k).z() * settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dz() * fld2(i, j, k+1+l) / fld1.mesh().dzdxi(k)
                       -0.5*fld1(i, j, k).z() * settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dz() * fld2(i, j, k-1-l) / fld1.mesh().dzdxi(k) 
                    );
                }

                for( int k=settings::m()/2; k<fld1.nk()-settings::m()/2; k++ )
                {
                    dptr[k] +=  
                    ( 
                        0.5*fld2(i, j, k) * (settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dy() * Vp[k]) //N
                       -0.5*fld2(i, j, k) * (settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dy() * Vm[k]) //S
                       +0.5*fld2(i, j, k) * (settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dx() * Up[k]) //E
                       -0.5*fld2(i, j, k) * (settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dx() * Um[k]) //W
                       +0.5*fld2(i, j, k) * (settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dz() * Wp[k]) / fld1.mesh().dzdxi(k)
                       -0.5*fld2(i, j, k) * (settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dz() * Wm[k]) / fld1.mesh().dzdxi(k)
                    );
                }
               
            }
        }
    }

    return divfld;
}


std::shared_ptr<Field<vector> > div
( 
    const Field<scalar>& fld1, 
    const std::shared_ptr<Field<tensor> > fld2 
)
{
    reuseTmp<vector> divfldt( fld1.mesh() );
    std::shared_ptr<Field<vector> > divfld( divfldt() );

    #pragma omp parallel for collapse(3)
    for( int i=settings::m()/2; i<fld1.ni()-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<fld1.nj()-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<fld1.nk()-settings::m()/2; k++ )
            {
                divfld->operator()(i, j, k)=0.0;

                for( int l=0; l<settings::m()/2; l++ )
                {
                    for( int c=0; c<3; c++ )
                    {
                        divfld->operator()(i, j, k)[c] +=  settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dy() * fld1(i, j+1+l, k) * fld2->operator()(i, j+1+l, k)[c*3+1]; //N
                        divfld->operator()(i, j, k)[c] += -settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dy() * fld1(i, j-1-l, k) * fld2->operator()(i, j-1-l, k)[c*3+1]; //S
                        divfld->operator()(i, j, k)[c] +=  settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dx() * fld1(i+1+l, j, k) * fld2->operator()(i+1+l, j, k)[c*3]; //E
                        divfld->operator()(i, j, k)[c] += -settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dx() * fld1(i-1-l, j, k) * fld2->operator()(i-1-l, j, k)[c*3]; //W
                        divfld->operator()(i, j, k)[c] +=  settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dz() * fld1(i, j, k+1+l) * fld2->operator()(i, j, k+1+l)[c*3+2]; //T
                        divfld->operator()(i, j, k)[c] += -settings::coef()[l] / 2.0 / (l+1) / fld1.mesh().dz() * fld1(i, j, k-1-l) * fld2->operator()(i, j, k-1-l)[c*3+2]; //B
                    }
                }
            }
        }
    }

    return divfld;
}



std::shared_ptr<Field<vector> > div
( 
    const Field<symmTensor>& fld 
)
{
    reuseTmp<vector> divfldt( fld.mesh() );
    std::shared_ptr<Field<vector> > divfld( divfldt() );

    #pragma omp parallel for collapse(3)
    for( int i=settings::m()/2; i<fld.ni()-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<fld.nj()-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<fld.nk()-settings::m()/2; k++ )
            {
                divfld->operator()(i, j, k)=0.0;

                for( int l=0; l<settings::m()/2; l++ )
                {
                    for( int c=0; c<3; c++ )
                    {
                        divfld->operator()(i, j, k)[c] +=  settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dy() * tensor( fld(i, j+1+l, k) )[c*3+1]; //N
                        divfld->operator()(i, j, k)[c] += -settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dy() * tensor( fld(i, j-1-l, k) )[c*3+1]; //S
                        divfld->operator()(i, j, k)[c] +=  settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dx() * tensor( fld(i+1+l, j, k) )[c*3]; //E
                        divfld->operator()(i, j, k)[c] += -settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dx() * tensor( fld(i-1-l, j, k) )[c*3]; //W
                        divfld->operator()(i, j, k)[c] +=  settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dz() * tensor( fld(i, j, k+1+l) )[c*3+2]; //T
                        divfld->operator()(i, j, k)[c] += -settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dz() * tensor( fld(i, j, k-1-l) )[c*3+2]; //B
                    }
                }
            }
        }
    }

    return divfld;
}


std::shared_ptr<Field<vector> > div
( 
    const Field<tensor>& fld 
)
{
    reuseTmp<vector> divfldt( fld.mesh() );
    std::shared_ptr<Field<vector> > divfld( divfldt() );

    #pragma omp parallel for collapse(3)
    for( int i=settings::m()/2; i<fld.ni()-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<fld.nj()-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<fld.nk()-settings::m()/2; k++ )
            {
                divfld->operator()(i, j, k)=0.0;

                for( int l=0; l<settings::m()/2; l++ )
                {
                    for( int c=0; c<3; c++ )
                    {
                        divfld->operator()(i, j, k)[c] +=  settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dy() * fld(i, j+1+l, k)[c*3+1]; //N
                        divfld->operator()(i, j, k)[c] += -settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dy() * fld(i, j-1-l, k)[c*3+1]; //S
                        divfld->operator()(i, j, k)[c] +=  settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dx() * fld(i+1+l, j, k)[c*3]; //E
                        divfld->operator()(i, j, k)[c] += -settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dx() * fld(i-1-l, j, k)[c*3]; //W
                        divfld->operator()(i, j, k)[c] +=  settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dz() * fld(i, j, k+1+l)[c*3+2]; //T
                        divfld->operator()(i, j, k)[c] += -settings::coef()[l] / 2.0 / (l+1) / fld.mesh().dz() * fld(i, j, k-1-l)[c*3+2]; //B
                    }
                }
            }
        }
    }

    return divfld;
}

std::shared_ptr<Field<vector> > div
( 
    const std::shared_ptr<Field<symmTensor> > fld 
)
{
    return ex::div( (*fld) );
}


std::shared_ptr<Field<vector> > div
( 
    const std::shared_ptr<Field<tensor> > fld 
)
{
    return ex::div( (*fld) );
}

std::shared_ptr<Field<scalar> > div
( 
    const std::shared_ptr<Field<vector> > fld
)
{
    return std::shared_ptr<Field<scalar> >( ex::div( (*fld) ) );
}

template <class T>
std::shared_ptr<Field<T> > div
( 
    const std::shared_ptr<Field<vector> > fld1,
    const std::shared_ptr<Field<T> > fld2
)
{
    return std::shared_ptr<Field<T> >( ex::div( (*fld1), (*fld2) ) );
}


template <class T>
std::shared_ptr<Field<T> > div
( 
    const std::shared_ptr<Field<vector> > fld1,
    const Field<T>& fld2
)
{
    return std::shared_ptr<Field<T> >( ex::div( (*fld1), fld2 ) );
}


}
