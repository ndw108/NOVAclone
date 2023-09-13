#include "poisson.h"
#include <memory>

poisson::poisson( std::shared_ptr<Field<scalar> > pptr )
:
    pptr_(pptr),
    ni_(pptr->mesh().glob_ni()-1),
    nj_(pptr->mesh().glob_nj()-1),
    nk_(pptr->mesh().glob_nk()-1),
    up_k2_( std::make_unique<double[]>( ni_*nj_*(std::floor(nk_/2.0)+1) ) ),
    k2_(up_k2_.get())
{
    for( int i=0; i<2; i++ )
    {
        phi_[i] = fftw_alloc_real(ni_*nj_*nk_);
        phihat_[i] = fftw_alloc_complex(ni_*nj_*(std::floor(nk_/2.0)+1));
    }
    r2c_ = fftw_plan_dft_r2c_3d(ni_, nj_, nk_, phi_[0], phihat_[0], FFTW_MEASURE);
    c2r_ = fftw_plan_dft_c2r_3d(ni_, nj_, nk_, phihat_[1], phi_[1], FFTW_MEASURE);

    wavenumbers_();
}

void poisson::rhs( std::shared_ptr<Field<scalar> > f )
{   
    int l=0;
    for( int i=settings::m()/2; i<f->ni()-settings::m()/2-1; i++ )
    {
        for( int j=settings::m()/2; j<f->nj()-settings::m()/2-1; j++ )
        {
            for( int k=settings::m()/2; k<f->nk()-settings::m()/2-1; k++ )
            { 
                phi_[0][l++] = f->operator()(i, j, k);
            }
        }
    }

    fftw_execute(r2c_);
}

void poisson::solve()
{
    int l=0;
    for( int i=0; i<ni_; i++ )
    {
        for( int j=0; j<nj_; j++ )
        {
            for( int k=0; k<std::floor(nk_/2.0)+1; k++ )
            { 
                if( k2_[l] < tools::eps )
                {
                    phihat_[1][l][0] = 0.0;
                    phihat_[1][l][1] = 0.0;
                }
                else
                {
                    phihat_[1][l][0] = -phihat_[0][l][0] / k2_[l];
                    phihat_[1][l][1] = -phihat_[0][l][1] / k2_[l];
                }
                l++;
            }
        }
    }
    phihat_[1][0][0] = 0.0;
    phihat_[1][0][1] = 0.0;


    fftw_execute(c2r_);

    l=0;
    for( int i=settings::m()/2; i<pptr_->ni()-1-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<pptr_->nj()-1-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<pptr_->nk()-1-settings::m()/2; k++ )
            { 
                pptr_->operator()(i,j,k) = phi_[1][l++] / (ni_*nj_*nk_);
            }
        }
    }

}

void poisson::wavenumbers_()
{
    Mesh& m = pptr_->mesh();
    auto pi = tools::pi;

    scalar kx, ky, kz; //TODO: will be complex for non-central schemes
    int l=0;
    for( int i=0; i<ni_; i++ )
    {
        for( int j=0; j<nj_; j++ )
        {
            for( int k=0; k<std::floor(nk_/2.0)+1; k++ )
            { 
                if (2*i<ni_)
                    kx = 2.0*pi*i/m.lx();
                else 
                    kx = 2.0*pi*(ni_ - i)/m.lx();

                if (2*j<nj_)
                    ky = 2.0*pi*j/m.ly();
                else 
                    ky = 2.0*pi*(nj_ - j)/m.ly();


                kz = 2.0*pi*k/m.lz();

                //modified wavenumbers
                scalar kpx=0.0;
                scalar kpy=0.0;
                scalar kpz=0.0;
                if(kx > tools::eps )
                {
                    for( int ll=0; ll<settings::m()/2; ll++ )
                    {
                        kpx+=kx*settings::coef()[ll]/(ll+1)*std::sin((ll+1)*kx*m.dx())/(kx*m.dx());
                    }
                }
                if(ky > tools::eps )
                {
                    for( int ll=0; ll<settings::m()/2; ll++ )
                    {
                        kpy+=ky*settings::coef()[ll]/(ll+1)*std::sin((ll+1)*ky*m.dy())/(ky*m.dy());
                    }
                }
                if(kz > tools::eps )
                {
                    for( int ll=0; ll<settings::m()/2; ll++ )
                    {
                        kpz+=kz*settings::coef()[ll]/(ll+1)*std::sin((ll+1)*kz*m.dz())/(kz*m.dz());
                    }
                }

                k2_[l++] = kpx*kpx + kpy*kpy + kpz*kpz;
            }
        }
    }
}


poisson::~poisson()
{
    for( int i=0; i<2; i++ )
    {
        fftw_free(phi_[i]);
        fftw_free(phihat_[i]);
    }
}
