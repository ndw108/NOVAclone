#include "poisson.h"

poisson::poisson( std::shared_ptr<Field<scalar> > pptr )
:
    pptr_(pptr)
{
    Mesh& m = pptr->mesh();

    for( int i=0; i<2; i++ )
    {
        phi_[i] = fftw_alloc_real(m.glob_ni()*m.glob_nj()*m.glob_nk());
        phihat_[i] = fftw_alloc_complex(m.glob_ni()*m.glob_nj()*std::floor(m.glob_nk()/2.0)+1);
    }
    r2c_ = fftw_plan_dft_r2c_3d(m.glob_ni(), m.glob_nj(), m.glob_nk(), phi_[0], phihat_[0], FFTW_MEASURE);
    c2r_ = fftw_plan_dft_c2r_3d(m.glob_ni(), m.glob_nj(), m.glob_nk(), phihat_[1], phi_[1], FFTW_MEASURE);
}

void poisson::rhs( std::shared_ptr<Field<scalar> > f )
{   
    int l=0;
    for( int i=settings::m()/2; i<f->ni()-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<f->nj()-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<f->nk()-settings::m()/2; k++ )
            { 
                phi_[0][l++] = f->operator()(i, j, k);
            }
        }
    }

    fftw_execute(r2c_);
}

void poisson::solve()
{
    Mesh& m = pptr_->mesh();
    int ninjnk = m.glob_ni()*m.glob_nj()*m.glob_nk();
    auto pi = tools::pi;

    scalar kx, ky, kz;

    int l=0;
    for( int i=0; i<m.glob_ni(); i++ )
    {
        for( int j=0; j<m.glob_nj(); j++ )
        {
            for( int k=0; k<std::floor(m.glob_nk()/2.0)+1; k++ )
            { 
                if (2*i<m.glob_ni())
                    kx = 2.0*pi*i/m.lx();
                else 
                    kx = 2.0*pi*(m.glob_ni() - i)/m.lx();

                if (2*j<m.glob_nj())
                    ky = 2.0*pi*j/m.ly();
                else 
                    ky = 2.0*pi*(m.glob_nj() - j)/m.ly();

                kz = 2.0*pi*k/m.lz();

                scalar k2 = kx*kx + ky*ky + kz*kz;

                if( k2 < tools::eps )
                {
                    phihat_[1][l][0] = 0.0;
                    phihat_[1][l][1] = 0.0;
                }
                else
                {
                    phihat_[1][l][0] = -phihat_[0][l][0] / k2;
                    phihat_[1][l][1] = -phihat_[0][l][1] / k2;
                }
                l++;
            }
        }
    }

    fftw_execute(c2r_);

    l=0;
    for( int i=settings::m()/2; i<pptr_->ni()-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<pptr_->nj()-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<pptr_->nk()-settings::m()/2; k++ )
            { 
                pptr_->operator()(i,j,k) = phi_[1][l++] / ninjnk;
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
