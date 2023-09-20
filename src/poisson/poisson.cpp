#include "poisson.h"
#include <memory>

poisson::poisson( std::shared_ptr<Field<scalar> > pptr )
:
    pptr_(pptr),
    n_(pptr->mesh().glob_n()),
    ni_(pptr->mesh().ni()-settings::m()-1),
    nj_(pptr->mesh().nj()-settings::m()-1),
    nk_(pptr->mesh().nk()-settings::m()-1),
    alloc_local_(pptr->mesh().alloc_local()),
    local_ni_(pptr->mesh().local_ni()),
    local_i_start_(pptr->mesh().local_i_start()),
    local_no_(pptr->mesh().local_no()),
    local_o_start_(pptr->mesh().local_o_start())
{
    for( int i=0; i<2; i++ )
    {
        phi_[i] = pfft_alloc_real(2*alloc_local_);
        phihat_[i] = pfft_alloc_complex(alloc_local_);
    }
    r2c_ = pfft_plan_dft_r2c_3d( n_, phi_[0], phihat_[0], parallelCom::pfftcomm(), PFFT_FORWARD, PFFT_TRANSPOSED_NONE | PFFT_MEASURE);
    c2r_ = pfft_plan_dft_c2r_3d( n_, phihat_[1], phi_[1], parallelCom::pfftcomm(), PFFT_BACKWARD, PFFT_TRANSPOSED_NONE | PFFT_MEASURE);

    up_k2_ = std::make_unique<double[]>( alloc_local_ );
    k2_ = up_k2_.get();
    
    wavenumbers_();
}

void poisson::rhs( std::shared_ptr<Field<scalar> > f )
{   
    int l=0;
    for( int i=settings::m()/2; i<pptr_->ni()-1-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<pptr_->nj()-1-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<pptr_->nk()-1-settings::m()/2; k++ )
            { 
                phi_[0][l++] = f->operator()(i, j, k);
            }
        }
    }

    pfft_execute(r2c_);
}

void poisson::solve()
{
    int l=0;
    for( int i=0; i<local_no_[0]; i++ )
    {
        for( int j=0; j<local_no_[1]; j++ )
        {
            for( int k=0; k<local_no_[2]; k++ )
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

    pfft_execute(c2r_);

    l=0;
    for( int i=settings::m()/2; i<pptr_->ni()-1-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<pptr_->nj()-1-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<pptr_->nk()-1-settings::m()/2; k++ )
            { 
                pptr_->operator()(i,j,k) = phi_[1][l++] / (n_[0]*n_[1]*n_[2]);
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
    for( int i=0; i<local_no_[0]; i++ )
    {
        for( int j=0; j<local_no_[1]; j++ )
        {
            for( int k=0; k<local_no_[2]; k++ )
            { 
                if (2*(i+local_o_start_[0])<n_[0])
                    kx = 2.0*pi*(i+local_o_start_[0])/m.lx();
                else 
                    kx = 2.0*pi*(n_[0] - i - local_o_start_[0])/m.lx();

                if (2*(j+local_o_start_[1])<n_[1])
                    ky = 2.0*pi*(j+local_o_start_[1])/m.ly();
                else 
                    ky = 2.0*pi*(n_[1] - j - local_o_start_[1])/m.ly();


                kz = 2.0*pi*(k+local_o_start_[2])/m.lz();

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
        pfft_free(phi_[i]);
        pfft_free(phihat_[i]);
    }
}
