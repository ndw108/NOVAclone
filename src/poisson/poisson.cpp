#include "poisson.h"
#include <stdint.h>
#include <memory>


poisson::poisson( std::shared_ptr<Field<scalar> > pptr )
:
    pptr_(pptr),
    n_(pptr->mesh().glob_n()),
    ni_(pptr->mesh().ni()-settings::m()-1),
    nj_(pptr->mesh().nj()-settings::m()-1),
    nk_(pptr->mesh().nk()-settings::m()-1),
//    alloc_local_(pptr->mesh().alloc_local()),
    local_ni_(pptr->mesh().local_ni()),
    local_i_start_(pptr->mesh().local_i_start()),
    local_no_(pptr->mesh().local_no()),
    local_o_start_(pptr->mesh().local_o_start())
{
fft = std::make_unique<FFT3d>(parallelCom::world(), 2);
//    fft options
//    tflag=0;
//    permute = 0;
//    memoryflag 1 for autoallocate, 0 for user allocate
//    fft->memoryflag = 0;
//    memorysize = 1; 

    //fft initialise grid 
     in_ilo = (int) 0;
     in_ihi = (int) (ni_)-1;
     in_jlo = (int) 0;
     in_jhi = (int) (nj_)-1;
     in_klo = (int) 0;
     in_khi = (int) (nk_)-1;
     out_ilo = (int) 0;
     out_ihi = (int) (ni_)-1;
     out_jlo = (int) 0;
     out_jhi = (int) (nj_)-1;
     out_klo = (int) 0;
     out_khi = (int) (nk_)-1;
//    nfft_in = (in_ihi-in_ilo+1) * (in_jhi-in_jlo+1) * (in_khi-in_klo+1);
  nfft_in = ni_*nj_*nk_;
  nfft_out = ni_*nj_*nk_;  
//  nfft_out = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) * (out_khi-out_klo+1);

std::cout<<ni_<<" "<<nj_<<" "<<nk_<<" ";
    //fft plan
    
    fft->setup(ni_, nj_, nk_,             // 3d verion
          in_ilo, in_ihi, in_jlo, 
          in_jhi, in_klo, in_khi,
          out_ilo, out_ihi, out_jlo, 
          out_jhi, out_klo, out_khi,
          permute, fftsize, sendsize, recvsize);


    //fft allocate memory

//    FFT_SCALAR *sendbuf = (FFT_SCALAR *) malloc(sendsize*sizeof(FFT_SCALAR));
//    FFT_SCALAR *recvbuf = (FFT_SCALAR *) malloc(recvsize*sizeof(FFT_SCALAR));
//    fft->setup_memory(sendbuf,recvbuf); 

//  bigint nbytes = ((bigint) sizeof(FFT_SCALAR)) * 2*fftsize;
//  phi_ = (FFT_SCALAR *) malloc(nbytes);
//  phihat_ = (FFT_SCALAR *) malloc(nbytes);
//  if (nbytes && phi_ == NULL) printf("Failed malloc for FFT grid");

std::cout<<fftsize<<" ";
//   for( int i=0; i<2; i++ )
//    {
        up_phi_ = std::make_unique_for_overwrite<FFT_SCALAR[]>(fftsize*2);
	up_phihat_ = std::make_unique_for_overwrite<FFT_SCALAR[]>(fftsize*2);
	phi_ = up_phi_.get();
	phihat_ = up_phihat_.get();
//    }


    up_k2_ = std::make_unique<double[]>(fftsize);
    k2_ = up_k2_.get();
    
    wavenumbers_();
}

//void poisson::rhs( std::shared_ptr<Field<scalar> > f )
//{   
//    int l=0;
//    for( int i=settings::m()/2; i<pptr_->ni()-1-settings::m()/2; i++ )
//    {
//        for( int j=settings::m()/2; j<pptr_->nj()-1-settings::m()/2; j++ )
//        {
//            for( int k=settings::m()/2; k<pptr_->nk()-1-settings::m()/2; k++ )
//            { 
//                phi_[l++]= f->operator()(i, j, k);
//                
//	    }
//        }
//    }
//  fft->compute(phi_,phihat_,1);
//}

void poisson::rhs( std::shared_ptr<Field<scalar> > f ) 
{   
    int l=0;
    for( int i=in_ilo; i<=in_ihi; i++ )
    {   
        for( int j=in_jlo; j<=in_jhi; j++ )
        {   
            for( int k=in_klo; k<=in_khi; k++ )
            {   
                phi_[l++]= f->operator()(i, j, k); 
    
            }   
        }   
    }   
  fft->compute(phi_,phihat_,1);
}


void poisson::solve()
{
    int l=0;
    for( int i=in_ilo; i<=in_ihi; i++ )
    {
        for( int j=in_jlo; j<=in_jhi; j++ )
        {
            for( int k=in_klo; k<=in_khi; k++ )
            { 
                if( k2_[l] < tools::eps )
                {
                    phihat_[l] = 0.0;
                }
                else
                {
                    phihat_[l] = -phihat_[l] / k2_[l];
                }
                l++;
		if( k2_[l] < tools::eps )
                {   
                    phihat_[l] = 0.0;
                }   
                else
                {   
                    phihat_[l] = -phihat_[l] / k2_[l];
                }   
                l++;
            }
        }
    }

   fft->compute(phihat_,phi_,-1);

    l=0;
    for( int i=settings::m()/2; i<pptr_->ni()-1-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<pptr_->nj()-1-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<pptr_->nk()-1-settings::m()/2; k++ )
            { 
                pptr_->operator()(i,j,k) = phi_[l++] / (n_[0]*n_[1]*n_[2]);
            
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
