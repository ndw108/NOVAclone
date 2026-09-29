#include "poisson.h"
#include <stdint.h>
#include <memory>


poisson::poisson( std::shared_ptr<Field<scalar>> pptr,
                  std::shared_ptr<heffte::fft3d<heffte::backend::fftw>> sharedFft)
:
    pptr_(pptr),
    n_(pptr->mesh().glob_n()),
    ni_(pptr->mesh().ni()-settings::m()-1),
    nj_(pptr->mesh().nj()-settings::m()-1),
    nk_(pptr->mesh().nk()-settings::m()-1),
    localStart_(pptr->mesh().localStart())
{
  //Define FFT grid.
  int NFAST = n_[0];
  int NMID = n_[1];
  int NSLOW = n_[2];

  nfast = NFAST;
  nmid = NMID;
  nslow = NSLOW;

  MPI_Comm world = parallelCom::world();

  npfast  = parallelCom::ni();
  npmid   = parallelCom::nj();
  npslow  = parallelCom::nk();
  ipfast  = parallelCom::i(); 
  ipmid   = parallelCom::j(); 
  ipslow  = parallelCom::k(); 
  //fft initialise grid 
  ilo = (int) 1.0*(ipfast)*nfast/npfast;
  ihi = (int) 1.0*(ipfast+1)*nfast/npfast - 1;
  jlo = (int) 1.0*(ipmid)*nmid/npmid;
  jhi = (int) 1.0*(ipmid+1)*nmid/npmid - 1;
  klo = (int) 1.0*(ipslow)*nslow/npslow;
  khi = (int) 1.0*(ipslow+1)*nslow/npslow - 1;

  //Define local coords for wavenumber calculations
  inxlo = static_cast<int> (1.0 * ipfast * nfast / npfast);
  inxhi = static_cast<int> (1.0 * (ipfast+1) * nfast / npfast) - 1;

  inylo = static_cast<int> (1.0 * ipmid * nmid / npmid);
  inyhi = static_cast<int> (1.0 * (ipmid+1) * nmid / npmid) - 1;

  inzlo = static_cast<int> (1.0 * ipslow * nslow / npslow);
  inzhi = static_cast<int> (1.0 * (ipslow+1) * nslow / npslow) - 1;

  nfft_in = (inxhi-inxlo+1) * (inyhi-inylo+1) * (inzhi-inzlo+1);


  localStart_[0] = inxlo;
  localStart_[1] = inylo;
  localStart_[2] = inzlo;

  //construct FFT
  up_N = std::make_unique<heffte::box3d<>>(std::array<int, 3> ({inxlo, inylo, inzlo})  ,std::array<int, 3> ({inxhi, inyhi, inzhi}));
  N = up_N.get();

  if( sharedFft )
  {   
      fft = sharedFft;
  }   
  else
  {   
      fft = std::make_shared<heffte::fft3d<heffte::backend::fftw>>(*N, *N, world);
  }
  
  //construct buffers and wavenumbers
  up_k2_ = std::make_unique<double[]>((fft->size_inbox()));
  k2_ = up_k2_.get();

  up_phi_ = std::make_unique<std::vector<double>>((fft->size_inbox()));
  phi_ = up_phi_.get();

  up_phihat_ = std::make_unique<std::vector<std::complex<double>>>((fft->size_outbox()));
  phihat_ = up_phihat_.get();

  wavenumbersFFT_();


//currently unused functions for DCTs and split FFTs
/*
    up_NDCT = std::make_unique<heffte::box3d<>>(std::array<int, 3> ({inxlo, inylo, inzlo})  ,std::array<int, 3> ({inxhi, inyhi, inzhi}));
    NDCT = up_NDCT.get();

    dct = std::make_unique<heffte::fft3d<heffte::backend::fftw_cos1>>(*NDCT, *NDCT, world);
    fft2d = std::make_unique<heffte::fft3d<heffte::backend::fftw>>(*NDCT, *NDCT, world);
    dct->nullify_executor(0); 
    dct->nullify_executor(2);
    fft2d->nullify_executor(1);

    fftx = std::make_unique<heffte::fft3d<heffte::backend::fftw>>(*N, *N, world);
    fftx->nullify_executor(1);
    fftx->nullify_executor(2);

    //plan 3d DCT-FFT-FFT
    //perform 1d r2r_DCT
    //translate domain from std::double to std::complex<double>
    //perform 2d r2c_DFT
    //convolute and solve
    //perform 2d c2r_IDFT
    //translate domain from std::complex<double> to std::double
    //perform 1d r2r_DCT
    //std::double cosphi
    //std::complex<double> dctphi
    //std::complex<double> dctphihat
    //heffte::fft1d<heffte::backend::fftw_cos>>
    //heffte::fft2d<heffte::backend::fftw>>

    //Create and allocate memory for FFT/DCT inputs & outputs.
    up_k2c_ = std::make_unique<double[]>((2*dct->size_inbox()));
    k2c_ = up_k2c_.get();

    up_k2cy_ = std::make_unique<double[]>((2*dct->size_inbox()));
    k2cy_ = up_k2cy_.get(); 

    up_dctworkspace_ = std::make_unique<std::vector<double>>(2*dct->size_workspace());  
    dctworkspace_ = up_dctworkspace_.get();

    up_cosphihat_ = std::make_unique<std::vector<double>>((2*dct->size_outbox()));
    cosphihat_ = up_cosphihat_.get();

    up_cosphi_ = std::make_unique<std::vector<double>>((2*dct->size_inbox()));
    cosphi_ = up_cosphi_.get();

    up_cosphicomplex_ = std::make_unique<std::vector<std::complex<double>>>((2*fft2d->size_outbox()));
    cosphicomplex_ = up_cosphicomplex_.get();
   
    up_psi_ = std::make_unique<std::vector<double>>((fftx->size_inbox()));
    psi_ = up_psi_.get();
   
    up_psihat_ = std::make_unique<std::vector<std::complex<double>>>((fftx->size_outbox()));
    psihat_ = up_psihat_.get();
  
    wavenumbersDCT_();
*/
}


void poisson::rhs( std::shared_ptr<Field<scalar> > f ) 
{   
    int l=0;
    for( int k=(settings::m()/2); k<(pptr_->nk()-1-settings::m()/2); k++ )
    {
        for( int j=(settings::m()/2); j<(pptr_->nj()-1-settings::m()/2); j++ )
        {
            for( int i=(settings::m()/2); i<(pptr_->ni()-1-settings::m()/2); i++ )
	    {
                phi_[0][l] = f->operator()(i, j, k) ; 
		l++;
            }   
        }
    }   
}

void poisson::solve()
{
    fft->forward(phi_->data(), phihat_->data());

    int l=0;
        for( int k=0; k<nk_; k++ )
        {
            for( int j=0; j<nj_; j++ )
            {
                for( int i=0; i<ni_; i++ )
                {
                if( k2_[l] < tools::eps )
                {
                    phihat_[0][l] = std::complex<double>(0.0, 0.0);
                }
                else
                {
                     phihat_[0][l] = std::complex<double>(( -std::real(phihat_[0][l]) / k2_[l]) , -std::imag(phihat_[0][l]) / k2_[l]);
                }

                l++;
            }
       }
    }
    fft->backward(phihat_->data(), phi_->data());

    l=0;
    for( int k=(settings::m()/2); k<(pptr_->nk()-1-settings::m()/2); k++ )
    {
        for( int j=(settings::m()/2); j<(pptr_->nj()-1-settings::m()/2); j++ )
        {
            for( int i=(settings::m()/2); i<(pptr_->ni()-1-settings::m()/2); i++ )
            {
                pptr_->operator()(i,j,k) = phi_[0][l] / (n_[0]*n_[1]*n_[2]);
                l++;

            }
        }
    }
}


/*
void poisson::rhsDCT( std::shared_ptr<Field<scalar> > f )
{ 

    int l=0;
    for( int k=(settings::m()/2); k<(pptr_->nk()-1-settings::m()/2); k++ )
    {   
        for( int j=(settings::m()/2); j<(pptr_->nj()-1-settings::m()/2); j++ )
        {   
            for( int i=(settings::m()/2); i<(pptr_->ni()-1-settings::m()/2); i++ )
            {                  
		cosphi_[0][l] = f->operator()(i, j, k) ; 
                l++;
            }   
        }   
    } 

}

void poisson::x_FFT( std::shared_ptr<Field<scalar> > f )
{  
    int l=0;
    for( int k=(settings::m()/2)+1; k<(settings::m()/2)+2; k++ )
    {
        for( int j=(settings::m()/2); j<(pptr_->nj()-1-settings::m()/2); j++ )
        {
            for( int i=(settings::m()/2); i<(pptr_->ni()-1-settings::m()/2); i++ )
            {
                psi_[0][l] = f->operator()(i, j, k) ;
                l++;
            }  
        }
    } 
    fftx->forward(psi_->data(), psihat_->data());
    
    std::ofstream data("fftxdat.dat");
    for(size_t ll=0; ll<((fftx->size_outbox())/nk_); ll++) {
    data<<std::real(psihat_[0][ll])<<" "<<std::setprecision(15)<<std::imag(psihat_[0][ll])<<" "<<std::setprecision(15)<<std::endl;    
    }
}

void poisson::solveDCT()
{
    dct->forward(cosphi_->data(), cosphihat_->data());

       int l=0;
        for( int k=0; k<nk_; k++ )
        {   
            for( int j=0; j<nj_; j++ )
            {   
                for( int i=0; i<ni_; i++ )
                {   
                {   
                     cosphihat_[0][l] =  (cosphihat_[0][l]) * k2cy_[l];
                }   
                l++;
                }   
            }   
        }
 
    fft2d->forward(cosphihat_->data(), cosphicomplex_->data());

	l=0;
        for( int k=0; k<nk_; k++ )
        {
            for( int j=0; j<nj_; j++ )
            {
                for( int i=0; i<ni_; i++ )
                {
		if( k2c_[l] < tools::eps )
                {
		    cosphicomplex_[0][l] = std::complex<double>(0.0 , 0.0);
		}
                else
                {
                     cosphicomplex_[0][l] =( std::complex<double>(( -std::real(cosphicomplex_[0][l]) / k2c_[l]) , -std::imag(cosphicomplex_[0][l]) / k2c_[l]  ) );
                }   
		l++;
                }   
            }   
        } 
  
    fft2d->backward(cosphicomplex_->data(), cosphihat_->data());

    dct->backward(cosphihat_->data(), cosphi_->data(), dctworkspace_->data());

    scalar pmax=0.0;
    l=0;
    for( int k=(settings::m()/2); k<(pptr_->nk()-1-settings::m()/2); k++ )
    {
        for( int j=(settings::m()/2); j<(pptr_->nj()-1-settings::m()/2); j++ )
        {
            for( int i=(settings::m()/2); i<(pptr_->ni()-1-settings::m()/2); i++ )
            {
//		pmax = std::max( pmax, (cosphi_[0][l] / (2.0*n_[0]*n_[1]*(n_[2]-1) )));
                pptr_->operator()(i,j,k) = cosphi_[0][l] / (2*(n_[0])*(n_[1]-1)*n_[2]);
                l++;

            }
        }
    }
//    std::cout<<" "<<dct->size_inbox()<<std::endl;
//    std::cout<<" "<<dct->size_outbox()<<std::endl;
//std::cout<<" "<<fft2d->size_outbox()<<std::endl;
//std::cout<<" "<<fft2d->size_outbox()<<std::endl;
//std::cout<<" "<<fft->size_outbox()<<std::endl;
//std::cout<<" "<<fft->size_inbox()<<std::endl;

}
*/

void poisson::wavenumbersFFT_()
{
    Mesh& m = pptr_->mesh();
    auto pi = tools::pi;

    
    scalar kx, ky, kz; //TODO: will be complex for non-central schemes
    int l=0;
	for( int k=0; k<nk_; k++ )
	{
	    for( int j=0; j<nj_; j++ )
	    {
		for( int i=0; i<ni_; i++ )
		{

		if (2*( i + localStart_[0] )<n_[0])
                    kx = 2.0*pi*( i + localStart_[0]  )/m.lx();
                else 
                    kx = 2.0*pi*(n_[0] - i - localStart_[0] )/m.lx();

		if (2*( j + localStart_[1] ) <n_[1])
                    ky = 2.0*pi*( j + localStart_[1] )/m.ly();
                else 
                    ky = 2.0*pi*(n_[1] - j - localStart_[1] )/m.ly();


                    kz = 2.0*pi*( k + localStart_[2] )/m.lz();

                //modified wavenumbers
                scalar kpx=0.0;
                scalar kpy=0.0;
                scalar kpz=0.0;
                if(kx > tools::eps )
                {
                    for( int ll=0; ll<(settings::m()/2); ll++ )
                    {
                        kpx+=kx*settings::coef()[ll]/(ll+1)*std::sin((ll+1)*kx*m.dx())/(kx*m.dx());
                    }
                }
                if(ky > tools::eps )
                {
                    for( int ll=0; ll<(settings::m()/2); ll++ )
                    {
                        kpy+=ky*settings::coef()[ll]/(ll+1)*std::sin((ll+1)*ky*m.dy())/(ky*m.dy());
                    }
                }
                if(kz > tools::eps )
                {
                    for( int ll=0; ll<(settings::m()/2); ll++ )
                    {
                        kpz+=kz*settings::coef()[ll]/(ll+1)*std::sin((ll+1)*kz*m.dz())/(kz*m.dz());
                    }
                }
                
                k2_[l] = kpx*kpx + kpy*kpy + kpz*kpz;
                l++;
	    }
        }
    }
}

/*
void poisson::wavenumbersDCT_()
{
    Mesh& m = pptr_->mesh();
    auto pi = tools::pi;


    scalar kx, ky, kz; 
    int l=0;
    kx=0;
    ky=0;
    kz=0;
        for( int k=0; k<nk_; k++ )
        {
            for( int j=0; j<nj_; j++ )
            {
                for( int i=0; i<ni_; i++ )
                {
                if (2*( i + localStart_[0] )<n_[0])
                    kx = 2.0*pi*( i + localStart_[0]  )/m.lx();
                else 
                    kx = pi*(n_[0] - i - localStart_[0] )/m.lx();

                if (2*( j + localStart_[1] ) <n_[1])
                    ky = pi*( j + localStart_[1] )/m.ly();
                else 
                    ky = pi*(n_[1] - j - localStart_[1] )/m.ly();


                    kz = 2.0*pi*( k + localStart_[2] )/m.lz();

                

		//modified wavenumbers
                scalar kpx=0.0;
                scalar kpy=0.0;
                scalar kpz=0.0;
                if(kx > tools::eps )
                {   
                    for( int ll=0; ll<(settings::m()/2); ll++ )
                    {   
                        kpx+=kx*settings::coef()[ll]/(ll+1)*std::sin((ll+1)*kx*m.dx())/(kx*m.dx());
                    }   
                }   
                if(ky > tools::eps )
                {   
                    for( int ll=0; ll<(settings::m()/2); ll++ )
                    {   
                        kpy+=ky*settings::coef()[ll]/(ll+1)*std::sin((ll+1)*ky*m.dy())/(ky*m.dy());
                    }   
                }   
                if(kz > tools::eps )
                {   
                    for( int ll=0; ll<(settings::m()/2); ll++ )
                    {   
                        kpz+=kz*settings::coef()[ll]/(ll+1)*std::sin((ll+1)*kz*m.dz())/(kz*m.dz());
                    }   
                }    

		k2cy_[l] = (kpy/(n_[1]-1)) * (kpy/(n_[1]-1) ); 
                k2c_[l] = kpx*kpx + kpz*kpz;
                l++;
            }   
        }   
    }   
}
*/

poisson::~poisson() = default;

