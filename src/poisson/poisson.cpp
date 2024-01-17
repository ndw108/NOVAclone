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
    localStart_(pptr->mesh().localStart())
{
fft = std::make_unique<FFT3d>(parallelCom::world(), 2);
  //fft options
  //tflag=0;
  permute = 0;

  NFAST = n_[0];
  NMID = n_[1];
  NSLOW = n_[2];

//  int MPI_Init( int *argc, char ***argv );
      MPI_Comm world = parallelCom::world();
  
  nfast = NFAST;
  nmid = NMID;
  nslow = NSLOW;

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

  //define local coords for wavenumbers
  inxlo = static_cast<int> (1.0 * ipfast * nfast / npfast);
  inxhi = static_cast<int> (1.0 * (ipfast+1) * nfast / npfast) - 1;

  inylo = static_cast<int> (1.0 * ipmid * nmid / npmid);
  inyhi = static_cast<int> (1.0 * (ipmid+1) * nmid / npmid) - 1;

  inzlo = static_cast<int> (1.0 * ipslow * nslow / npslow);
  inzhi = static_cast<int> (1.0 * (ipslow+1) * nslow / npslow) - 1;

  nfft_in = (inxhi-inxlo+1) * (inyhi-inylo+1) * (inzhi-inzlo+1);

{
  for (int iproc = 0; iproc < parallelCom::worldSize(); iproc++) {
    if (parallelCom::myProcNo() != iproc) continue;
    if (parallelCom::myProcNo() >= 1) MPI_Recv(&tmp,0,MPI_INT,parallelCom::myProcNo()-1,0,world,MPI_STATUS_IGNORE);
    {
        int nxlocal = inxhi - inxlo + 1; 
        int nylocal = inyhi - inylo + 1; 

        for (int m = 0; m < nfft_in; m++) {
          ilocal = m % nxlocal;
          jlocal = (m/nxlocal) % nylocal;
          klocal = m / (nxlocal*nylocal);
          iglobal = inxlo + ilocal;
          jglobal = inylo + jlocal;
          kglobal = inzlo + klocal;
          printf("Value (%d,%d,%d) on proc %d \n",
               iglobal,jglobal,kglobal,
               parallelCom::myProcNo());
	if (m == 0)
	{
	localStart_[0] = iglobal;
        localStart_[1] = jglobal;
        localStart_[2] = kglobal;
	}
      }
    }
  if (parallelCom::myProcNo() < parallelCom::worldSize()) MPI_Send(&tmp,0,MPI_INT,parallelCom::myProcNo()+1,0,world);
  } 
}

std::cout<<ni_<<" "<<nj_<<" "<<nk_<<" ";
std::cout<<ilo<<" "<<ihi<<" "<<jlo<<" "<<jhi<<" "<<klo<<" "<<khi<<" "<<npfast<<" "<<npmid<<" "<<npslow<<" "<<nfast<<" "<<nmid<<" "<<nslow<<" ";

    //fft plan
    fft->setup(nfast, nmid, nslow,             // 3d verion
          ilo, ihi, jlo, 
          jhi, klo, khi,
          ilo, ihi, jlo, 
          jhi, klo, khi,
          permute, fftsize, sendsize, recvsize);


    //fft allocate memory
    std::cout<<fftsize<<" ";
        
    up_phi_ = std::make_unique_for_overwrite<FFT_SCALAR[]>(2*fftsize);
    phi_ = up_phi_.get();

    up_k2_ = std::make_unique<double[]>(2*fftsize);
    k2_ = up_k2_.get();
    
    wavenumbers_();
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
                phi_[l] = (double) f->operator()(i, j, k); 
		l++;
		phi_[l] = 0.0;
		l++;
            }   
        }   
    }
    fft->compute(phi_,phi_,1);
}


void poisson::solve()
{
   int l=0;
    for( int k=0; k<nk_; k++ )
    {   
        for( int j=0; j<nj_; j++ )
        {   
	    for( int i=0; i<ni_; i++ )
	    {
		if( k2_[l] < tools::eps )
                {
                    phi_[l] = 0.0;
		}
                else
                {
                    phi_[l] = -phi_[l] / k2_[l];
		}
		l++;	
           
		if( k2_[l] < tools::eps )
                {   
                    phi_[l] = 0.0;
                }   
                else
                {   
                    phi_[l] = -phi_[l] / k2_[l];
                }   
                l++; 

	    }
        }
    }

    fft->compute(phi_,phi_,-1);
    
    l=0;
    for( int k=(settings::m()/2); k<(pptr_->nk()-1-settings::m()/2); k++ )
    {
        for( int j=(settings::m()/2); j<(pptr_->nj()-1-settings::m()/2); j++ )
        {
            for( int i=(settings::m()/2); i<(pptr_->ni()-1-settings::m()/2); i++ )
	    {
		pptr_->operator()(i,j,k) = phi_[l];
    		l++;
		l++;
	
	    }
        }
    }

}

void poisson::wavenumbers_()
{
    Mesh& m = pptr_->mesh();
    auto pi = tools::pi;

    printf("Assigned local start coordinates (%d,%d,%d) on proc %d \n",
                localStart_[0],localStart_[1],localStart_[2],parallelCom::myProcNo());
    
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
		k2_[l] = kpx*kpx + kpy*kpy + kpz*kpz;
		l++;
	    }
        }
    }
}

poisson::~poisson()
{
	free(phi_);
}
