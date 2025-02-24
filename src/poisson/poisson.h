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

#ifndef POISSON_H
#define POISSON_H

#include <stdint.h>
#include "heffte.h"
#include <heffte_common.h>

#include "Mesh/Mesh.h"
#include "Field/Field.h"
#include "Tools/Tools.h"
#include <mpi.h>
#include <memory>
#include <complex>

using namespace heffte;

class poisson
{
    protected:
    std::unique_ptr<heffte::fft3d<heffte::backend::fftw>> fft;
    std::unique_ptr<heffte::fft3d<heffte::backend::fftw>> fft2d;
    std::unique_ptr<heffte::fft3d<heffte::backend::fftw_cos1>> dct;
    std::unique_ptr<heffte::fft3d<heffte::backend::fftw>> fftx;
    std::shared_ptr<Field<scalar> > pptr_; 
    int fftsize;
    unsigned ni_, nj_, nk_;
    ptrdiff_t* n_;
    ptrdiff_t* localStart_;
    int npfast,npmid,npslow;
    int ipfast,ipmid,ipslow;
    int ilocal,jlocal,klocal,iglobal,jglobal,kglobal;
    int inxlo,inxhi,inylo,inyhi,inzlo,inzhi;
    typedef int64_t bigint;
    int cflag,eflag,pflag,mflag,sflag,rflag; 
    int NFAST,NMID,NSLOW;
    int nfast,nmid,nslow;
    int* tmp;
    int nfft_in;
    int nfft_out;
    int ilo,ihi,jlo,jhi,klo,khi;
    int permute,exchange,sendsize,recvsize;
    int flag,niter,tflag;
    double tmax; 
    int memorysize;
 
    std::unique_ptr<std::vector<std::complex<double>>> up_phihat_;
    std::vector<std::complex<double>>* phihat_; 

    std::unique_ptr<std::vector<double>> up_cosphi_;
    std::vector<double>* cosphi_;

    std::unique_ptr<std::vector<std::complex<double>>> up_cosphicomplex_;
    std::vector<std::complex<double>>* cosphicomplex_;

    std::unique_ptr<std::vector<double>> up_cosphihat_;
    std::vector<double>* cosphihat_;

    std::unique_ptr<std::vector<double>> up_dctworkspace_;
    std::vector<double>* dctworkspace_;

    std::unique_ptr<std::vector<double>> up_phi_;
    std::vector<double>* phi_;

    std::unique_ptr<std::vector<std::complex<double>>> up_psihat_;
    std::vector<std::complex<double>>* psihat_;

    std::unique_ptr<std::vector<double>> up_psi_;
    std::vector<double>* psi_;

    std::unique_ptr<double[]> up_k2cy_;
    double* k2cy_;

    std::unique_ptr<double[]> up_k2_;
    double* k2_;

    std::unique_ptr<double[]> up_k2c_;
    double* k2c_;

    heffte::box3d<>* N;
    std::unique_ptr<heffte::box3d<>> up_N;

    heffte::box3d<>* NDCT;
    std::unique_ptr<heffte::box3d<>> up_NDCT;

    public:

    poisson( std::shared_ptr<Field<scalar> > );
    ~poisson();
    void wavenumbersFFT_();
    void wavenumbersDCT_();
    void x_FFT( std::shared_ptr<Field<scalar> > f );    
    void rhs(std::shared_ptr<Field<scalar> > f );
    void rhsDCT(std::shared_ptr<Field<scalar> > f );
    void solve();
    void solveDCT();
};

#endif


