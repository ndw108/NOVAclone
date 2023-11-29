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
#include "fft3d.h"
#include "Mesh/Mesh.h"
#include "Field/Field.h"
#include "Tools/Tools.h"
#include <memory>

using namespace FFTMPI_NS;

class poisson
{
    protected:
    std::unique_ptr<FFT3d> fft;
    std::shared_ptr<Field<scalar> > pptr_; 
    int fftsize;
    unsigned ni_, nj_, nk_;
    ptrdiff_t* n_;
    ptrdiff_t* local_ni_;
    ptrdiff_t* local_i_start_;
    ptrdiff_t* local_no_;
    ptrdiff_t* local_o_start_;
//    ptrdiff_t alloc_local_;

    typedef int64_t bigint;
    int cflag,eflag,pflag,mflag,sflag,rflag; 
//    double* sendbuf;
//    double* recvbuf;
    int nfft_in;
    int nfft_out;
    int nfast;
    int nmid;
    int nslow;
    int in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi;
    int out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi;
    int permute,sendsize,recvsize;
    int flag,niter,tflag;
    double tmax; 
    int memorysize;

    std::unique_ptr<FFT_SCALAR[]> up_phihat_;
    std::unique_ptr<FFT_SCALAR[]> up_phi_;
    FFT_SCALAR* phi_;
    FFT_SCALAR* phihat_;
//    std::array<double*, 2> phihat_; 
//    std::array<double*, 2> phi_;
    std::unique_ptr<double[]> up_k2_;
    double* k2_;


    public:
    poisson( std::shared_ptr<Field<scalar> > );
    
    void wavenumbers_();
    void rhs(std::shared_ptr<Field<scalar> > f );
    void solve();
};

#endif



