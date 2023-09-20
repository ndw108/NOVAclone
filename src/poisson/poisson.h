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

#include "pfft.h"
#include "Mesh/Mesh.h"
#include "Field/Field.h"
#include "Tools/Tools.h"
#include <memory>

class poisson
{
    protected:
    pfft_plan r2c_;
    pfft_plan c2r_;
    std::shared_ptr<Field<scalar> > pptr_; 
    unsigned ni_, nj_, nk_;
    ptrdiff_t* n_;
    ptrdiff_t* local_ni_;
    ptrdiff_t* local_i_start_;
    ptrdiff_t* local_no_;
    ptrdiff_t* local_o_start_;
    ptrdiff_t alloc_local_;

    std::array<pfft_complex*, 2> phihat_; 
    std::array<double*,2> phi_;
    std::unique_ptr<double[]> up_k2_;
    double* k2_;

    void wavenumbers_();

    public:
    poisson( std::shared_ptr<Field<scalar> > );
    ~poisson();

    void rhs(std::shared_ptr<Field<scalar> > f );
    void solve();
};

#endif



