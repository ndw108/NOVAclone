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
#include "Field/Field.h"

std::shared_ptr<Field<symmTensor> > twoSymm
(
    std::shared_ptr<Field<tensor> > fld
)
{
    reuseTmp<symmTensor> rest( fld->mesh() );
    std::shared_ptr<Field<symmTensor> > res( rest() );

    int nn = fld->ni()*fld->nj()*fld->nk();
    tensor* ptrf = fld->ptr();
    symmTensor* ptrr = res->ptr();

    #pragma omp parallel for 
    for( int i=0; i<nn; i++ )
    {
       ptrr[i] = twoSymm( ptrf[i] );
    }
    
    return res;
}

std::shared_ptr<Field<symmTensor> > symm
(
    std::shared_ptr<Field<tensor> > fld
)
{
    reuseTmp<symmTensor> rest( fld->mesh() );
    std::shared_ptr<Field<symmTensor> > res( rest() );

    int nn = fld->ni()*fld->nj()*fld->nk();
    tensor* ptrf = fld->ptr();
    symmTensor* ptrr = res->ptr();

    #pragma omp parallel for 
    for( int i=0; i<nn; i++ )
    {
       ptrr[i] = symm( ptrf[i] );
    }
    
    return res;
}

std::shared_ptr<Field<scalar> > operator&&
(
    std::shared_ptr<Field<tensor> > fld1,
    std::shared_ptr<Field<tensor> > fld2
)
{
    reuseTmp<scalar> rest( fld1->mesh() );
    std::shared_ptr<Field<scalar> > res( rest() );

    int nn = fld1->ni()*fld1->nj()*fld1->nk();
    scalar* ptrr = res->ptr();
    tensor* ptr1 = fld1->ptr();
    tensor* ptr2 = fld2->ptr();

    #pragma omp parallel for 
    for( int i=0; i<nn; i++ )
    {  
       ptrr[i] = ptr1[i] && ptr2[i]; 
    }

    return res;
}

