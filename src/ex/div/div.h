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
#ifndef EXPDIV_H
#define EXPDIV_H

#if defined(_OPENMP)
#include <omp.h>
#endif

#include<memory>
#include<iostream>
#include <stdlib.h>

#include "Types/vector/vector.h"
#include "Types/tensor/tensor.h"

#include "Field/Field.h"
#include "reuseTmp/reuseTmp.h"

namespace ex
{

std::shared_ptr<Field<scalar> > div
(
    const Field<vector>&
);



template<class T>
std::shared_ptr<Field<T> > div
(
    const Field<vector>&,
    const Field<T>&
);



std::shared_ptr<Field<vector> > div
(
    const Field<scalar>&,
    const std::shared_ptr<Field<tensor> >
);

std::shared_ptr<Field<vector> > div
(
    const Field<scalar>&,
    const std::shared_ptr<Field<symmTensor> >
);

std::shared_ptr<Field<vector> > div
( 
    const Field<symmTensor>& fld 
);

std::shared_ptr<Field<vector> > div
( 
    const std::shared_ptr<Field<symmTensor> > fld 
);

std::shared_ptr<Field<vector> > div
( 
    const Field<tensor>& fld 
);

std::shared_ptr<Field<vector> > div
( 
    const std::shared_ptr<Field<tensor> > fld 
);



}

#include "ex/div/div.hxx"

#endif
