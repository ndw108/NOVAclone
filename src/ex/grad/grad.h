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
#ifndef EXPGRAD_H
#define EXPGRAD_H

#if defined(_OPENMP)
#include <omp.h>
#endif

#include<memory>
#include<iostream>
#include <stdlib.h>

#include "Types/vector/vector.h"
#include "Types/tensor/tensor.h"
#include "Field/Field.h"

namespace ex
{

std::shared_ptr<Field<vector> > grad
(
    const Field<scalar>&
);


std::shared_ptr<Field<tensor> > grad
(
    const Field<vector>&
);


std::shared_ptr<Field<tensor> > grad
(
    const std::shared_ptr<Field<vector> >
);

std::shared_ptr<Field<vector> > grad
(
    const std::shared_ptr<Field<scalar> >
);

}

#endif
