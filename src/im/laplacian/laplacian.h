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
#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include "im/fdMatrix/fdMatrix.h"
#include "Types/vector/vector.h"
#include "im/laplacian/laplacian.h"

namespace im
{

template <class T>
std::shared_ptr<fdMatrix<T> > laplacian
(
    const scalar,
    std::shared_ptr<Field<T> >
);


template <class T>
std::shared_ptr<fdMatrix<T> > laplacian
(
    std::shared_ptr<Field<T> >
);

}

#include "im/laplacian/laplacian.hxx"

#endif
