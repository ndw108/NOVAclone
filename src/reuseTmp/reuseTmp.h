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
#ifndef REUSETMP_H
#define REUSETMP_H

#include <memory>
#include "Field/Field.h"
#include "BC/BC.h"

template <class T>
class Field;

template <class T>
class reuseTmp
{
    static std::vector<std::shared_ptr<Field<T> > > tmp_;
    static int index_;

public:

    reuseTmp( Mesh ); 

    std::shared_ptr<Field<T> > operator()()
    {
        return tmp_[index_];
    } 
};


#include "reuseTmp/reuseTmp.hxx"

#endif
