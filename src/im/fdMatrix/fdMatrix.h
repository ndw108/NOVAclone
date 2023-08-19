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
#ifndef FD_H
#define FD_H

#include <memory>
#include<vector>
#include "Types/scalar/scalar.h"
#include "pTraits/pTraits.h"
#include "Field/Field.h"
#include "Tools/Tools.h"

namespace im
{

template<class T>
class fdMatrix : public Field<T> 
{
protected:
    std::vector<std::shared_ptr<Field<T> > > al_;
    std::shared_ptr<Field<T> > ap_;
    std::shared_ptr<Field<T> > transportedFld_;

public:
    //constructor
    fdMatrix( std::shared_ptr<Field<T> > );
    fdMatrix( std::shared_ptr<fdMatrix<T> > ); 

    auto& al( int i )
    {
        return (*al_[i]);
    }

    const auto& al( int i ) const
    {
        return (*al_[i]);
    }

    auto& ap()
    {
        return (*ap_);
    }
 

    const auto& ap() const
    {
        return (*ap_);
    }
   
    void solve( scalar, int );
    void deferredCorrection( const std::shared_ptr<Field<T> > );
    void relax( const scalar& );
    T residual();
    void setReference( int, int, int, scalar );

    fdMatrix<T>& operator=( const fdMatrix<T>& );
    fdMatrix<T>& operator+=( const fdMatrix<T>& ); 
    fdMatrix<T>& operator-=( const fdMatrix<T>& ); 
    fdMatrix<T>& operator*=( const fdMatrix<T>& );
    void operator==( const Field<T>& ); 
};


//---------
//global operators
template<class T>
std::shared_ptr<fdMatrix<T> > operator+
(
    const std::shared_ptr<fdMatrix<T> >,
    const std::shared_ptr<fdMatrix<T> >
);

template<class T>
std::shared_ptr<fdMatrix<T> > operator==
(
    const std::shared_ptr<fdMatrix<T> >,
    const std::shared_ptr<Field<T> >
);

template<class T>
std::shared_ptr<fdMatrix<T> > operator-
(
    const std::shared_ptr<fdMatrix<T> >,
    const std::shared_ptr<fdMatrix<T> >
);

template<class T>
std::shared_ptr<fdMatrix<T> > operator-
(
    const std::shared_ptr<fdMatrix<T> >
);

#include "im/fdMatrix/fdMatrix.hxx"

}

#endif
 
