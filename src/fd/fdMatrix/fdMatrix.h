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

#include<vector>
#include "Types/scalar/scalar.h"
#include "pTraits/pTraits.h"
#include "Field/Field.h"
#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

namespace fd
{

template<class T>
class fdMatrix :
    public Field<T>
{
protected:
    std::vector<Field<double> > al_;
    Field<double> ap_;
    Field<T>& transportedFld_;
 
public:
    //constructor
    fdMatrix( Field<T>& );
    fdMatrix( const fdMatrix<T>& );
    fdMatrix( const std::shared_ptr<fdMatrix<T> >& ); 

    std::vector<Field<double> >& al()
    {
        return al_;
    }
    const std::vector<Field<double> >& al() const
    {
        return al_;
    }

    Field<double>& ap()
    {
        return ap_;
    }
 

    const Field<double>& ap() const
    {
        return ap_;
    }
   
    void solve();
    void deferredCorrection( const std::shared_ptr<Field<T> >& );
    void relax( const double& );
    T residual();
    void setReference( int, int, int, double );

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

#include "fd/fdMatrix/fdMatrix.hxx"

}

#endif
 
