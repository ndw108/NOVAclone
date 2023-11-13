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
#ifndef FIELD_H
#define FIELD_H

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <memory>
#include "Mesh/Mesh.h"
#include "Types/vector/vector.h"
#include "Types/scalar/scalar.h"
#include "Types/tensor/tensor.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <thread>
#include "boost/format.hpp"
#include "reuseTmp/reuseTmp.h"
#include "parallelCom/parallelCom.h"

#include <boost/align/aligned_alloc.hpp>
#include <boost/align/assume_aligned.hpp>

#ifdef HAVE_MPI
#include <boost/mpi.hpp>
#endif

#ifdef HAVE_MUI
#include "MUI/mui.h"
#endif



template <class T> class Field;
template <class T> class BC;

template <class T>
class Field
{
    protected:
    Mesh m_;
    std::string FieldName_;
    T* v_;
    std::unique_ptr<BC<T> > bc_[54];

    int curTime_;


    public:

    //constructors
    Field( Mesh&, std::string );
    Field( Mesh&, const T&, std::string );
    Field( Field<T>& );
    Field( std::shared_ptr<Field<T> > );
    Field( Mesh& );
    Field( Mesh&, const T& );

    //destructor
    ~Field();

    //access
    inline int index(int i, int j, int k) const
    {
        return k+j*nk()+i*(nk()*nj());
    }
    
    inline T& operator()(int i, int j, int k) const
    {
        return this->v_[this->index(i,j,k)];
    }

    inline T& operator()(int i) const
    {
        return this->v_[i];
    }


    T* ptr() const
    {
        return v_;
    }

    const Mesh& mesh() const
    {
        return m_;
    }
    Mesh& mesh()
    {
        return m_;
    }
    int ni() const
    {
        return m_.ni();
    }
    int nj() const
    {
        return m_.nj();
    }
    int nk() const
    {
        return m_.nk();
    }

    std::string& name()
    {
        return FieldName_;
    }

    void correctBoundaryConditions();  
    void write(std::string, std::string);
    void read(std::string, std::string);
    void read(std::string, std::string, int);
    void setBC( std::string );
 
    //manipulation

    virtual void operator+=( const Field<T>& ); 
    virtual void operator+=( const scalar );
    virtual void operator+=( const std::shared_ptr<Field<T> > ); 
    virtual void operator-=( const Field<T>& ); 
    virtual void operator-=( const scalar );
    virtual void operator-=( const std::shared_ptr<Field<T> > ); 
    virtual void operator*=( const Field<T>& ); 
    virtual void operator*=( const scalar ); 
    virtual void operator*=( const std::shared_ptr<Field<T> > );
    virtual void operator/=( const scalar ); 
    virtual void operator=( const Field<T>& ); 
    virtual void operator=( const scalar ); 
    virtual void operator=( const std::shared_ptr<Field<T> > ); 

    void updateMUI();  

    void bound( T, T );

    T weightedAverage();
};


template <class T>
std::shared_ptr<Field<T> > operator+
( 
    const Field<T>&,
    const Field<T>&
);

template <class T>
std::shared_ptr<Field<T> > operator+
( 
    const Field<T>&,
    const std::shared_ptr<Field<T> >
);

template <class T>
std::shared_ptr<Field<T> > operator-
( 
    const std::shared_ptr<Field<T> >
);

template <class T>
std::shared_ptr<Field<T> > operator-
( 
    const Field<T>&,
    const Field<T>&
);

template <class T>
std::shared_ptr<Field<T> > operator-
( 
    const Field<T>&,
    const T&
);

template <class T>
std::shared_ptr<Field<T> > operator-
( 
    const std::shared_ptr<Field<T> >,
    const T&
);



template <class T>
std::shared_ptr<Field<T> > operator-
( 
    const Field<T>&,
    const std::shared_ptr<Field<T> >
);


template <class T>
std::shared_ptr<Field<T> > operator-
( 
    const std::shared_ptr<Field<T> >,
    const std::shared_ptr<Field<T> >
);

template <class T>
std::shared_ptr<Field<T> > operator+
( 
    const std::shared_ptr<Field<T> >,
    const std::shared_ptr<Field<T> >
);


template <class T>
std::shared_ptr<Field<T> > operator*
( 
    scalar,
    const Field<T>&
);

template <class T>
std::shared_ptr<Field<T> > operator*
( 
    scalar,
    const std::shared_ptr<Field<T> >&
);


template <class T1, class T2>
std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > operator*
( 
    const Field<T1>&,
    const Field<T2>&
);

template <class T1, class T2>
std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > operator*
( 
    const Field<T1>&,
    const T2&
);

template <class T1, class T2>
std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > operator*
( 
    const std::shared_ptr<Field<T1> >&,
    const T2&
);

template <class T1, class T2>
std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > operator*
( 
    const std::shared_ptr<Field<T1> >&,
    const Field<T2>&
);

template <class T1, class T2>
std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > operator*
( 
    const Field<T1>&,
    const std::shared_ptr<Field<T2> >&
);



template <class T1, class T2>
std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > operator/
( 
    const Field<T1>&,
    const Field<T2>&
);

template <class T1, class T2>
std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > operator/
( 
    const std::shared_ptr<Field<T1> >,
    const Field<T2>&
);

template <class T1, class T2>
std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > operator/
( 
    const Field<T1>&,
    const T2&
);


template <class T>
std::shared_ptr<Field<T> > dev2
(
    std::shared_ptr<Field<T> >
);


template <class T> 
std::shared_ptr<Field<T> > transpose
(
    std::shared_ptr<Field<T> >
);

std::shared_ptr<Field<symmTensor> > twoSymm
(
    std::shared_ptr<Field<tensor> >
);

std::shared_ptr<Field<symmTensor> > symm
(
    std::shared_ptr<Field<tensor> >
);

std::shared_ptr<Field<scalar> > operator&&
(
    std::shared_ptr<Field<tensor> >,
    std::shared_ptr<Field<tensor> >
);

std::shared_ptr<Field<scalar> > sqrt
(
    std::shared_ptr<Field<scalar> >
);

std::shared_ptr<Field<scalar> > sqrt
(
    Field<scalar>&
);


std::shared_ptr<Field<vector> > operator^
(
    const Field<vector>&,
    const Field<vector>&
);


template <class T1, class T2>
std::shared_ptr<Field<typename innerProductType<T1, T2>::type> > operator&
( 
    const Field<T1>& f1,
    const Field<T2>& f2
);


template <class T1, class T2>
std::shared_ptr<Field<typename innerProductType<T1, T2>::type> > operator&
( 
    const Field<T1>& f1,
    const std::shared_ptr<Field<T2> >& f2
);

template <class T1, class T2>
std::shared_ptr<Field<typename innerProductType<T1, T2>::type> > operator&
( 
    const std::shared_ptr<Field<T1> >& f1,
    const std::shared_ptr<Field<T2> >& f2
);

#include "Field/Field.hxx"
#include "Field/updateMUI.hxx"

#endif
