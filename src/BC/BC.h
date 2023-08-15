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
#ifndef BC_H
#define BC_H
#include "BCenums.h"
#include "parallelCom/parallelCom.h"

#ifdef HAVE_MUI
#include "MUI/mui.h"
#endif

#include "shapes/box.h"
#include "Field/Field.h"

template <class T> class Field;
template <class T> class BC;

template <class T>
class BC
{

    BC() {}
    protected:
    int curTime_;
    int dir_, comp_, type_;
    scalar val_;

    public:
    BC( std::string& l ) : curTime_(-1) 
    {
        std::stringstream strProcessor(l);
        std::string dir, type;

        strProcessor>>dir>>type>>val_>>comp_;

        dir_ = BCdirmap[dir];
        type_ = BCcodemap[type];
    }

    //access
    scalar& val()
    {
        return val_;
    }

    int dir() const
    {
        return dir_;
    }

    int comp() const
    {
        return comp_;
    }

    int type() const
    {
        return type_;
    }

    virtual void initUpdate( Field<T>* f ) {}

    virtual void update( Field<T>* ) = 0;

    virtual ~BC() = default;

private:
    void updateMUI( Field<T>* f );
};

template <class T>
class periodicBC : public BC<T>
{
    public:
    periodicBC( std::string& l ) : BC<T>( l ) {} 
    void update( Field<T>* );
};

template <class T>
class neumannBC : public BC<T>
{
    public:
    neumannBC( std::string& l ) : BC<T>( l ) {} 
    void update( Field<T>* );
};

template <class T>
class dirichletBC : public BC<T>
{
    public:
    dirichletBC( std::string& l ) : BC<T>( l ) {} 
    void update( Field<T>* );
};


template <class T>
class muiBC : public BC<T>
{
    int inf_;
    bool neumann_;
    bool dirichlet_;

    public:


    muiBC( std::string& l ) : BC<T>( l )
    {
        std::stringstream ss( l );
        std::string buf;
        ss>>buf>>buf>>buf>>buf>>inf_>>std::boolalpha>>neumann_>>std::boolalpha >>dirichlet_;

        if( neumann_ )
        {
            assert( !dirichlet_ );
        }

        if( dirichlet_ )
        {
            assert( !neumann_ );
        }
    }

    void update( Field<T>* );

    int inf() const
    {
        return inf_;
    }
};

template <class T>
class mpiBC : public BC<T>
{
    static std::array<std::vector<boost::mpi::request>, 6> reqs_;
    std::vector<T> rbuf_;
    std::vector<T> sbuf_;
    
public:
    mpiBC( std::string& l ) : BC<T>( l ) {}
    void initUpdate( Field<T>* );
    void update( Field<T>* );
};



template <class T1, class T2>
std::unique_ptr<BC<T1> > make(std::string l)
{
    return std::unique_ptr<BC<T1> >( new T2( l ) );
}

template <class T>
class BCmaker
{
  
public:
    static std::unique_ptr<BC<T> > (*makeBC[])(std::string);
};

template <class T>
std::unique_ptr<BC<T> > (* BCmaker<T>::makeBC [])(std::string) = 
{ 
    make<T, periodicBC<T> >, 
    make<T, neumannBC<T> >, 
    make<T, dirichletBC<T> >, 
    make<T, muiBC<T> >,
    make<T, mpiBC<T> >
};

#include "BC/periodicBC.hxx"
#include "BC/dirichletBC.hxx"
#include "BC/neumannBC.hxx"
#include "BC/muiBC.hxx"
#include "BC/mpiBC.hxx"

#endif
