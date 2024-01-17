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
#ifndef PARALLELCOM_H
#define PARALLELCOM_H

#ifdef HAVE_MPI
#include <boost/mpi.hpp>
#endif

#include <functional>
#include <vector>
#include <cmath> 
#include <iostream>

#include "Types/vector/vector.h"
#include "settings/settings.h"
#include "shapes/shape.h"
#include "shapes/box.h"
#include "Mesh/Mesh.h"
#include "BC/BCenums.h"

#ifdef HAVE_MUI
#include "MUI/mui.h"
#endif

#ifdef HAVE_FFTMPI
#include "fft3d.h"
#endif

template <class T>
class Field;

class parallelCom
{
    static int myRank_;
    static int worldSize_;
    static int ni_, nj_, nk_;
    static int i_, j_, k_;


#ifdef HAVE_MPI
    static boost::mpi::environment env_;
    static boost::mpi::communicator world_;
#endif

#ifdef HAVE_MUI
    static std::vector<mui::uniface<mui::default_config>*> muiInterfaces_;
    static std::array<std::vector<std::unique_ptr<shape> >, 10> muiShapesSend_;
    static std::array<std::vector<std::unique_ptr<shape> >, 10> muiShapesRecv_;
    static std::array<vector, 10> muiTrans_;
#endif



private:
    parallelCom() {};

public:
    //access
    static int myProcNo()
    {
        return myRank_;
    }
    
    static int worldSize()
    {
        return worldSize_;
    }

    static bool master()
    {
        return myProcNo() == 0;
    }

    static int ni()
    {
        return ni_;
    }

    static int nj()
    {
        return nj_;
    }

    static int nk()
    {
        return nk_;
    }

    static int i()
    {
        return i_;
    }

    static int j()
    {
        return j_;
    }

    static int k()
    {
        return k_;
    }

    static int procNo( int i, int j, int k )
    {
        return k+j*nk()+i*(nk()*nj());
    }

    static void decompose(int, int, int);
    static void decompose( std::string );

    static bool internal( int dir, int type )
    {

        if( type == periodic )
        {
            if( ( dir==north || dir==south ) && nj_>1 )
            {
                return true;
            }
            if( ( dir==east || dir==west ) && ni_>1 )
            {
                return true;
            }
            if( ( dir==top || dir==bottom ) && nk_>1 )
            {
                return true;
            }
        }


        if( dir == north && j_==nj_-1 )
        {
            return false;
        }
        else if( dir == south && j_==0 )
        {
            return false;
        }
        else if( dir == east && i_==ni_-1 )
        {
            return false;
        }
        else if( dir == west && i_==0 )
        {
            return false;
        }
        else if( dir == top && k_==nk_-1 )
        {
            return false;
        }
        else if( dir == bottom && k_==0 )
        {
            return false;
        }

        return true;
    }

#ifdef HAVE_MPI
    static boost::mpi::communicator& world()
    {
        return world_;
    }
#endif

#ifdef HAVE_MUI
    static mui::uniface<mui::default_config>* muiInterface(int i)
    {
        return muiInterfaces_[i];
    }

    static std::vector<std::unique_ptr<shape> >& muiShapesSend( int i )
    {
        return muiShapesSend_[i];
    }

    static std::vector<std::unique_ptr<shape> >& muiShapesRecv( int i )
    {
        return muiShapesRecv_[i];
    }

    static vector& muiTrans( int i )
    {
        return muiTrans_[i];
    }

    static int numInterfaces()
    {
        return muiInterfaces_.size();
    }
    
    static void initMUI( std::string, Box );  

#endif

};

template <class T>
class plusOp
{
public:
    T operator()(T a, T b) {
       return a+b;
    }
};

template <class T>
class maxOp
{
public:
    T operator()(T a, T b) {
       return std::max(a,b);
    }
};

template <class T>
class minOp
{
public:
    T operator()(T a, T b) {
       return std::min(a,b);
    }
};


template <class T>
void reduce( T& val, auto f )
{
#ifdef HAVE_MPI
    T sumVal;
    boost::mpi::reduce(parallelCom::world(), val, sumVal, f, 0);
    val = sumVal;
#endif
}

template <class T>
void all_reduce( T& val, auto f, boost::mpi::communicator& comm=parallelCom::world() )
{
#ifdef HAVE_MPI
    T sumVal;
    boost::mpi::all_reduce(comm, val, sumVal, f);
    val = sumVal;
#endif
}

template <class T>
void reduce( std::vector<T>& vals, std::plus<T> )
{
#ifdef HAVE_MPI
    std::vector<T> res(vals.size(), T(0) );
    boost::mpi::reduce(parallelCom::world(), &vals[0], vals.size(), &res[0], std::plus<T>(), 0);
    vals = res;
#endif
}

template <class T>
void all_reduce( std::vector<T>& vals, std::plus<T> )
{
#ifdef HAVE_MPI
    std::vector<T> res(vals.size(), T(0) );
    boost::mpi::all_reduce(parallelCom::world(), &vals[0], vals.size(), &res[0], std::plus<T>());
    vals = res;
#endif
}

template <class T>
void broadcast( T& val )
{
#ifdef HAVE_MPI
    boost::mpi::broadcast(parallelCom::world(), val, 0);
#endif
}
#endif
