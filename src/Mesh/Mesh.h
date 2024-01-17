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
#ifndef MESH_H
#define MESH_H

#include "Time/Time.h"
#include "Types/vector/vector.h"
#include "parallelCom/parallelCom.h"
#include "boost/format.hpp"
#include "shapes/box.h"

#include <fstream>
#include <iostream>
#include <iomanip>


class Time;

class Mesh
{
    int ni_;
    int nj_;
    int nk_;

    ptrdiff_t glob_n_[3];

    ptrdiff_t cx;
    ptrdiff_t cy;
    ptrdiff_t cz;

    scalar dx_;
    scalar dy_;
    scalar dz_; 

    scalar lx_;
    scalar ly_;
    scalar lz_;

    scalar* dzdxi_;
    scalar* d2zdxi2_;

    vector origin_;

    Time& runTime_;

    Box bbox_;
    Box bboxHalo_;
    
    scalar s1_;
    scalar s2_;

    bool hyperbollic_;

#ifdef HAVE_FFTMPI
    ptrdiff_t localStart_[3]; 
#endif

    public:
    Mesh( int, int, int, scalar, scalar, scalar, Time& );
    Mesh( std::string, Time& );

    //access
    int ni() const
    {
        return ni_;
    }
    int nj() const
    {
        return nj_;
    }
    int nk() const
    {
        return nk_;
    }
    scalar dx() const
    {
        return dx_;
    }
    scalar dy() const
    {
        return dy_;
    }
    scalar dz() const
    {
        return dz_;
    }
    scalar lx() const
    {
        return lx_;
    }
    scalar ly() const
    {
        return ly_;
    }
    scalar lz() const
    {
        return lz_;
    }


    scalar z( int k ) const;

    scalar dzdxi( int k ) const
    {
        return dzdxi_[k];
    }

    scalar d2zdxi2( int k ) const
    {
        return d2zdxi2_[k];
    }

    Time& time() const
    {
        return runTime_;
    }

    vector origin() const
    {
        return origin_;
    }

    vector loc( int i, int j, int k ) const
    {
        return vector
        ( 
            origin_.x() + (i-settings::m()/2)*dx(), 
            origin_.y() + (j-settings::m()/2)*dy(), 
            this->z(k)
        );
    }

    Box& bbox()
    {
        return bbox_;
    } 
 
    Box& bboxHalo()
    {
        return bboxHalo_;
    } 

#ifdef HAVE_FFTMPI
    ptrdiff_t* glob_n()
    {
        return glob_n_;
    }
    ptrdiff_t* localStart()
    {
        return localStart_;
    }
#endif

    void sendMUIPoints();
    void write( std::string );
    void writeCase( std::string );
    void writeSOS();
    void writeRestartInfo( std::string );
};

#endif
