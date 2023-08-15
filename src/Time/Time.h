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
#ifndef TIME_H
#define TIME_H

#include <vector>
#include <fstream>
#include <sstream>
#include "Types/scalar/scalar.h"
#include "settings/settings.h"
#include "parallelCom/parallelCom.h"

class Time
{

  protected:
    std::vector<scalar> steps_;
    scalar dt_;
    scalar curTime_;
    scalar endTime_;
    int writeEvery_;
    int timeStep_;
    int fileNum_;
    bool writeEvery_timebased_;
    scalar writeInterval_;
    bool write_;
    bool log_;

  public:
    Time( scalar, scalar, int );
    Time( scalar, scalar, scalar );
    Time( std::string );
   
    scalar& dt()
    {
        return dt_;
    }

    const std::vector<scalar>& steps()
    {
        return steps_;
    }

    bool run() const
    {
        return curTime_<endTime_; 
    }

    bool write() const
    {
        return write_;
    }

    bool writelog() const
    {
        return log_;
    }


    scalar curTime() const
    {
        return curTime_;
    }

    int timeStep() const
    {
        return timeStep_;
    }

    int fileNum() const
    {
        return fileNum_;
    }


    void operator++(int);
};

#endif
