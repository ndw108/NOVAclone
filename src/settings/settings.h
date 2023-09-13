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
#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
#include <cassert>
#include <vector>
#include <array>
#include "Types/scalar/scalar.h"

class settings
{
    static bool debug_;
    static bool restart_;
    static std::string zoneName_;
    static int m_; 
    static const std::array< std::array<scalar, 4>, 4 > coef_;
    static const std::array< std::array<scalar, 7>, 6 > coefBD_;
    static const std::array<std::array<scalar, 4 >, 4 > RK4a_;
    static const std::array<scalar, 4 > RK4c_;

    public:
    static bool& debug()
    {
        return debug_;
    }

    static bool& restart()
    {
        return restart_;
    }

    static std::string& zoneName()
    {
        return zoneName_;
    }

    inline static int& m()
    {
        return m_;
    }

    inline static const std::array<scalar, 4>& coef()
    {
        return coef_[(m_-2)/2];
    }

    inline static const std::array<scalar, 4>& coef(int o)
    {
        return coef_[(o-2)/2];
    }

    inline static const std::array<scalar, 7>& coefBD(int o)
    {
        return coefBD_[o-1];
    }
    
    inline static const std::array<std::array<scalar, 4>, 4>& RK4a()
    {
        return RK4a_;
    }

    inline static const std::array<scalar, 4>& RK4c()
    {
        return RK4c_;
    }


    static void process( int argc, char* argv[] );
};

#endif
