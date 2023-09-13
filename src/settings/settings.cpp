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
#include "settings/settings.h"

bool settings::debug_ = false;
bool settings::restart_ = false;
std::string settings::zoneName_(".");
int settings::m_(6);
const std::array<std::array<scalar, 4>, 4 > settings::coef_ = 
{ 
    std::array<scalar, 4>{ { 1.0, 0.0, 0.0, 0.0 } }, //2nd
    std::array<scalar, 4>{ { 4.0/3.0, -1.0/3.0, 0.0, 0.0 } }, //4th       
    std::array<scalar, 4>{ { 6.0/4.0, -12.0/20.0, 1.0/10.0, 0.0 } }, //6th 
    std::array<scalar, 4>{ { 8.0/5.0, -4.0/5.0, 24.0/105.0, -8.0/280.0 } } //8th
};

const std::array<std::array<scalar, 7>, 6 > settings::coefBD_ = 
{ 
    std::array<scalar, 7>{ { -1.0,       1.0,  0.0,      0.0,       0.0,      0.0,      0.0     } }, //1st
    std::array<scalar, 7>{ { -3.0/2.0,   2.0, -1.0/2.0,  0.0,       0.0,      0.0,      0.0     } }, //2nd       
    std::array<scalar, 7>{ { -11.0/6.0,  3.0, -3.0/2.0,  1.0/3.0,   0.0,      0.0,      0.0     } }, //3rd
    std::array<scalar, 7>{ { -11.0/6.0,  3.0, -3.0/2.0,  1.0/3.0,   0.0,      0.0,      0.0     } }, //4th
    std::array<scalar, 7>{ { -11.0/6.0,  3.0, -3.0/2.0,  1.0/3.0,   0.0,      0.0,      0.0     } }, //5th
    std::array<scalar, 7>{ { -49.0/20.0, 6.0, -15.0/2.0, 20.0/3.0, -15.0/4.0, 6.0/5.0, -1.0/6.0 } } //6th
};


const std::array<std::array<scalar, 4>, 4 > settings::RK4a_ =
{
    std::array<scalar, 4>{ { 0.5, 0  , 0  , 0 } },
    std::array<scalar, 4>{ { 0  , 0.5, 0  , 0 } },
    std::array<scalar, 4>{ { 0  , 0  , 1.0, 0 } },
    std::array<scalar, 4>{ { 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 } }
};
const std::array<scalar, 4 > settings::RK4c_ = { 0.5, 0.5, 1.0, 1.0 };

void settings::process
(
    int argc, 
    char* argv[]
)
{
    for( int i=1; i<argc; i++ )
    {
        if( argv[i] == std::string("-restart") )
        {
            settings::restart() = true;
        }

        if( argv[i] == std::string("-zone") )
        {
            settings::zoneName() = std::string( argv[i+1] );
        }

        if( argv[i] == std::string("-m") )
        {
            settings::m() = atoi( argv[i+1] );
            assert( settings::m()%2 == 0 );
        }

    }
}
