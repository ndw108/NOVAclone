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

#include "ex/curl/curl.h"

namespace ex
{

std::shared_ptr<Field<vector> > curl
( 
    const Field<vector>& fld 
)
{
    reuseTmp<vector> curlfldt( fld.mesh() );
    std::shared_ptr<Field<vector> > curlfld( curlfldt() ); 
    
    reuseTmp<tensor> gradt( fld.mesh() );
    std::shared_ptr<Field<tensor> > grad( gradt() );
    
    (*grad) = ex::grad( fld );
    
    for( int i=settings::m()/2; i<fld.ni()-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<fld.nj()-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<fld.nk()-settings::m()/2; k++ )
            {
                curlfld->operator()(i, j, k) = 
                vector
                (
                    (*grad)(i, j, k).zy() - (*grad)(i, j, k).yz(),
                    (*grad)(i, j, k).xz() - (*grad)(i, j, k).zx(),
                    (*grad)(i, j, k).yx() - (*grad)(i, j, k).xy()
                );
            }
        }
    }

    return curlfld;

}

std::shared_ptr<Field<vector> > curl
( 
    const std::shared_ptr<Field<vector> > fld
)
{
    return ex::curl( (*fld) );
}


}
