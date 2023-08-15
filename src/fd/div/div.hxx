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
namespace fd
{

template<class T>
std::shared_ptr<fdMatrix<T> > div
( 
    const Field<vector>& fld1, 
    Field<T>& fld2 
)
{
    std::shared_ptr<fdMatrix<T> > matrix( new fdMatrix<T>( fld2 ) );

    //comput CDS coefs
    for( int i=settings::m()/2; i<fld1.ni()-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<fld1.nj()-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<fld1.nk()-settings::m()/2; k++ )
            {
                matrix->al()[0]()[i][j][k] =  0.5 / fld1.mesh().dy() * fld1()[i][j+1][k].y(); //N
                matrix->al()[1]()[i][j][k] = -0.5 / fld1.mesh().dy() * fld1()[i][j-1][k].y(); //S
                matrix->al()[2]()[i][j][k] =  0.5 / fld1.mesh().dx() * fld1()[i+1][j][k].x(); //E
                matrix->al()[3]()[i][j][k] = -0.5 / fld1.mesh().dx() * fld1()[i-1][j][k].x(); //W
                matrix->al()[4]()[i][j][k] =  0.5 / fld1.mesh().dz() * fld1()[i][j][k+1].z(); //T
                matrix->al()[5]()[i][j][k] = -0.5 / fld1.mesh().dz() * fld1()[i][j][k-1].z(); //B 
            }
        }
    }

//    //apply deferred corrention
    std::shared_ptr<Field<T> > hoa = fdc::div( fld1, fld2 );
    matrix->deferredCorrection( hoa );   
    
    return matrix; 
}

}
