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
std::shared_ptr<fdMatrix<T> > ddt
(
    Field<T>& fld
)
{
    std::shared_ptr<fdMatrix<T> > matrix( new fdMatrix<T>( fld ) );

    double dt=0.0001;
    matrix->ap() = 1.0/dt;

    (*static_cast<Field<T>* >(matrix.get())) = fld/dt;

    return matrix;
}


}
