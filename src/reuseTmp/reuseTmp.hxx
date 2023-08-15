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
template <class T>
std::vector<std::shared_ptr<Field<T> > > reuseTmp<T>::tmp_( 15, std::shared_ptr<Field<T> >( NULL ) );

template <class T>
int reuseTmp<T>::index_ = -1;

template <class T>
reuseTmp<T>::reuseTmp
(
    Mesh m
)
{
    #pragma omp master
    for( int i=0; i<=15; i++ )
    {
        assert( i!=15 );
        if( tmp_[i].unique() || tmp_[i] == NULL)
        {
            if( tmp_[i] == NULL )
            {
                tmp_[i] = std::shared_ptr<Field<T> >( new Field<T>( m, 0.0 ) );
            }
            index_ = i;
            break;
        }
    }

    memset( tmp_[index_]->ptr(), 0, sizeof(T)*m.ni()*m.nj()*m.nk() );

    #pragma omp barrier
}













