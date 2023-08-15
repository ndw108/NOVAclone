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


#ifdef HAVE_MPI
template <class T>
std::array<std::vector<boost::mpi::request>, 6> mpiBC<T>::reqs_;

template <class T>
void mpiBC<T>::initUpdate
(
    Field<T>* fld_ptr
)
{
    //send all components together
    if( this->comp() != 0 )
    {
        return;
    }

    int procI = parallelCom::i();
    int procJ = parallelCom::j();
    int procK = parallelCom::k();

    int ni=fld_ptr->ni();
    int nj=fld_ptr->nj();
    int nk=fld_ptr->nk();

   
    if( this->dir() == east && parallelCom::ni()>1 )
    {
        //clear send buffers
        sbuf_.clear();

        //pack
        for( int j=0; j<fld_ptr->nj(); j++ )
        {
            for( int k=0; k<fld_ptr->nk(); k++ )
            {
                for( int i=fld_ptr->ni()-settings::m()-1; i<fld_ptr->ni()-settings::m()/2-1; i++ )
                {
                    sbuf_.push_back( fld_ptr->operator()(i, j, k) );
                }
            }
        }

        rbuf_.resize( nj*nk*(settings::m()/2+1) ); 

        if( procI != parallelCom::ni()-1 )
        {
            reqs_[west].push_back( parallelCom::world().isend( parallelCom::procNo( procI+1, procJ, procK ), 4, &sbuf_[0], nj*nk*settings::m()/2 ) );
        }
        else if( this->type() == periodic )
        {
            reqs_[west].push_back( parallelCom::world().isend( parallelCom::procNo( 0, procJ, procK ), 4, &sbuf_[0], nj*nk*settings::m()/2 ) );
        }

        if( procI != parallelCom::ni()-1 )
        {
            reqs_[east].push_back( parallelCom::world().irecv( parallelCom::procNo( procI+1, procJ, procK ), 5, &rbuf_[0], nj*nk*(settings::m()/2+1) ) );
        }
        else
        {
            reqs_[east].push_back( parallelCom::world().irecv( parallelCom::procNo( 0, procJ, procK ), 5, &rbuf_[0], nj*nk*(settings::m()/2+1) ) );
        }
    }
    else if( this->dir() == west && parallelCom::ni()>1 )
    {
        //clear send buffers
        sbuf_.clear();

        //pack
        for( int j=0; j<fld_ptr->nj(); j++ )
        {
            for( int k=0; k<fld_ptr->nk(); k++ )
            {
                for( int i=settings::m()/2; i<settings::m()+1; i++ )
                {
                    sbuf_.push_back( fld_ptr->operator()(i, j, k) );
                }
            }
        }

        rbuf_.resize( nj*nk*settings::m()/2 ); 

        if( procI != 0 )
        {
            reqs_[west].push_back( parallelCom::world().irecv( parallelCom::procNo( procI-1, procJ, procK ), 4, &rbuf_[0], nj*nk*settings::m()/2 ) );
        }
        else
        { 
            reqs_[west].push_back( parallelCom::world().irecv( parallelCom::procNo( parallelCom::ni()-1, procJ, procK ), 4, &rbuf_[0], nj*nk*settings::m()/2 ) );
        }

        if( procI != 0 )
        {
            reqs_[east].push_back( parallelCom::world().isend( parallelCom::procNo( procI-1, procJ, procK ), 5, &sbuf_[0], nj*nk*(settings::m()/2+1) ) );
        }
        else if( this->type() == periodic )
        {
            reqs_[east].push_back( parallelCom::world().isend( parallelCom::procNo( parallelCom::ni()-1, procJ, procK ), 5, &sbuf_[0], nj*nk*(settings::m()/2+1) ) );
        }
    }
    else if( this->dir()==top && parallelCom::nk()>1 )
    {      
        sbuf_.clear(); 
        
        //pack
        for( int i=0; i<ni; i++ )
        {
            for( int j=0; j<nj; j++ )
            {
                for( int k=nk-settings::m()-1; k<nk-settings::m()/2-1; k++ )
                {
                    sbuf_.push_back( fld_ptr->operator()(i, j, k) );
                }
            }
        }

        rbuf_.resize( ni*nj*(settings::m()/2+1) );

        if( procK != parallelCom::nk()-1 )
        {
            reqs_[bottom].push_back( parallelCom::world().isend( parallelCom::procNo( procI, procJ, procK+1 ), 0, &sbuf_[0], ni*nj*settings::m()/2 ) );
        }
        else if( this->type() == periodic )
        {
            reqs_[bottom].push_back( parallelCom::world().isend( parallelCom::procNo( procI, procJ, 0 ), 0, &sbuf_[0], ni*nj*settings::m()/2 ) );
        }

        if( procK != parallelCom::nk()-1 )
        {
            reqs_[top].push_back( parallelCom::world().irecv( parallelCom::procNo( procI, procJ, procK+1 ), 1, &rbuf_[0], ni*nj*(settings::m()/2+1) ) );
        }
        else
        {
            reqs_[top].push_back( parallelCom::world().irecv( parallelCom::procNo( procI, procJ, 0 ), 1, &rbuf_[0], ni*nj*(settings::m()/2+1) ) );
        }
    }
    else if( this->dir()==bottom && parallelCom::nk()>1 )
    {      
        sbuf_.clear(); 
 
        for( int i=0; i<ni; i++ )
        {
            for( int j=0; j<nj; j++ )
            {
                for( int k=settings::m()/2; k<settings::m()+1; k++ )
                {
                    sbuf_.push_back( fld_ptr->operator()(i, j, k) );
                }
            }
        } 

        rbuf_.resize( ni*nj*settings::m()/2 );

        if( procK != 0 )
        {
            reqs_[bottom].push_back( parallelCom::world().irecv( parallelCom::procNo( procI, procJ, procK-1 ), 0, &rbuf_[0], ni*nj*settings::m()/2 ) );
        }
        else
        { 
            reqs_[bottom].push_back( parallelCom::world().irecv( parallelCom::procNo( procI, procJ, parallelCom::nk()-1 ), 0, &rbuf_[0], ni*nj*settings::m()/2 ) );
        }

        if( procK != 0 )
        {
            reqs_[top].push_back( parallelCom::world().isend( parallelCom::procNo( procI, procJ, procK-1 ), 1, &sbuf_[0], ni*nj*(settings::m()/2+1) ) );
        }
        else if( this->type() == periodic )
        {
            reqs_[top].push_back( parallelCom::world().isend( parallelCom::procNo( procI, procJ, parallelCom::nk()-1 ), 1, &sbuf_[0], ni*nj*(settings::m()/2+1) ) );
        }
    }
    else if( this->dir() == north && parallelCom::nj()>1 )
    {
        sbuf_.clear();

        //pack
        for( int i=0; i<ni; i++ )
        {
            for( int k=0; k<nk; k++ )
            {
                for( int j=nj-settings::m()-1; j<nj-settings::m()/2-1; j++ )
                {
                    sbuf_.push_back( fld_ptr->operator()(i, j, k) );
                }
            }
        }

        rbuf_.resize( ni*nk*(settings::m()/2+1) );

        if( procJ != parallelCom::nj()-1 )
        {
            reqs_[south].push_back( parallelCom::world().isend( parallelCom::procNo( procI, procJ+1, procK ), 2, &sbuf_[0], ni*nk*settings::m()/2 ) );
        }
        else if ( this->type() == periodic) 
        {
            reqs_[south].push_back( parallelCom::world().isend( parallelCom::procNo( procI, 0, procK ), 2, &sbuf_[0], ni*nk*settings::m()/2 ) );
        }

        if( procJ != parallelCom::nj()-1 )
        {
            reqs_[north].push_back( parallelCom::world().irecv( parallelCom::procNo( procI, procJ+1, procK ), 3, &rbuf_[0], ni*nk*(settings::m()/2+1) ) );
        }
        else
        {
            reqs_[north].push_back( parallelCom::world().irecv( parallelCom::procNo( procI, 0, procK ), 3, &rbuf_[0], ni*nk*(settings::m()/2+1) ) );
        }
    }
    else if( this->dir() == south && parallelCom::nj()>1 )
    {
        sbuf_.clear();

        //pack

        for( int i=0; i<ni; i++ )
        {
            for( int k=0; k<nk; k++ )
            {
                for( int j=settings::m()/2; j<settings::m()+1; j++ )
                {
                    sbuf_.push_back( fld_ptr->operator()(i, j, k) );
                }
            }
        }

        rbuf_.resize( ni*nk*settings::m()/2 );

        if( procJ != 0 )
        {
            reqs_[south].push_back( parallelCom::world().irecv( parallelCom::procNo( procI, procJ-1, procK ), 2, &rbuf_[0], ni*nk*settings::m()/2 ) );
        }
        else
        { 
            reqs_[south].push_back( parallelCom::world().irecv( parallelCom::procNo( procI, parallelCom::nj()-1, procK ), 2, &rbuf_[0], ni*nk*settings::m()/2 ) );
        }

        if( procJ != 0 )
        {
            reqs_[north].push_back( parallelCom::world().isend( parallelCom::procNo( procI, procJ-1, procK ), 3, &sbuf_[0], ni*nk*(settings::m()/2+1) ) );
        }
        else if( this->type() == periodic)
        {
            reqs_[north].push_back( parallelCom::world().isend( parallelCom::procNo( procI, parallelCom::nj()-1, procK ), 3, &sbuf_[0], ni*nk*(settings::m()/2+1) ) );
        }
    }
}


template <class T>
void mpiBC<T>::update
(
    Field<T>* fld_ptr
)
{
    //we'll send all components together
    if( this->comp() != 0 )
    {
        return;
    }

    int ni=fld_ptr->ni();
    int nj=fld_ptr->nj();
    int nk=fld_ptr->nk();

    boost::mpi::wait_all( reqs_[this->dir()].begin(), reqs_[this->dir()].end() );
    reqs_[this->dir()].clear();


    //unpack
    if( this->dir() == west && parallelCom::ni()>1 )
    {
        int loc=0; 
        for( int j=0; j<nj; j++ )
        {
            for( int k=0; k<nk; k++ )
            {
                for( int i=0; i<settings::m()/2; i++ )
                {
                    fld_ptr->operator()(i, j, k) = rbuf_[loc++];
                }
            }
        }
    }
    else if( this->dir() == east && parallelCom::ni()>1 )
    {
        int loc=0;
        for( int j=0; j<nj; j++ )
        {
            for( int k=0; k<nk; k++ )
            {
                for( int i=ni-settings::m()/2-1; i<ni; i++ )
                {
                    fld_ptr->operator()(i, j, k) = rbuf_[loc++];
                }
            }
        } 
    }
    else if( this->dir() == south && parallelCom::nj()>1 )
    {
        int loc=0;
        for( int i=0; i<ni; i++ )
        {
            for( int k=0; k<nk; k++ )
            {
                for( int j=0; j<settings::m()/2; j++ )
                { 
                    fld_ptr->operator()(i, j, k) = rbuf_[loc++];
                }
            }
        } 
    }
    else if( this->dir() == north && parallelCom::nj()>1 )
    {
        int loc=0;
        for( int i=0; i<ni; i++ )
        {
            for( int k=0; k<nk; k++ )
            {
                for( int j=nj-settings::m()/2-1; j<nj; j++ )
                {
                    fld_ptr->operator()(i, j, k) = rbuf_[loc++];
                }
            }
        } 
    }
    else if( this->dir() == bottom && parallelCom::nk()>1 )
    {
        int loc=0;
        for( int i=0; i<ni; i++ )
        {
            for( int j=0; j<nj; j++ )
            {
                for( int k=0; k<settings::m()/2; k++ )
                { 
                    fld_ptr->operator()(i, j, k) = rbuf_[loc++];
                }
            }
        }
    }
    else if( this->dir() == top && parallelCom::nk()>1 )
    {
        int loc=0;
        for( int i=0; i<ni; i++ )
        {
            for( int j=0; j<nj; j++ )
            {
                for( int k=nk-settings::m()/2-1; k<nk; k++ )
                {
                    fld_ptr->operator()(i, j, k) = rbuf_[loc++];
                }
            }
        }
    }
}
#endif
