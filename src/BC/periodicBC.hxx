template <class T>
void periodicBC<T>::update
(
    Field<T>* fld_ptr
)
{ 
    if( this->dir() == east && parallelCom::ni() == 1 )
    {
        #pragma omp parallel for collapse(2)
        for( int k=0; k<fld_ptr->nk(); k++ )
        {
            for( int j=0; j<fld_ptr->nj(); j++ )
            {
                for( int i=0; i<settings::m()/2; i++ )
                { 
                    fld_ptr->operator()(i, j, k) = fld_ptr->operator()(i+fld_ptr->ni()-settings::m()-1, j, k);
                }
                for( int i=fld_ptr->ni()-settings::m()/2-1; i<fld_ptr->ni(); i++ )
                {
                    fld_ptr->operator()(i, j, k) = fld_ptr->operator()(i-fld_ptr->ni()+settings::m()+1, j, k);
                }           
            }
        }
    }


    if( this->dir() == south && parallelCom::nj() == 1 )
    {
        #pragma omp parallel for collapse(2)
        for( int k=0; k<fld_ptr->nk(); k++ )
        {
            for( int i=0; i<fld_ptr->ni(); i++ )
            {
                for( int j=0; j<settings::m()/2; j++ )
                { 
                    fld_ptr->operator()(i, j, k) = fld_ptr->operator()(i, j+fld_ptr->nj()-settings::m()-1, k);
                }           
                for( int j=fld_ptr->nj()-settings::m()/2-1; j<fld_ptr->nj(); j++ )
                {
                    fld_ptr->operator()(i, j, k) = fld_ptr->operator()(i, j-fld_ptr->nj()+settings::m()+1, k);
                }
            }
        }
    }


    if( this->dir() == bottom && parallelCom::nk() == 1 )
    {
        #pragma omp parallel for collapse(2)
        for( int i=0; i<fld_ptr->ni(); i++ )
        {
            for( int j=0; j<fld_ptr->nj(); j++ )
            {
                for( int k=0; k<settings::m()/2; k++ )
                { 
                    fld_ptr->operator()(i, j, k) = fld_ptr->operator()(i, j, k+fld_ptr->nk()-settings::m()-1);
                }
                for( int k=fld_ptr->nk()-settings::m()/2-1; k<fld_ptr->nk(); k++ )
                {
                    fld_ptr->operator()(i, j, k) = fld_ptr->operator()(i, j, k-fld_ptr->nk()+settings::m()+1);
                }
            }
        }
    }

}
