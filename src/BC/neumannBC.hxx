template <class T>
void neumannBC<T>::update
(
    Field<T>* fld_ptr
)
{
    if( this->dir() == bottom && parallelCom::k() == 0 )
    {
        for( int i=0; i<fld_ptr->ni(); i++ )
        {   
            for( int j=0; j<fld_ptr->nj(); j++ )
            {  
                int k=settings::m()/2;
                fld_ptr->operator()(i,j,k) =  this->val();
            }
        }
    }
    else if( this->dir() == top && parallelCom::k() == parallelCom::nk()-1 )
    {
        for( int i=0; i<fld_ptr->ni(); i++ )
        {   
            for( int j=0; j<fld_ptr->nj(); j++ )
            {  
                int k=fld_ptr->nk()-settings::m()/2-1;
                fld_ptr->operator()(i,j,k) =  this->val(); 
            }
        }
    }
}

