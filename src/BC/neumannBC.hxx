template <class T>
void neumannBC<T>::update
(
    Field<T>* fld_ptr
)
{
        if( this->dir()== bottom && parallelCom::k() == 0 ) 
        {   
            for( int i=0; i<fld_ptr->ni(); i++ )
            {   
                for( int j=0; j<fld_ptr->nj(); j++ )
                {   
                    int k=settings::m()/2; 

                    component(fld_ptr->operator()(i, j, k), this->comp()) = component(fld_ptr->operator()(i, j, k+1), this->comp())
                    -this->val()*fld_ptr->mesh().dz();

                    for( int l=0; l<settings::m()/2; l++ )
                    {
                        component(fld_ptr->operator()(i, j, k-(l+1)), this->comp()) = component(fld_ptr->operator()(i, j, k+(l+1)), this->comp())
                       -2.0*(l+1)*this->val()*fld_ptr->mesh().dz();
		    }
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
                
		    component(fld_ptr->operator()(i, j, k), this->comp()) = component(fld_ptr->operator()(i, j, k-1), this->comp()) 
                    + this->val()*fld_ptr->mesh().dz();

                    for( int l=0; l<settings::m()/2; l++ )
                    {   
                        component(fld_ptr->operator()(i, j, k+(l+1)), this->comp()) = component(fld_ptr->operator()(i, j, k-(l+1)), this->comp())
                       -2.0*(l+1)*this->val()*fld_ptr->mesh().dz();
                    } 
                }   
            }   
        }      
}
