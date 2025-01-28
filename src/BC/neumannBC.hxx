template <class T>
void neumannBC<T>::update
(
    Field<T>* fld_ptr
)
{
	//bottom
        if( this->dir() == bottom && parallelCom::k() == 0 ) 
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

	//top
        if( this->dir() == top && parallelCom::k() == parallelCom::nk()-1 )
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


	//east
        if( this->dir() == east && parallelCom::i() == 0 ) 
        {   
            for( int j=0; j<fld_ptr->nj(); j++ )
            {   
                for( int k=0; k<fld_ptr->nk(); k++ )
                {   
                    int i=settings::m()/2; 

                    component(fld_ptr->operator()(i, j, k), this->comp()) = component(fld_ptr->operator()(i+1, j, k), this->comp())
                    -this->val()*fld_ptr->mesh().dx();

                    for( int l=0; l<settings::m()/2; l++ )
                    {   
                        component(fld_ptr->operator()(i-(l+1), j, k), this->comp()) = component(fld_ptr->operator()(i+(l+1), j, k), this->comp())
                       -2.0*(l+1)*this->val()*fld_ptr->mesh().dx();
                    }   
                }   
            }   
        }

        //west
        if( this->dir() == west && parallelCom::i() == parallelCom::ni()-1 )
        {   
            for( int j=0; j<fld_ptr->nj(); j++ )
            {   
                for( int k=0; k<fld_ptr->nk(); k++ )
                {   
                    int i=fld_ptr->ni()-settings::m()/2-1;
    
                    component(fld_ptr->operator()(i, j, k), this->comp()) = component(fld_ptr->operator()(i-1, j, k), this->comp()) 
                    + this->val()*fld_ptr->mesh().dx();

                    for( int l=0; l<settings::m()/2; l++ )
                    {   
                        component(fld_ptr->operator()(i+(l+1), j, k), this->comp()) = component(fld_ptr->operator()(i-(l+1), j, k), this->comp())
                       -2.0*(l+1)*this->val()*fld_ptr->mesh().dx();
                    }   
                }   
            }   
        }

        //south
        if( this->dir() == south && parallelCom::j() == 0 ) 
        {   
            for( int i=0; i<fld_ptr->ni(); i++ )
            {   
                for( int k=0; k<fld_ptr->nk(); k++ )
                {   
                    int j=settings::m()/2; 

                    component(fld_ptr->operator()(i, j, k), this->comp()) = component(fld_ptr->operator()(i, j+1, k), this->comp())
                    -this->val()*fld_ptr->mesh().dy();

                    for( int l=0; l<settings::m()/2; l++ )
                    {
                        component(fld_ptr->operator()(i, j-(l+1), k), this->comp()) = component(fld_ptr->operator()(i, j+(l+1), k), this->comp())
                       -2.0*(l+1)*this->val()*fld_ptr->mesh().dy();
                    }
                }
            }
        }

	//north
        if( this->dir() == north && parallelCom::j() == parallelCom::nj()-1 )
        {
            for( int i=0; i<fld_ptr->ni(); i++ )
            {
                for( int k=0; k<fld_ptr->nk(); k++ )
                {
                    int j=fld_ptr->nj()-settings::m()/2-1;

                    component(fld_ptr->operator()(i, j, k), this->comp()) = component(fld_ptr->operator()(i, j-1, k), this->comp())
                    + this->val()*fld_ptr->mesh().dy();

                    for( int l=0; l<settings::m()/2; l++ )
                    {
                        component(fld_ptr->operator()(i, j+(l+1), k), this->comp()) = component(fld_ptr->operator()(i, j-(l+1), k), this->comp())
                       -2.0*(l+1)*this->val()*fld_ptr->mesh().dy();
                    }
                }
            }
        }

}
