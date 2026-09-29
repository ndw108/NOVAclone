template <class T>
void muiBC<T>::update
(
    Field<T>* fld_ptr
)
{
#ifdef HAVE_MUI
    bool updated = ( this->curTime_ == fld_ptr->mesh().time().timeStep() ) ? true : false;
    //if( updated ) return;
    this->curTime_ = fld_ptr->mesh().time().timeStep();    

    std::stringstream muiName;
    if( neumann_ )
    {
        muiName << "gradz";
    }

    muiName << fld_ptr->name() << this->comp();

    if( settings::debug() )
    {
        std::cout<< "For field "<< fld_ptr->name()<< ", component "<<this->comp()<<", Interface "<<inf_<<", zone "<<settings::zoneName();
        if( parallelCom::muiInterface(inf_)->is_ready( muiName.str(), fld_ptr->mesh().time().timeStep() ) )
        {
            std::cout<< "is ready"<<std::endl;
        }
        else
        {
            std::cout<< "is NOT ready!"<<std::endl;
        }
    }

    scalar r = std::pow
    (
        ( 
            std::pow( fld_ptr->mesh().dx(), 2) 
          + std::pow( fld_ptr->mesh().dz(), 2) 
          + std::pow( fld_ptr->mesh().dz(), 2)
        ), 0.5
    );

    //mui::sampler_gauss<double,double,mui::default_config> s1( 2.0*r, 0.05*r ); 
    //mui::sampler_nearest_neighbor<double,double,mui::default_config> s1;
    mui::sampler_exact<double,double,mui::default_config> s1;
    mui::chrono_sampler_exact<> s2;

    if( this->dir() == bottom && parallelCom::k() == 0 )
    {
        for( int i=0; i<fld_ptr->ni(); i++ )
        {
            for( int j=0; j<fld_ptr->nj(); j++ )
            {
                int k = settings::m()/2;
                vector loc( fld_ptr->mesh().loc(i, j, k) );
                loc = loc + parallelCom::muiTrans(inf_);
                scalar val =  
                parallelCom::muiInterface(inf_)->fetch
                ( 
                    muiName.str(), 
                    mui::point3d( loc.x(), loc.y(), loc.z() ), 
                    fld_ptr->mesh().time().timeStep(),
                    s1,
                    s2
                ); 

                if( dirichlet_ )
                {
                    this->dirichlet( fld_ptr, i, j, k, val );
                }
                else
                {
                    scalar c = 
                    parallelCom::muiInterface(inf_)->fetch<double>( "alpha_f" ) /
                    parallelCom::muiInterface(inf_)->fetch<double>( "alpha_s" );

                    this->neumann( fld_ptr, i, j, k, val*c );
                }
            }
        }
    }

    if( this->dir() == top && parallelCom::k() == parallelCom::nk()-1 )
    {
        for( int i=0; i<fld_ptr->ni(); i++ )
        {
            for( int j=0; j<fld_ptr->nj(); j++ )
            {
                int k = fld_ptr->nk()-settings::m()/2-1; 
                vector loc( fld_ptr->mesh().loc(i, j, k) );
                loc = loc + parallelCom::muiTrans(inf_);
                scalar val = 
                parallelCom::muiInterface(inf_)->fetch
                ( 
                    muiName.str(), 
                    mui::point3d( loc.x(), loc.y(), loc.z() ), 
                    fld_ptr->mesh().time().timeStep(),
                    s1,
                    s2
                );

                if( dirichlet_ )
                {
                    this->dirichlet( fld_ptr, i, j, k, val );
                }
                else
                {
                    scalar c = 
                    parallelCom::muiInterface(inf_)->fetch<double>( "alpha_f" ) /
                    parallelCom::muiInterface(inf_)->fetch<double>( "alpha_s" );

                    this->neumann( fld_ptr, i, j, k, val*c );
                }
            }
        }
    }


#endif
}

