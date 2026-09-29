template <class T>
void Field<T>::updateMUI()
{
#ifdef HAVE_MUI
    if( parallelCom::numInterfaces() == 0 || this->name() != "T" )
    {
        return;
    }


    bool updated = ( curTime_ == mesh().time().timeStep() );

    if( !updated )
    {    
        curTime_ = mesh().time().timeStep();
        for( int inf=0; inf<parallelCom::numInterfaces(); inf++ )
        {
            parallelCom::muiInterface(inf)->forget(mesh().time().timeStep());  

            for( int comp=0; comp<pTraits<T>::nComp; comp++ )
            {
                std::stringstream muiName;
                muiName << this->name() << comp;
        
                int N=0;
    
                for( int i=0; i<this->ni(); i++ )
                {
                    for( int j=0; j<this->nj(); j++ )
                    {
                        for( int k=0; k<this->nk(); k++ )
                        {
                            vector loc( this->mesh().loc(i, j, k) );                       
                            vector locTrans( loc+parallelCom::muiTrans(inf) );
     
                            bool inside=false;
                            for( unsigned s=0; s<parallelCom::muiShapesSend(inf).size(); s++ )
                            {
                                if( parallelCom::muiShapesSend(inf)[s]->inside( loc ) )
                                {
                                    inside = true;
                                    break;
                                }
                            }
        
                            if( inside )
                            {
                                N++;
                                parallelCom::muiInterface(inf)->push
                                ( 
                                    muiName.str(),
                                    component( this->operator()(i, j, k), comp ) 
                                );

            
                                scalar gradz = 0;
                                for( int ll=0; ll<settings::m()/2-1; ll++ )
                                {
                                    int order = settings::m();
                                    gradz +=
                
                                    settings::coef(order)[ll] / (ll+1) / 2.0 *
                                    (
                                        component(this->operator()(i, j, k+ll+1), comp)
                                       -component(this->operator()(i, j, k-ll-1), comp)
                                    )/ this->mesh().dz() / this->mesh().dzdxi(k);
                                }
                            
                                parallelCom::muiInterface(inf)->push
                                ( 
                                    "gradz"+muiName.str(),
                                    gradz 
                                );
                            }
                        }
                    }
                }
 

                if( settings::debug() )
                {
                    std::cout<<"For interface "<<inf<<" , zone "<<settings::zoneName()<<", component "<<comp <<
                        " of field "<<this->name()<<" "<<N<<" points were pushed for rank "<< parallelCom::myProcNo()<<std::endl;
                }
            }
        
       
            parallelCom::muiInterface(inf)->commit(mesh().time().timeStep()); 
            parallelCom::muiInterface(inf)->barrier(mesh().time().timeStep()-1);
            parallelCom::muiInterface(inf)->forget(mesh().time().timeStep()-1);  
        }
    }
#endif
}

