{
    reuseTmp<scalar> drhoTt( mesh );
    std::shared_ptr<Field<scalar> > drhoT( drhoTt() );
  
    reuseTmp<scalar> rhoTct( mesh );
    std::shared_ptr<Field<scalar> > rhoTc( rhoTct() );

    reuseTmp<scalar> rhoTOldt( mesh );
    std::shared_ptr<Field<scalar> > rhoTOld( rhoTOldt() );

    *rhoTc = rho*T;
    *rhoTOld = rho*T;

    for( int rk=0; rk<4; rk++ )
    {
        *drhoT =
        time.dt() *
        (
           -ex::div( rho*U, T )
           +ex::laplacian( ka, T )
        );

//        drhoT->correctBoundaryConditions();

        *rhoTc = rhoTc+RK4[rk]*drhoT;

        if( rk<3 )
        {
            T = rhoTOld/rho + (int(rk/2)+1)/2.0 * (drhoT/rho);
        }
        else
        {
            T=rhoTc/rho;
        }
    
        T.correctBoundaryConditions();
    }
}
