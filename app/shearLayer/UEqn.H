{
    reuseTmp<vector> drhoUt( mesh );
    std::shared_ptr<Field<vector> > drhoU( drhoUt() );

    reuseTmp<vector> rhoUOldt( mesh );
    std::shared_ptr<Field<vector> > rhoUOld( rhoUOldt() );

    reuseTmp<vector> rhoUct( mesh );
    std::shared_ptr<Field<vector> > rhoUc( rhoUct() );

    *rhoUc = rho*U;
    *rhoUOld = rho*U;

//    std::shared_ptr<Field<tensor> > gradU( fdc::grad( U ) );
//    gradU->correctBoundaryConditions();

    for( int rk=0; rk<4; rk++ )
    {
        *drhoU =
        time.dt() *
        (
           -fdc::div( rho*U, U )
           -fdc::grad( p )
           +fdc::laplacian( mu, U ) 
//           +mu*fdc::div( ( dev2(transpose(gradU) ) ) )
           -( (rho*T)*g )  
        );

        *rhoUc = rhoUc+RK4[rk]*drhoU;

        if( rk<3 )
        {
            U = rhoUOld/rho + (int(rk/2)+1)/2.0 * (drhoU/rho);
        }
        else
        {
            U=rhoUc/rho;
        }    

        U.correctBoundaryConditions();
    }
}
