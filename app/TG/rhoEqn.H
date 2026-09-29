{
    reuseTmp<scalar> drhot( mesh );
    std::shared_ptr<Field<scalar> > drho( drhot() );
  
    reuseTmp<scalar> rhoct( mesh );
    std::shared_ptr<Field<scalar> > rhoc( rhoct() );

    reuseTmp<scalar> rhoOldt( mesh );
    std::shared_ptr<Field<scalar> > rhoOld( rhoOldt() );

    *rhoOld = rho;
    *rhoc = rho;

    for( int rk=0; rk<4; rk++ )
    {
        *drho =
        time.dt() *
        (
           -ex::div( rho*U )
        );

        //drho->correctBoundaryConditions();

        *rhoc = rhoc+RK4[rk]*drho;

        if( rk<3 )
        {
            rho = rhoOld + (int(rk/2)+1)/2.0 * drho;
        }
        else
        {
            rho=rhoc;
        }

        rho.correctBoundaryConditions();        
    }
}
