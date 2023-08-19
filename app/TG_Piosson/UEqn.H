{
    reuseTmp<vector> dUt( mesh );
    std::shared_ptr<Field<vector> > dU( dUt() );

    reuseTmp<vector> UOldt( mesh );
    std::shared_ptr<Field<vector> > UOld( UOldt() );

    reuseTmp<vector> Uct( mesh );
    std::shared_ptr<Field<vector> > Uc( Uct() );

    *Uc = U;
    *UOld = U;

    for( int rk=0; rk<4; rk++ )
    {
        *dU =
        time.dt() *
        (
           -ex::div( U, U )
           +ex::laplacian( mu, U ) 
        );

        *Uc = Uc+RK4[rk]*dU;

        if( rk<3 )
        {
            U = UOld + (int(rk/2)+1)/2.0 * dU;
        }
        else
        {
            U=Uc;
        }    
    
        U.correctBoundaryConditions();
    }

}
