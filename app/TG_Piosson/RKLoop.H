{

    std::array<std::shared_ptr<Field<vector> >, 4> dU;
    for( int i=0; i<4; i++ )
    {
        reuseTmp<vector> dUt( mesh );
        dU[i]= dUt();
    }

    reuseTmp<vector> UOldt( mesh );
    std::shared_ptr<Field<vector> > UOld( UOldt() );

    *UOld = U;

    for( int rk=0; rk<4; rk++ )
    {
        dU[rk] =
        time.dt() *
        (
           -ex::div( U, U )
           +ex::laplacian( mu, U ) 
        );

        Ustar = UOld;
        for( int j=0; j<=rk; j++ )
        {
            Ustar = Ustar+settings::RK4a()[rk][j]*dU[j];
        }
  }
        Ustar.correctBoundaryConditions();

        pEqn.rhs(ex::div(Ustar) / time.dt() );

	pEqn.solve();
        p.correctBoundaryConditions();
        U = Ustar - time.dt() * ex::grad(p);

        U.correctBoundaryConditions();
  
  

    if( time.writelog() )
    {
        auto divU = ex::div(U);
        auto maxdivU = tools::maximum(divU);

        if( parallelCom::master() )
        {
            std::cout<<"Velocity divergence: " << maxdivU<<std::endl;
        }
    }


}
