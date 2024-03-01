{
//    std::array<std::shared_ptr<Field<vector> >, 4> dU; 
//    for( int i=0; i<4; i++ )
//    {   
        reuseTmp<vector> dUt( mesh );
	std::shared_ptr<Field<vector> > dU( dUt() );
//        dU[i]= dUt();
//    }  

    reuseTmp<vector> UOldt( mesh );
    std::shared_ptr<Field<vector> > UOld( UOldt() );

    *UOld = U;

//    std::shared_ptr<Field<tensor> > gradU( fdc::grad( U ) );
//    gradU->correctBoundaryConditions();

    for( int rk=0; rk<4; rk++ )
    {
        *dU =
        time.dt() *
        (
           -ex::div( U, U )
           +ex::laplacian( mu, U ) 
//           +mu*ex::div( ( dev2(transpose(gradU) ) ) )
           -( (rho0*T)*g )  
        );

        Ustar = UOld;
//        for( int j=0; j<=rk; j++ )
//        {   
//            Ustar = Ustar+settings::RK4[rk]*dU;
            Ustar = Ustar+RK4[rk]*dU;
//        }   
   
	Ustar.correctBoundaryConditions();
//        pEqn.rhs(ex::div(Ustar) / time.dt() / settings::RK4c()[rk]);
//        pEqn.solve();
        p.correctBoundaryConditions();
//        U = Ustar - settings::RK4c()[rk]*time.dt() * ex::grad(p); //(p_rgh?)

        U.correctBoundaryConditions();
   } 







}
