{
    reuseTmp<vector> dUt( mesh );
    std::shared_ptr<Field<vector> > dU( dUt() );

    reuseTmp<vector> UOldt( mesh );
    std::shared_ptr<Field<vector> > UOld( UOldt() );

    reuseTmp<vector> Uct( mesh );
    std::shared_ptr<Field<vector> > Uc( Uct() );

    *UOld = U;


    for( int rk=0; rk<4; rk++ )
    {
        *dU =
        time.dt() *
        (
           -ex::div( U, U )
           +ex::laplacian( mu, U ) 
//           +mu*ex::div( ( dev2(transpose(gradU) ) ) )
//           -( (rho0*T)*(g) )  
           +ex::div( B, B ) 
           -ex::grad( 0.5*(B&B) )
        );


        Ustar = UOld;
        Ustar = Ustar+RK4[rk]*dU;
   
	Ustar.correctBoundaryConditions();
        pEqn.rhs(ex::div(Ustar) / time.dt() / settings::RK4c()[rk]);
        pEqn.solve();
        p.correctBoundaryConditions();
        U = Ustar - settings::RK4c()[rk]*time.dt() * ex::grad(p); 

        U.correctBoundaryConditions();
   } 

//    if( time.writelog() )
//    {   
        auto divU = ex::div(U);
        auto maxdivU = tools::maximum(divU);

        if( parallelCom::master() )
        {   
            std::cout<<"Velocity divergence: " << maxdivU<<std::endl;
        }   
//    }


}
