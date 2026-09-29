{
    std::array<std::shared_ptr<Field<vector> >, 4> dB; 
    for( int i=0; i<4; i++ )
    {   
        reuseTmp<vector> dBt( mesh );
        dB[i]= dBt();
    }


    reuseTmp<vector> Boldt( mesh );
    std::shared_ptr<Field<vector> > Bold( Boldt() );

    *Bold = B;

    for( int rk=0; rk<4; rk++ )
    {

	dB[rk] =
        time.dt() *
        (
           -ex::div( U, B )
	   +ex::div( B, U ) 
           +ex::laplacian( eta, B ) 
        );
	
	Bstar = Bold;
	for( int j=0; j<=rk; j++ )
        {   
	      Bstar = Bstar + settings::RK4a()[rk][j]*dB[j]; 
	} 

	B = Bstar;
	B.correctBoundaryConditions();

    }

    pBEqn.rhs(ex::div(B));
    pBEqn.solve();
    pB.correctBoundaryConditions();
    B -= ex::grad(pB);
    B.correctBoundaryConditions();

    J = ex::curl(B);
    J.correctBoundaryConditions();

    if( time.writelog() )
    {   
       auto divB = ex::div(B);
       auto maxdivB = tools::maximum(divB);
        if( parallelCom::master() )
        {   
            std::cout<<"Magnetic field divergence: " << maxdivB<<std::endl;
        }   
    }

}
