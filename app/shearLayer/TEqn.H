{
    reuseTmp<scalar> dTt( mesh );
    std::shared_ptr<Field<scalar> > dT( dTt() );
  
    reuseTmp<scalar> Tct( mesh );
    std::shared_ptr<Field<scalar> > Tc( Tct() );

    reuseTmp<scalar> TOldt( mesh );
    std::shared_ptr<Field<scalar> > TOld( TOldt() );

    *Tc = T;
    *TOld = T;

    for( int rk=0; rk<4; rk++ )
    {
        *dT =
        time.dt() *
        (
           -ex::div( U, T )
           +ex::laplacian( ka, T )
        );

//        drhoT->correctBoundaryConditions();

        *Tc = Tc+RK4[rk]*dT;

        if( rk<3 )
        {
            T = TOld + (int(rk/2)+1)/2.0 * (dT);
        }
        else
        {
            T=Tc;
        }
    
        T.correctBoundaryConditions();
    }
}
