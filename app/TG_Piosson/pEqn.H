{    
    pEqn.rhs(ex::div(U) / time.dt());
    pEqn.solve();
    p.correctBoundaryConditions();
    U -= time.dt() * ex::grad(p);
    U.correctBoundaryConditions();


    if( time.writelog() )
    {
        reuseTmp<scalar> divUt( mesh );
        auto divU( divUt() );
        divU = ex::div(U);
        auto maxdivU = tools::maximum(divU);

        if( parallelCom::master() )
        {
            std::cout<<"Velocity divergence: " << maxdivU<<std::endl;
        }
    }
}
