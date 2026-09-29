//scalar for Helmholtz–Hodge div cleaning is chi
//
{
    auto divU = ex::div(U);
    auto maxdivU = tools::maximum(divU);

    if( parallelCom::master() )
    {   
        std::cout<<"Max velocity divergence pre-HH: " << maxdivU<<std::endl;
    }

    Ustar = U;

    Ustar.correctBoundaryConditions();
    chiEqn.rhs(ex::div(Ustar) );
    chiEqn.solve();
    chi.correctBoundaryConditions();

    U = Ustar - ex::grad(chi);
    U.correctBoundaryConditions();

    divU = ex::div(U);
    maxdivU = tools::maximum(divU);

    if( parallelCom::master() )
    {   
        std::cout<<"Max velocity divergence after-HH: " << maxdivU<<std::endl;
    }

}
