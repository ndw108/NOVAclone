{
    for( int i=settings::m()/2; i<mesh.ni()-settings::m()/2-1; i++ )
    {
        for( int j=settings::m()/2; j<mesh.nj()-settings::m()/2; j++ )
        {
            double pHydc=0.0;

            for( int k=settings::m()/2; k<mesh.nk()-settings::m()/2-1; k++ )
            {
                pHydc += T(i, j, k)*mesh.dz()*g[2]*rho0;
            }

            std::vector<scalar> cumTot( parallelCom::worldSize(), scalar(0.0) );
            cumTot[parallelCom::myProcNo()] = pHydc;

            all_reduce( cumTot, std::plus<scalar>() );

            pHydc = 0.0;
            double aboveProcs=0.0;

            for( int pk=parallelCom::k()+1; pk<parallelCom::nk(); pk++ )
            {
                aboveProcs += cumTot[parallelCom::procNo( parallelCom::i(), parallelCom::j(), pk )];
            } 

            for( int k=mesh.nk()-2-settings::m()/2; k>=settings::m()/2; k-- )
            {
                pHydc += T(i, j, k)*mesh.dz()*g[2]*rho0;
                pHyd(i, j, k) = pHydc; // + aboveProcs;
            }

        }
    }
//add field and implement null boundary condition
//    pHyd.correctBoundaryConditions();

    p += pHyd;
}

