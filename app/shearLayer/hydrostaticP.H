{
    reuseTmp<scalar> pHydt( mesh );
    std::shared_ptr<Field<scalar> > pHyd( pHydt() );

    for( int i=fd::m/2; i<mesh.ni()-fd::m/2-1; i++ )
    {
        for( int j=fd::m/2; j<mesh.nj()-fd::m/2; j++ )
        {
            long double pHydc=0.0L;

            for( int k=fd::m/2; k<mesh.nk()-fd::m/2-1; k++ )
            {
                pHydc += T(i, j, k)*mesh.dz()*g[2]*rho0;
            }

            std::vector<scalar> cumTot( parallelCom::worldSize(), scalar(0.0) );
            cumTot[parallelCom::myProcNo()] = pHydc;

            all_reduce( cumTot, std::plus<scalar>() );

            pHydc = 0.0;
            long double aboveProcs=0.0L;

            for( int pk=parallelCom::k()+1; pk<parallelCom::nk(); pk++ )
            {
                aboveProcs += cumTot[parallelCom::procNo( parallelCom::i(), parallelCom::j(), pk )];
            } 

            for( int k=mesh.nk()-2-fd::m/2; k>=fd::m/2; k-- )
            {
                pHydc += T(i, j, k)*mesh.dz()*g[2]*rho0;
                (*pHyd)(i, j, k) += pHydc + aboveProcs;
            }

        }
    }

    pHyd->correctBoundaryConditions();

    p += pHyd;
}

