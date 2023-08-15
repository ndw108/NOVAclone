{
    for( int k=0; k<mesh.nk(); k++ )
    {
        std::vector<vector> proc_avgU( parallelCom::worldSize(), 0.0 );
        std::vector<scalar> proc_avgT( parallelCom::worldSize(), 0.0 );
        std::vector<int> proc_N( parallelCom::worldSize(), 0 );
        
        for( int j=fd::m/2; j<mesh.nj()-fd::m/2-1; j++ )
        {
            for( int i=fd::m/2; i<mesh.ni()-fd::m/2-1; i++ )
            {
                proc_avgU[parallelCom::myProcNo()] += U(i, j, k);
                proc_avgT[parallelCom::myProcNo()] += T(i,j,k);
                proc_N[parallelCom::myProcNo()] += 1;
            }
        }

        all_reduce(proc_avgU, std::plus<vector>());
        all_reduce(proc_avgT, std::plus<scalar>());
        all_reduce(proc_N, std::plus<int>());
        
        vector avgU(0,0,0);
        scalar avgT=0.0;
        int N=0;

        for( int proci=0; proci<parallelCom::ni(); proci++ )
        {
            for( int procj=0; procj<parallelCom::nj(); procj++ )
            {
                avgU+=proc_avgU[parallelCom::procNo(proci, procj, parallelCom::k())];
                avgT+=proc_avgT[parallelCom::procNo(proci, procj, parallelCom::k())];
                N+=proc_N[parallelCom::procNo(proci, procj, parallelCom::k())];
            }
        }

        avgU/=N;
        avgT/=N;
 
        for( int j=0; j<mesh.nj(); j++ )
        {   
            for( int i=0; i<mesh.ni(); i++ )
            {
                (*up)(i, j, k) = U(i, j, k)-avgU;
                (*Tp)(i, j, k) = T(i, j, k)-avgT;
            }
        }

        //up->correctBoundaryConditions();
        //Tp->correctBoundaryConditions();
    }
}

