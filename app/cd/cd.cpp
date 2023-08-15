/*  This file is part of NOVA.

    NOVA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    NOVA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with NOVA.  If not, see <http://www.gnu.org/licenses/>.

    Author: Alex Skillen. alex.skillen@manchester.ac.uk

*/
#include "Field/Field.h"
#include <iostream>
#include "Types/vector/vector.h"
#include "fdc/grad/grad.h"
#include "fdc/curl/curl.h"
#include "fdc/div/div.h"
#include "fdc/laplacian/laplacian.h"
#include "Time/Time.h"
#include "parallelCom/parallelCom.h"
#include "settings/settings.h"
#include "BC/BC.h"

#include <boost/timer/timer.hpp>

#include <cmath>
#include <random>

#include "gperftools/profiler.h"

int main(int argc, char* argv[])
{
    settings::process( argc, argv );

    double L=1;
    double nu=0.01;
    double pe = L/nu;

    Time time( 1e-5, 0.7, 0.1 );

    parallelCom::decompose( settings::zoneName()+"/"+"mesh" );

    Mesh mesh( settings::zoneName()+"/"+"mesh", time );

    mesh.write(settings::zoneName()+"/data");

    Field<vector> U( mesh, vector(0,0,1), settings::zoneName()+"/"+"UBC" );
    Field<scalar> p( mesh,0, settings::zoneName()+"/"+"pBC" );

    p.correctBoundaryConditions();
    p.write(settings::zoneName()+"/data", "p"); 
   
    scalar RK4[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 };

    boost::timer::cpu_timer timer;

    while( time.run() )
    {

        time++;

        if( parallelCom::master() && time.timeStep() % 1000 == 0 )
        {
            std::cout<< "Time: " << time.curTime() <<std::endl;
        }
        
        {
            reuseTmp<scalar> dpt( mesh );
            std::shared_ptr<Field<scalar> > dp( dpt() );
  
            reuseTmp<scalar> pct( mesh );
            std::shared_ptr<Field<scalar> > pc( pct() );

            reuseTmp<scalar> pOldt( mesh );
            std::shared_ptr<Field<scalar> > pOld( pOldt() );

            *pc = p;
            *pOld = p;

            for( int rk=0; rk<4; rk++ )
            {
                *dp =
                time.dt() *
                (
                 -fdc::div( U, p ) 
                 +fdc::laplacian( nu, p )
                );

                *pc = pc+RK4[rk]*dp;

                if( rk<3 )
                {
                    p = pOld + (int(rk/2)+1)/2.0 * dp;
                }
                else
                {
                    p=pc;
                }
            
                p.correctBoundaryConditions();
            }
        }
        
        mesh.write(settings::zoneName()+"/data");
        p.write(settings::zoneName()+"/data", "p");

    }
        
    std::cout<<"Pe: "<<pe<<std::endl;

    double rmsE=0.0;

    for( int k=settings::m()/2; k<mesh.nk()-settings::m()/2; k++ )
    {
        double x=mesh.z(k);
        double p_exac = ( std::exp( x*pe/L ) -1 ) / (std::exp(pe)-1 );
        double err = (p(settings::m()/2, settings::m()/2, k)-p_exac)/(p_exac+1e-10);
        rmsE += err*err;
    }

    reduce( rmsE, plusOp<scalar>() );
    int N = (mesh.nk()-settings::m());
    reduce( N, plusOp<int>() );

    if( parallelCom::master() )
    {
        rmsE/=N;
        rmsE = std::pow( rmsE, 0.5 );

        std::cout<<"RMS Error: "<< rmsE <<std::endl;

        std::cout<< timer.elapsed().wall / 1e9 <<std::endl;
    }    

    mesh.write(settings::zoneName()+"/data");

    return 0;
}
