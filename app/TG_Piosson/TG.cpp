#include "Field/Field.h"
#include <iostream>
#include "Types/vector/vector.h"
#include "ex/grad/grad.h"
#include "ex/div/div.h"
#include "ex/laplacian/laplacian.h"
#include "Time/Time.h"
#include "parallelCom/parallelCom.h"
#include "settings/settings.h"
#include "BC/BC.h"
#include <functional>
#include "Tools/Tools.h"
#include "poisson/poisson.h"

#include <boost/timer/timer.hpp>

#include <cmath>
#include <memory>
#include <random>


int main(int argc, char* argv[])
{
    settings::process( argc, argv ); 
    Time time( 0.01, 20, 1.0 ); //args: dt, endT, write interval / steps

    const scalar pi = tools::pi;
    parallelCom::decompose( settings::zoneName()+"/"+"mesh" ); 

    Mesh mesh( settings::zoneName()+"/"+"mesh", time );

    mesh.write(settings::zoneName()+"/data");
    
    Field<vector> U( mesh, vector(0,0,0), "U" );
    Field<vector> Ustar( mesh, vector(0,0,0), "U" );

    std::shared_ptr<Field<scalar> > p_ptr( std::make_shared<Field<scalar> >( mesh, 0, "p" ) );
    auto& p = (*p_ptr); 

    scalar mu = 1.0/1600.0;

    poisson pEqn(p_ptr);

    //initial conditions
    for( int i=settings::m()/2; i<mesh.ni()-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<mesh.nj()-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<mesh.nk()-settings::m()/2; k++ )
            {
                scalar x = -pi+(i-settings::m()/2)*mesh.dx()+mesh.origin().x();
                scalar y = -pi+(j-settings::m()/2)*mesh.dy()+mesh.origin().y();
                scalar z = -pi+(k-settings::m()/2)*mesh.dz()+mesh.origin().z();

                U(i, j, k) = vector(std::sin(x)*std::cos(y)*std::cos(z),-std::cos(x)*std::sin(y)*std::cos(z),0);
	    }
        }
    }

    U.correctBoundaryConditions();
    p.correctBoundaryConditions();


    std::ofstream data( settings::zoneName()+"/"+"dat.dat");
    scalar EkOld = 0.0;

    boost::timer::cpu_timer timer;

    while( time.run() )
    {
        double start = timer.elapsed().wall / 1e9;

        time++;

        #include "RKLoop.H" 
       
        mesh.write(settings::zoneName()+"/data");
        U.write(settings::zoneName()+"/data", "U");
        p.write(settings::zoneName()+"/data", "p"); 


        if( time.writelog() )
        {
            if( parallelCom::master() )
            { 
                std::cout << "Step: " << time.timeStep() << ".              Time: " << time.curTime() << std::endl;
                std::cout << "Elapsed wall time for timestep: " <<  timer.elapsed().wall / 1e9 - start <<std::endl;
            }

            tools::CFL( U, mesh );
        }

        scalar Ek=0.0;
        int n=0;
    
        for( int i=settings::m()/2; i<mesh.ni()-settings::m()/2-1; i++ )
        {
            for( int j=settings::m()/2; j<mesh.nj()-settings::m()/2-1; j++ )
            {
                for( int k=settings::m()/2; k<mesh.nk()-settings::m()/2-1; k++ )
                {
                    Ek += 0.5 * (U(i, j, k).x() * U(i, j, k).x() + U(i, j, k).y() * U(i, j, k).y() + U(i, j, k).z() * U(i, j, k).z() );
                    n++;
                }
            }
        }
    
        reduce( Ek, plusOp<scalar>() );
        reduce( n, plusOp<int>() );

        Ek /= n;

        if( time.curTime() > time.dt() && parallelCom::master() )
        {
            data<<time.curTime()<<" "<<std::setprecision(15)<<Ek<<" "<<std::setprecision(15)<<-(Ek-EkOld)/time.dt()<<std::endl;
        }

        EkOld = Ek;
 
    }

    std::cout<< timer.elapsed().wall / 1e9 <<std::endl;
    
    mesh.write(settings::zoneName()+"/data");

    return 0;
}
