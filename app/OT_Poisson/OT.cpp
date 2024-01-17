#include "Field/Field.h"
#include <iostream>
#include "Types/vector/vector.h"
#include "ex/grad/grad.h"
#include "ex/div/div.h"
#include "ex/curl/curl.h"
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
    Time time( 0.001, 8, 50 ); //args: dt, endT, write interval / steps

    const scalar pi = tools::pi;
    parallelCom::decompose( settings::zoneName()+"/"+"mesh" ); 

    Mesh mesh( settings::zoneName()+"/"+"mesh", time );

    mesh.write(settings::zoneName()+"/data");
    
    Field<vector> U( mesh, vector(0,0,0), "U" );
    Field<vector> Ustar( mesh, vector(0,0,0), "U" );
    Field<vector> B( mesh, vector(0,0,0), "B" );
    Field<vector> Bstar( mesh, vector(0,0,0), "B" );
    Field<vector> J( mesh, vector(0,0,0), "B" );
    Field<vector> omegav( mesh, vector(0,0,0), "U" );

    std::shared_ptr<Field<scalar> > p_ptr( std::make_shared<Field<scalar> >( mesh, 0, "p" ) );
    auto& p = (*p_ptr);

    std::shared_ptr<Field<scalar> > pB_ptr( std::make_shared<Field<scalar> >( mesh, 0, "pB" ) );
    auto& pB = (*pB_ptr);

    scalar mu = 0.01;
    scalar eta = 0.01;

    poisson pEqn(p_ptr);
    poisson pBEqn(pB_ptr);

    //initial conditions
    for( int i=settings::m()/2; i<mesh.ni()-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<mesh.nj()-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<mesh.nk()-settings::m()/2; k++ )
            {
                scalar x = (i-settings::m()/2)*mesh.dx()+mesh.origin().x();
                scalar y = (j-settings::m()/2)*mesh.dy()+mesh.origin().y();
                scalar z = (k-settings::m()/2)*mesh.dz()+mesh.origin().z();
		
		U(i, j, k) = vector((-2.0)*std::sin(y),2.0*std::sin(x),0.0);
                B(i, j, k) = 0.818*vector( ((-2.0)*std::sin(y*2.0))+std::sin(z),((2.0*std::sin(x))+std::sin(z)),( std::sin(x)+std::sin(y) ));
            }
        }
    }

    U.correctBoundaryConditions();
    p.correctBoundaryConditions();
    B.correctBoundaryConditions();
    pB.correctBoundaryConditions();

    std::ofstream data( settings::zoneName()+"/"+"dat.dat");
    scalar EkOld = 0.0;

    boost::timer::cpu_timer timer;

    while( time.run() )
    {
        double start = timer.elapsed().wall / 1e9;

        time++;

        #include "RKLoop.H" 
	#include "BEqn.H"	

        mesh.write(settings::zoneName()+"/data");
        U.write(settings::zoneName()+"/data", "U");
        p.write(settings::zoneName()+"/data", "p"); 
	B.write(settings::zoneName()+"/data", "B");
	J.write(settings::zoneName()+"/data", "J");
	pB.write(settings::zoneName()+"/data", "pB");

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
        scalar Em=0.0;
        scalar Jmax=0.0;
	omegav = ex::curl(U);
        scalar Jav=0.0;
        scalar omega=0.0;
        int n=0;
    
        for( int i=settings::m()/2; i<mesh.ni()-settings::m()/2-1; i++ )
        {
            for( int j=settings::m()/2; j<mesh.nj()-settings::m()/2-1; j++ )
            {
                for( int k=settings::m()/2; k<mesh.nk()-settings::m()/2-1; k++ )
                {
                    Ek += 0.5 * (U(i, j, k).x() * U(i, j, k).x() + U(i, j, k).y() * U(i, j, k).y() + U(i, j, k).z() * U(i, j, k).z() );
                    Em += 0.5 * (B(i, j, k).x() * B(i, j, k).x() + B(i, j, k).y() * B(i, j, k).y() + B(i, j, k).z() * B(i, j, k).z() );
		    Jmax = std::max( Jmax, sqrt(J(i, j, k).x() * J(i, j, k).x() + J(i, j, k).y() * J(i, j, k).y() + J(i, j, k).z() * J(i, j, k).z())); 
		    omega += sqrt(omegav(i, j, k).x() * omegav(i, j, k).x() + omegav(i, j, k).y() * omegav(i, j, k).y() + omegav(i, j, k).z() * omegav(i, j, k).z());
                    Jav += sqrt(J(i, j, k).x() * J(i, j, k).x() + J(i, j, k).y() * J(i, j, k).y() + J(i, j, k).z() * J(i, j, k).z());
		    n++;
                }
            }
        }
    
        reduce( Ek, plusOp<scalar>() );
	reduce( Jmax, maxOp<scalar>() );
        reduce( Em, plusOp<scalar>() );
        reduce( Jav, plusOp<scalar>() );
        reduce( omega, plusOp<scalar>() );
	reduce( n, plusOp<int>() );

        Ek /= n;
	Em /= n;
        Jav /= n;
        omega /= n;

	scalar epsilon=0.0;
        epsilon = -mu*(omega*omega) - eta*(Jav*Jav);

        if( time.curTime() > time.dt() && parallelCom::master() )
        {
            data<<time.curTime()<<" "<<std::setprecision(15)<<Ek<<" "<<std::setprecision(15)<<Em<<" "<<std::setprecision(15)<<Jav<<" "<<std::setprecision(15)<<Jmax<<" "<<std::setprecision(15)<<omega<<" "<<std::setprecision(15)<<epsilon<<" "<<std::setprecision(15)<<-(Ek-EkOld)/time.dt()<<std::endl;
	}

        EkOld = Ek;
 
    }

    std::cout<< timer.elapsed().wall / 1e9 <<std::endl;
    
    mesh.write(settings::zoneName()+"/data");

    return 0;
}
