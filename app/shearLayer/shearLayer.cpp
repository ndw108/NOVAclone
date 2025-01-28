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
    Time time( 0.001, 420, 1 ); //args: dt, endT, write interval / steps

    const scalar pi = tools::pi;
    parallelCom::decompose( settings::zoneName()+"/"+"mesh" );  

    Mesh mesh( settings::zoneName()+"/"+"mesh", time );

    mesh.write(settings::zoneName()+"/data"); 
/*--------------------
    scalar Mach=0.2;
    scalar Re0=1967.0;
    scalar Ri0=0.08;
    scalar Pr=7;

    scalar rho0=999.57;
    scalar c=0.0417;

    scalar u0=c*Mach;
    scalar h0=0.2;//Re0*mu/rho0/u0;
    scalar t0=0.1;//-Ri0*u0*u0/g[2]/h0;
    scalar mu = rho0*u0*h0/Re0;
    scalar ka = mu/Pr;
    scalar eta = mu;
//    vector g( 0, 0, -Ri0*u0*u0/t0/h0 );
    vector g( 0, 0, -0.000011129 );

    Field<vector> U( mesh, vector(1,0,0), "U"  );  
    Field<scalar> T( mesh, 0, "T" );
    Field<scalar> pHyd( mesh,0, "pHyd" );
    Field<vector> Ustar( mesh, vector(0,0,0), "U" );    
//    Field<vector> B( mesh, vector(0,0,0), "B" );
//    Field<vector> Bstar( mesh, vector(0,0,0), "B" );

----------------------*/

	Field<vector> U( mesh, vector(0,0,0), "U"  ); 
	Field<vector> Ustar( mesh, vector(0,0,0), "U" );
	Field<vector> B( mesh, vector(0,0,0), "B" );
	Field<vector> Bstar( mesh, vector(0,0,0), "B" );


	scalar mu = 1.0/28000.0;
	scalar u0=1.0;
	scalar delta = 1.0/28.0;
	scalar cn = 0.001;
	scalar eta = mu*10.0;	


    std::default_random_engine eng( parallelCom::myProcNo() );
    std::uniform_real_distribution<double> dist(-1.0,1.0);

    std::shared_ptr<Field<scalar> > p_ptr( std::make_shared<Field<scalar> >( mesh, (0.87), "p" ) );
    auto& p = (*p_ptr);

    std::shared_ptr<Field<scalar> > pB_ptr( std::make_shared<Field<scalar> >( mesh, 0, "pB" ) );
    auto& pB = (*pB_ptr);

    poisson pEqn(p_ptr);
    poisson pBEqn(pB_ptr);

    {
        //initial conditions
        for( int i=settings::m()/2; i<mesh.ni()-settings::m()/2; i++ )
        {
          for( int k=settings::m()/2; k<mesh.nk()-settings::m()/2; k++ )
          {

            for( int j=settings::m()/2; j<(mesh.nj()-settings::m()/2); j++ )
            { 


//		    scalar z = (k-settings::m()/2)*mesh.dz()+mesh.origin().z();
                    scalar x = (i-settings::m()/2)*mesh.dx()+mesh.origin().x();
		    scalar y = (j-settings::m()/2)*mesh.dy()+mesh.origin().y();

		    scalar psi = u0*std::exp(-(((y-0.5)*(y-0.5))/(delta*delta)))*(std::cos(8.0*pi*x)+std::cos(20.0*pi*x));
		    scalar psi2 = u0*std::exp(-(((y-1.5)*(y-1.5))/(delta*delta)))*(std::cos(8.0*pi*x)+std::cos(20.0*pi*x));
//                    scalar k0 = 2.0*3.1415926536 / (14.0*h0/2.0);
//                    scalar a=0.005;
//                    scalar b=0.4177;

		if ( y < 1.0 )
		{
		    U(i, j, k) = vector((u0*std::tanh(((4.0 * y) - 2.0)/delta)), 0, 0);
                    U(i, j, k).x() += ( cn * mesh.dy() * psi);
                    U(i, j, k).y() += -cn * mesh.dx() * psi;
		    B(i, j, k) = vector((u0*std::tanh(((4.0 * y) - 2.0)/delta)), 0, 0);
		}
		else
		{  
		    U(i, j, k) = vector(-u0*std::tanh(( (4.0 * y) - 6.0)/delta), 0, 0);
                    U(i, j, k).x() += ( cn * mesh.dy() * psi2);
                    U(i, j, k).y() += -cn * mesh.dx() * psi2;
		    B(i, j, k) = vector(-u0*std::tanh(( (4.0 * y) - 6.0)/delta), 0, 0);
		}
	    }

//                    U(i, j, k).x() += ( cn * mesh.dy() * psi);
//		    U(i, j, k).y() += -cn * mesh.dx() * psi;

//                    U(i, j, k).z() += a*u0/2.0 *dist(eng) / std::cosh(2.0*z/h0) / std::cosh(2.0*z/h0);
//		    B(i, j, k) = vector( u0/2.0*std::tanh( 2.0*z / h0 ), 0, 0 );

//                    T(i, j, k) = t0/2.0*std::tanh( 2.0*z / h0 ); //  + a*t0/2.0 *dist(eng) / std::cosh(2.0*z/h0) / std::cosh(2.0*z/h0);
	    }
        }
        


        U.correctBoundaryConditions();
        p.correctBoundaryConditions();
//        T.correctBoundaryConditions();
	B.correctBoundaryConditions();

    	U.write(settings::zoneName()+"/data", "U");
  //      T.write(settings::zoneName()+"/data", "T");
        p.write(settings::zoneName()+"/data", "p");
	B.write(settings::zoneName()+"/data", "B");
    }

   
    scalar RK4[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 };

    std::ofstream data( settings::zoneName()+"/"+"dat.dat");

//    if( parallelCom::master() )
//    {   
//        scalar epsPlus = 6.5e-4*std::pow(u0,3)/h0;
//        scalar LbMinus = std::pow( std::pow( mu/rho0,3 )/epsPlus, 0.25 )*std::pow(Pr,-0.5);
//
//        std::cout<<"Suggested mesh spacing: " << 1.5*LbMinus << " . Actual grid spacing: " << std::max(mesh.dx(), std::max( mesh.dy(), mesh.dz() ) )<<std::endl;
//    } 

    boost::timer::cpu_timer timer;

    while( time.run() )
    {
        double start = timer.elapsed().wall / 1e9;
//        scalar cfl = settings::restart() ? 0.5 : std::min( (0.5-0.0001)/1000.0 * time.timeStep(), 0.5 ); 
//        time.dt() = cfl*std::min(mesh.dx(), std::min( mesh.dy(), mesh.dz() ) )/(u0/Mach+u0/2.0);

        time++;

        if( parallelCom::master() )
            {   
                std::cout << "Step: " << time.timeStep() << ".              Time: " << time.curTime() << std::endl;
                std::cout << "Elapsed wall time for timestep: " <<  timer.elapsed().wall / 1e9 - start <<std::endl;
  //              std::cout << "drhoz: " << drhodz << "pBbottom: " << (t0/2.0*g[2]*rho0) <<std::endl;
	    }   

            tools::CFL( U, mesh );
        
        
        #include "UEqn.H" 
//        #include "TEqn.H"
//	#include "BEqn.H"

	mesh.write(settings::zoneName()+"/data");
        U.write(settings::zoneName()+"/data", "U");
//	T.write(settings::zoneName()+"/data", "T");
	p.write(settings::zoneName()+"/data", "p");
//	B.write(settings::zoneName()+"/data", "B");

//        if( time.writelog() )
//        {
//	    reuseTmp<vector> upt( mesh );
//            std::shared_ptr<Field<vector> > up( upt() );
//            reuseTmp<scalar> Tpt( mesh );
//            std::shared_ptr<Field<scalar> > Tp( Tpt() );
//            #include "layerAverage.H"
//
//            reuseTmp<scalar> epst( mesh );
//            std::shared_ptr<Field<scalar> > eps( epst() );
//            *eps = 2.0*( ex::grad( up ) && ex::grad( up ) )*mu/rho0;
//            (*eps).write("e");
//
//            reuseTmp<scalar> xit( mesh );
//            std::shared_ptr<Field<scalar> > xi( xit() );
//            *xi = 2.0*( ex::grad( Tp ) & ex::grad( Tp ) )*ka/rho0;
//            (*xi).write("x");
//
//            reuseTmp<vector> vortt( mesh );
//            std::shared_ptr<Field<vector> > vort( vortt() );
//            *vort = ex::curl( U );
//            (*vort).write("v");
//     
//            reuseTmp<vector> gpt( mesh );
//            std::shared_ptr<Field<vector> > gp( gpt() );
//            *gp = ex::grad(p);
//            (*gp).write("z");
//
//            mesh.writeSOS();
//	}
    }
        

    std::cout<< timer.elapsed().wall / 1e9 <<std::endl;
    
    mesh.write(settings::zoneName()+"/data");

    return 0;
}
