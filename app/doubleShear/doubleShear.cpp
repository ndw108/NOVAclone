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

#include "H5Cpp.h"

#include <cmath>
#include <memory>
#include <random>
#include <filesystem>


//this solver runs the MHD double-shear layer case
//oblique tearing modes develop if the mesh is refined enough

int main(int argc, char* argv[])
{
    settings::process( argc, argv ); 
    Time time( 0.00016, 6.01, 6250 ); //args: dt, endT, write interval / steps
    
    boost::timer::cpu_timer timer;

    const scalar pi = tools::pi;
    parallelCom::decompose( settings::zoneName()+"/"+"mesh" ); 

    Mesh mesh( settings::zoneName()+"/"+"mesh", time );

    Field<vector> U( mesh, vector(0,0,0), "U" );
    Field<vector> Ustar( mesh, vector(0,0,0), "U" );
    Field<vector> Ustart( mesh, vector(0,0,0), "U" );
    Field<vector> B( mesh, vector(0,0,0), "B" );
    Field<vector> Bstar( mesh, vector(0,0,0), "B" );
    Field<vector> J( mesh, vector(0,0,0), "B" );
    Field<vector> omegav( mesh, vector(0,0,0), "U" );
    Field<vector> psi(mesh, vector(0,0,0), "B" );
    Field<vector> curlPsi(mesh, vector(0,0,0), "B" );
    Field<vector> divCL( mesh, vector(0,0,0), "B");

    std::shared_ptr<Field<scalar> > Umag_ptr(std::make_shared<Field<scalar> >(mesh, 0, "U") );
    auto& Umag = (*Umag_ptr);

    std::shared_ptr<Field<scalar> > p_ptr( std::make_shared<Field<scalar> >( mesh, 0, "p" ) );
    auto& p = (*p_ptr);

    std::shared_ptr<Field<scalar> > pB_ptr( std::make_shared<Field<scalar> >( mesh, 0, "pB" ) );
    auto& pB = (*pB_ptr);

    scalar mu = 0.00005;
    scalar eta = 0.00005;

    poisson pEqn(p_ptr);
    poisson pBEqn(pB_ptr, pEqn.getFft());

    if( parallelCom::master() )
    {
    std::cout << "Timer after poisson creation: " <<  timer.elapsed().wall / 1e9 <<std::endl;
    std::cout << "   " << std::endl;
    }

    std::default_random_engine eng( parallelCom::myProcNo() );
    std::uniform_real_distribution<double> dist(-1.0,1.0);

    //initial conditions
    if( settings::restart() == false )
    { 
      for( int k=settings::m()/2; k<mesh.nk()-settings::m()/2; k++ )
      {
          for( int j=settings::m()/2; j<mesh.nj()-settings::m()/2; j++ )
          {
              for( int i=settings::m()/2; i<mesh.ni()-settings::m()/2; i++ )
              {
                  scalar x = (i-settings::m()/2)*mesh.dx()+mesh.origin().x();
                  scalar y = (j-settings::m()/2)*mesh.dy()+mesh.origin().y();
                  scalar z = (k-settings::m()/2)*mesh.dz()+mesh.origin().z();
			
                  //shear profile
                  double delta = 0.01;    

                  double y1 = pi / 4.0;          // pi/2
                  double y2 = 3.0 * pi / 4.0;    // 3pi/2
    
                  B(i, j, k) = vector(0.0, 0.0, 0.0);
                  B(i, j, k).z() += (1.0);
                  B(i, j, k).x() += 1.0 * ( 
                    std::tanh((y - 2.0 * pi / 4.0) / delta)
                  - std::tanh((y - 3.0 * 2.0 * pi / 4.0) / delta)
                  - 1.0 
                  );  

		  //primary instability
                  double psi_0 = 0.1;

		  // Streamfunction
                  double psi = psi_0 * std::cos(x) * std::cos(2.0 * y);

                  // u' = dψ/dy
                  double up = -2.0 * psi_0 * std::cos(x) * std::sin(2.0 * y);

                  // v' = -dψ/dx
                  double vp = psi_0 * std::sin(x) * std::cos(2.0 * y);

		  //add the magnetic instability
		  B(i,j,k).x() += up;
		  B(i,j,k).y() += vp;

		  //Random seeded noise
		  U(i, j, k) += (0.005*vector(dist(eng) , dist(eng), dist(eng)))*std::exp(-((y-y1)*(y-y1)));
                  U(i, j, k) += (0.005*vector(dist(eng) , dist(eng), dist(eng)))*std::exp(-((y-y2)*(y-y2)));

	      }
          }
      }
    }
    
    //initial conditions from restart
    if( settings::restart() == true )
    {       
    #include "H5restart.H"
    }

    U.correctBoundaryConditions();
    p.correctBoundaryConditions();
    B.correctBoundaryConditions();
    pB.correctBoundaryConditions();

    scalar EkOld = 0.0;
    scalar EmOld = 0.0;

    std::ofstream data( settings::zoneName()+"/"+"dat.dat");

    if (parallelCom::master()) {
        tools::checkOutputDir();
    }   
    MPI_Barrier(MPI_COMM_WORLD);

    if( parallelCom::master() )
    {
    std::cout << "Timer after completing initialisation: " <<  timer.elapsed().wall / 1e9 <<std::endl;
    std::cout << "Time loop now starting " << std::endl;
    std::cout << "   " << std::endl;
    }

    while( time.run() )
    {
        double start = timer.elapsed().wall / 1e9;

	time++;
         #include "UEqn.H" 
         #include "BEqn.H"


       if( (time.writelog()) || (time.write()) )
        {
            if( parallelCom::master() )
            {
                std::cout << "Step: " << time.timeStep() << ".              Time: " << time.curTime() << std::endl;
                std::cout << "Elapsed wall time for timestep: " <<  timer.elapsed().wall / 1e9 - start <<std::endl;
                std::cout << "Kinetic CFL: " << std::endl;
            }

            tools::CFL( U, mesh );

            if( parallelCom::master() )
            {
            std::cout << "Magnetic CFL: " << std::endl;
            }

            tools::CFL( B, mesh );
        }

	//HDF5 data writing
	//wrties t=0, t=1, then every t=0.1
	bool isMultipleOf0p1 = (std::abs(std::round((time.curTime()-time.dt())* 10.0) - (time.curTime()-time.dt()) * 10.0) < 1e-6);
	if( time.write() 
	    || time.writePlusOne() 
	    || time.writeMinusOne() 
	    || ((time.curTime()-time.dt()) > 1.0 + 0.5*time.dt() && isMultipleOf0p1) )
	{
	    #include "H5write.H"
            if( parallelCom::master() )
            {   
            std::cout << "Elapsed wall time after writing data: " <<  timer.elapsed().wall / 1e9 - start <<std::endl;
            std::cout << "   " << std::endl;
            }
	}

        scalar Ek=0.0;
        scalar Em=0.0;
        scalar Jmax=0.0;
        scalar Jx_max=0.0;
        scalar Jy_max=0.0;
        scalar Jz_max=0.0;
        scalar Umax=0.0;
	scalar Bmax=0.0;
        divCL = ex::grad(pB);
        scalar max_divCL=0.0;
        omegav = ex::curl(U);
        scalar Jav=0.0;
        scalar Jx_av=0.0;
        scalar Jy_av=0.0;
        scalar Jz_av=0.0;
        scalar omega=0.0;
        scalar Omega_x_max=0.0;
        scalar Omega_y_max=0.0;
        scalar Omega_z_max=0.0;
        scalar Omega_x_av=0.0;
        scalar Omega_y_av=0.0;
        scalar Omega_z_av=0.0;
        scalar hydroHelicity=0.0;
        scalar crossHelicity=0.0;

        int n=0;
    
        for( int i=settings::m()/2; i<mesh.ni()-settings::m()/2-1; i++ )
        {
            for( int j=settings::m()/2; j<mesh.nj()-settings::m()/2-1; j++ )
            {
                for( int k=settings::m()/2; k<mesh.nk()-settings::m()/2-1; k++ )
                {
                    Ek += 0.5 * (U(i, j, k).x() * U(i, j, k).x() + U(i, j, k).y() * U(i, j, k).y() + U(i, j, k).z() * U(i, j, k).z() );
                    Em += 0.5 * (B(i, j, k).x() * B(i, j, k).x() + B(i, j, k).y() * B(i, j, k).y() + B(i, j, k).z() * B(i, j, k).z() );
		    hydroHelicity += (U(i, j, k).x() * omegav(i, j, k).x()) + (U(i, j, k).y() * omegav(i, j, k).y()) + (U(i, j, k).z() * omegav(i, j, k).z());
                    crossHelicity += (U(i, j, k).x() * B(i, j, k).x()) + (U(i, j, k).y() * B(i, j, k).y()) + (U(i, j, k).z() * B(i, j, k).z());
                    max_divCL = std::max( max_divCL, sqrt((divCL(i,j,k).x() * divCL(i,j,k).x())));
                    max_divCL = std::max( max_divCL, sqrt((divCL(i,j,k).y() * divCL(i,j,k).y())));
                    max_divCL = std::max( max_divCL, sqrt((divCL(i,j,k).z() * divCL(i,j,k).z())));
                    Jx_av += sqrt(J(i, j, k).x() * J(i, j, k).x());
                    Jy_av += sqrt(J(i, j, k).y() * J(i, j, k).y());
                    Jz_av += sqrt(J(i, j, k).z() * J(i, j, k).z());
                    Jx_max = std::max( Jx_max, sqrt(J(i, j, k).x() * J(i, j, k).x()));
                    Jy_max = std::max( Jy_max, sqrt(J(i, j, k).y() * J(i, j, k).y()));
                    Jz_max = std::max( Jz_max, sqrt(J(i, j, k).z() * J(i, j, k).z()));
                    Omega_x_av += sqrt(omegav(i, j, k).x() * omegav(i, j, k).x());
                    Omega_y_av += sqrt(omegav(i, j, k).y() * omegav(i, j, k).y());
                    Omega_z_av += sqrt(omegav(i, j, k).z() * omegav(i, j, k).z());
                    Omega_x_max = std::max( Omega_x_max, sqrt(omegav(i, j, k).x() * omegav(i, j, k).x()));
                    Omega_y_max = std::max( Omega_y_max, sqrt(omegav(i, j, k).y() * omegav(i, j, k).y()));
                    Omega_z_max = std::max( Omega_z_max, sqrt(omegav(i, j, k).z() * omegav(i, j, k).z()));
                    Jmax = std::max( Jmax, sqrt(J(i, j, k).x() * J(i, j, k).x() + J(i, j, k).y() * J(i, j, k).y() + J(i, j, k).z() * J(i, j, k).z()));
                    Umax = std::max( Umax, sqrt(U(i, j, k).x() * U(i, j, k).x() + U(i, j, k).y() * U(i, j, k).y() + U(i, j, k).z() * U(i, j, k).z()));
                    Bmax = std::max( Bmax, sqrt(B(i, j, k).x() * B(i, j, k).x() + B(i, j, k).y() * B(i, j, k).y() + B(i, j, k).z() * B(i, j, k).z()));
		    omega += sqrt(omegav(i, j, k).x() * omegav(i, j, k).x() + omegav(i, j, k).y() * omegav(i, j, k).y() + omegav(i, j, k).z() * omegav(i, j, k).z());
                    Jav += sqrt(J(i, j, k).x() * J(i, j, k).x() + J(i, j, k).y() * J(i, j, k).y() + J(i, j, k).z() * J(i, j, k).z());
		    n++;
                }
            }
        }

	reduce( Ek, plusOp<scalar>() );
        reduce( Jmax, maxOp<scalar>() );
        reduce( Jx_max, maxOp<scalar>() );
        reduce( Jy_max, maxOp<scalar>() );
        reduce( Jz_max, maxOp<scalar>() );
        reduce( Umax, maxOp<scalar>() );
	reduce( Bmax, maxOp<scalar>() );
        reduce( Em, plusOp<scalar>() );
        reduce( Jav, plusOp<scalar>() );
        reduce( Jx_av, plusOp<scalar>() );
        reduce( Jy_av, plusOp<scalar>() );
        reduce( Jz_av, plusOp<scalar>() );
        reduce( omega, plusOp<scalar>() );
        reduce( Omega_x_av, plusOp<scalar>() );
        reduce( Omega_y_av, plusOp<scalar>() );
        reduce( Omega_z_av, plusOp<scalar>() );
        reduce( Omega_x_max, plusOp<scalar>() );
        reduce( Omega_y_max, plusOp<scalar>() );
        reduce( Omega_z_max, plusOp<scalar>() );
        reduce( max_divCL, plusOp<scalar>() );
        reduce( hydroHelicity, plusOp<scalar>() );
        reduce( crossHelicity, plusOp<scalar>() );
        reduce( n, plusOp<int>() );

        Ek /= n;
        Em /= n;
        Jav /= n;
        Jx_av /= n;
        Jy_av /= n;
        Jz_av /= n;
        omega /= n;
        Omega_x_av /= n;
        Omega_y_av /= n;
        Omega_z_av /= n;
        hydroHelicity /= n;
        crossHelicity /= n;

	scalar epsilon=0.0;
        epsilon = -mu*(omega*omega) - eta*(Jav*Jav);

        if( time.curTime() > time.dt() && parallelCom::master() )
        {
            data<<time.curTime()<<" "<<std::setprecision(15)<<Ek<<" "<<std::setprecision(15)<<Em<<" "<<std::setprecision(15)<<Jav<<" "<<std::setprecision(15)
                    <<Jmax<<" "<<std::setprecision(15)<<Jx_max<<" "<<std::setprecision(15)<<Jy_max<<" "<<std::setprecision(15)<<Jz_max<<" "<<std::setprecision(15)<<
                    Jx_av<<" "<<std::setprecision(15)<<Jy_av<<" "<<std::setprecision(15)<<Jz_av<<" "<<std::setprecision(15)<<omega<<" "<<std::setprecision(15)
                    <<Omega_x_av<<" "<<std::setprecision(15)<<Omega_y_av<<" "<<std::setprecision(15)<<Omega_z_av<<" "<<std::setprecision(15)
                    <<Omega_x_max<<" "<<std::setprecision(15)<<Omega_y_max<<" "<<std::setprecision(15)<<Omega_z_max<<" "<<std::setprecision(15)
                    <<hydroHelicity<<" "<<std::setprecision(15)<<crossHelicity<<" "<<std::setprecision(15)
                    <<epsilon<<" "<<std::setprecision(15)<<-(Ek-EkOld)/time.dt()<<std::setprecision(15)
                    <<" "<<(Em-EmOld)/time.dt()<<std::setprecision(15)<<" "<<Umax<<" "<<Bmax<<" "<<max_divCL<<" "<<std::endl;
        }

        EkOld = Ek;
	EmOld = Em;

    }

    std::cout<< timer.elapsed().wall / 1e9 <<std::endl;
    
    //mesh.write(settings::zoneName()+"/data");

    return 0;
}
