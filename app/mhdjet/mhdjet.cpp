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

#include "H5cpp.h"

#include <cmath>
#include <memory>
#include <random>

    //H5 Create
    //H5 Write dataset using MPI
    void writeVectorFieldToHDF5( MPI_Comm comm, const std::string &fileName, const std::vector<std::vector<std::vector<double>>> &vectorField) {
    
    MPI_Comm world = parallelCom::world();
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Dimensions of the vector field
    hsize_t hni = vectorField.size();
    hsize_t hnj = vectorField[0].size();
    hsize_t hnk = vectorField[0][0].size();

    //prepare data
    std::vector<double> data(hni * hnj * hnk );  
    int n=0;
    for (int i = 0; i < hni; ++i) {
        for (int j = 0; j < hnj; ++j) {
            for (int k = 0; k < hnk; ++k) {
            data[n] = vectorField[i][j][k];
            n++;
            }
        }
    }

    // Create an HDF5 file
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, parallelCom::world(), MPI_INFO_NULL);
    hid_t file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

    // Create a data space
    hsize_t dims[] = {   static_cast<hsize_t>((hni)*parallelCom::ni() ),
			 static_cast<hsize_t>((hnj)*parallelCom::nj() ),
			 static_cast<hsize_t>((hnk)*parallelCom::nk() )  };  // 3 for the vector components (x, y, z)
    hid_t dataspace_id = H5Screate_simple(3, dims, NULL);

    //Create a dataset
    hid_t dataset_id = H5Dcreate2(file_id, "/3DArray", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Create the local dataspace for each MPI rank
    hsize_t local_dims[3] = {hni, hnj, hnk};
    hsize_t start[3] = {static_cast<hsize_t>((hni*parallelCom::i()) ) , static_cast<hsize_t>((hnj*parallelCom::j()) ) , static_cast<hsize_t>((hnk*parallelCom::k()) ) };
    hsize_t count[3] = {hni, hnj, hnk};
    hid_t memspace_id = H5Screate_simple(3, local_dims, NULL);


    //Sanity Checks for data indexing
//    std::cout << "  ni:  " << parallelCom::ni() << "  nj:  " << parallelCom::nj() << "  nk:  " << parallelCom::nk() << "  i:  " <<parallelCom::i()  << "  j:  " <<parallelCom::j()<< "  k:  " <<parallelCom::k() << "  hni  " << hni << "  hnj  " << hnj << "  hnk  " << hnk << "  dims_x  " << dims[0] << "  dims_y  " << dims[1] << "  dims_z  " << dims[2] << "  LocalStart_x  " << start[0] << "  LocalStart_y  " << start[1]<< "  LocalStart_z  " << start[2] << std::endl;
//    std::cout << " myProNo_x " << parallelCom::myProcNo() << std::endl;


    // Select the hyperslab to write
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL);

    // Write the data from each processor into the HDF5 dataset
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, data.data());

    // Close the dataset and file
    // Close resources
    H5Pclose(plist_id);
    H5Sclose(dataspace_id);
    H5Sclose(memspace_id);
    H5Dclose(dataset_id);
    H5Fclose(file_id);

}


/*
//H5 Write dataset in serial
    void writeVectorFieldToHDF5(const std::string &fileName, const Field<vector> &vectorField) {
    // Dimensions of the vector field
    hsize_t hni = vectorField[0]-settings::m();
    hsize_t hnj = vectorField[1]-settings::m();
    hsize_t hnk = vectorField[2]-settings::m();

    // Create an HDF5 file
    H5::H5File file(fileName, H5F_ACC_TRUNC);

    // Create a data space
    hsize_t dims[3] = {hni, hnj, hnk };  // 3 for the vector components (x, y, z)
    H5::DataSpace dataspace(3, dims);

   // Create a dataset
    H5::DataSet dataset = file.createDataSet("VectorField", H5::PredType::IEEE_F64LE, dataspace);

    // Prepare data to write
    std::vector<double> data(hni * hnj * hnk * 3); 
    for (int i = 0; i < hni; ++i) {
        for (int j = 0; j < hnj; ++j) {
            for (int k = 0; k < hnk; ++k) {
                data[(i * hnj * hnk + j * hnk + k)  + 0] = vectorField[0];  // x component
                data[(i * hnj * hnk + j * hnk + k)  + 1] = vectorField[1];  // y component
                data[(i * hnj * hnk + j * hnk + k)  + 2] = vectorField[2];  // z component
                
            }
        }
    }

    // Write data to the dataset
    dataset.write(data.data(), H5::PredType::IEEE_F64LE);

    // Close the dataset and file
    dataset.close();
    file.close();
}
*/



int main(int argc, char* argv[])
{
    settings::restart() = true;

    settings::process( argc, argv ); 
    Time time( 0.02, 160.02, 50 ); //args: dt, endT, write interval / steps

    const scalar pi = tools::pi;
    parallelCom::decompose( settings::zoneName()+"/"+"mesh" ); 

    Mesh mesh( settings::zoneName()+"/"+"mesh", time );
    
    mesh.write(settings::zoneName()+"/data");


    Field<vector> U( mesh, vector(0,0,0), "U" );
    Field<vector> Ustar( mesh, vector(0,0,0), "U" );
    Field<vector> Ustart( mesh, vector(0,0,0), "U" );
    Field<vector> B( mesh, vector(0,0,0), "B" );
    Field<vector> Bstar( mesh, vector(0,0,0), "B" );
    Field<vector> J( mesh, vector(0,0,0), "B" );
    Field<vector> omegav( mesh, vector(0,0,0), "U" );

    std::shared_ptr<Field<scalar> > Umag_ptr(std::make_shared<Field<scalar> >(mesh, 0, "U") );
    auto& Umag = (*Umag_ptr);

    std::shared_ptr<Field<scalar> > p_ptr( std::make_shared<Field<scalar> >( mesh, 0, "p" ) );
    auto& p = (*p_ptr);

    std::shared_ptr<Field<scalar> > pB_ptr( std::make_shared<Field<scalar> >( mesh, 0, "pB" ) );
    auto& pB = (*pB_ptr);

    std::shared_ptr<Field<scalar> > fftx_ptr( std::make_shared<Field<scalar> >( mesh, 0, "U" ) );
    auto& fftx = (*fftx_ptr);

    scalar mu = 0.000285714;
    scalar eta = 0.000285714;

    poisson fftxDo(fftx_ptr);
    poisson pEqn(p_ptr);
    poisson pBEqn(pB_ptr);

    std::default_random_engine eng( parallelCom::myProcNo() );
    std::uniform_real_distribution<double> dist(-1.0,1.0);

    //initial conditions
    for( int i=settings::m()/2; i<mesh.ni()-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<mesh.nj()-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<mesh.nk()-settings::m()/2; k++ )
            {
                scalar x = (i-settings::m()/2)*mesh.dx()+mesh.origin().x();
                scalar y = -10+(j-settings::m()/2)*mesh.dy()+mesh.origin().y();
                scalar z = (k-settings::m()/2)*mesh.dz()+mesh.origin().z();
		
		int perU=0.0;
                for( int px=0; px<2; px++ )
                {   
                    perU += dist(eng)*sin((px*x)-(dist(eng)*13.96));
                }

		for( int px=3; px<mesh.ni()/3; px++ )
		{
		    perU += dist(eng)*sin((px*x)-(dist(eng)*13.96));
		}

                int perZ=0.0;
                for( int pz=0; pz<2; pz++ )
                {   
                    perZ += dist(eng)*sin((pz*z)-(dist(eng)*13.96));
                }   

                for( int pz=3; pz<mesh.nk()/3; pz++ )
                {   
                    perZ += dist(eng)*sin((pz*z)-(dist(eng)*13.96));
                } 

		U(i, j, k) = vector(1.0/(std::cosh(y)*std::cosh(y)), 0.0, 0.0);
		U(i, j, k).y() += (0.0001*dist(eng))*std::exp(-(y*y));
		U(i, j, k).x() += (0.001*std::sin(x*0.9) + 0.0001*dist(eng))*std::exp(-(y*y));
		U(i, j, k).z() += (0.0001*dist(eng))*std::exp(-(y*y));		

		B(i, j, k) = vector(0.015, 0.0, 0.0);

	    }
        }
    }

    U.correctBoundaryConditions();
    p.correctBoundaryConditions();
    B.correctBoundaryConditions();
    pB.correctBoundaryConditions();

    std::ofstream data( settings::zoneName()+"/"+"dat.dat");
    scalar EkOld = 0.0;
    scalar EmOld = 0.0;

    boost::timer::cpu_timer timer;


    while( time.run() )
    {
        double start = timer.elapsed().wall / 1e9;


        mesh.write(settings::zoneName()+"/data");
        U.write(settings::zoneName()+"/data", "U");
        p.write(settings::zoneName()+"/data", "p"); 
	B.write(settings::zoneName()+"/data", "B");
	J.write(settings::zoneName()+"/data", "J");
	pB.write(settings::zoneName()+"/data", "pB");

        time++;
        #include "UEqn.H" 
        #include "BEqn.H"

        if( time.writelog() )
        {
            if( parallelCom::master() )
            { 
                std::cout << "Step: " << time.timeStep() << ".              Time: " << time.curTime() << std::endl;
                std::cout << "Elapsed wall time for timestep: " <<  timer.elapsed().wall / 1e9 - start <<std::endl;
            }

            tools::CFL( U, mesh );
        }


//	HDF5 data writing
	if( time.write() )
        {   
        int hni = mesh.glob_n()[0]/parallelCom::ni();
        int hnj = mesh.glob_n()[1]/parallelCom::nj();
        int hnk = mesh.glob_n()[2]/parallelCom::nk();

	std::vector<std::vector<std::vector<double>>> vectorField(hni, std::vector<std::vector<double>>(hnj, std::vector<double>(hnk)));
	//write Ux
	int n=0;
	for (int k = 0; k < hnk; ++k) {
        	for (int j = 0; j < hnj; ++j) {
            	    for (int i = 0; i < hni; ++i) {
                	vectorField[i][j][k] = U( (i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2) ).x();  // x component
            		n++;
			}
        	}
    	}
	std::string H5title = "./h5data/Ux" +  std::to_string( (time.curTime()-time.dt()) ) + ".h5";
	writeVectorFieldToHDF5(parallelCom::world(), H5title, vectorField);

	//write Uy
        n=0;
        for (int k = 0; k < hnk; ++k) {
                for (int j = 0; j < hnj; ++j) {
                    for (int i = 0; i < hni; ++i) {
                        vectorField[i][j][k] = U( (i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2) ).y();  // x component
                        n++;
                        }   
                }   
        }	
        H5title = "./h5data/Uy_" +  std::to_string( (time.curTime()-time.dt()) ) + ".h5";
        writeVectorFieldToHDF5(parallelCom::world(), H5title, vectorField);

        //write Uz
        n=0;
        for (int k = 0; k < hnk; ++k) {
                for (int j = 0; j < hnj; ++j) {
                    for (int i = 0; i < hni; ++i) {
                        vectorField[i][j][k] = U( (i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2) ).z();  // x component
                        n++;
                        }   
                }   
        }           
        H5title = "./h5data/Uz_" +  std::to_string( (time.curTime()-time.dt()) ) + ".h5";
        writeVectorFieldToHDF5(parallelCom::world(), H5title, vectorField);

        //write p
        n=0;
        for (int k = 0; k < hnk; ++k) {
                for (int j = 0; j < hnj; ++j) {
                    for (int i = 0; i < hni; ++i) {
                        vectorField[i][j][k] = p( (i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2) );  // x component
                        n++;
                        }   
                }   
        }           
        H5title = "./h5data/p_" +  std::to_string( (time.curTime()-time.dt()) ) + ".h5";
        writeVectorFieldToHDF5(parallelCom::world(), H5title, vectorField);

        //write Bx
        n=0;
        for (int k = 0; k < hnk; ++k) {
                for (int j = 0; j < hnj; ++j) {
                    for (int i = 0; i < hni; ++i) {
                        vectorField[i][j][k] = B( (i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2) ).x();  // x component
                        n++;
                        }   
                }   
        }    
        H5title = "./h5data/Bx_" +  std::to_string( (time.curTime()-time.dt()) ) + ".h5";
        writeVectorFieldToHDF5(parallelCom::world(), H5title, vectorField);

        //write By
        n=0;
        for (int k = 0; k < hnk; ++k) {
                for (int j = 0; j < hnj; ++j) {
                    for (int i = 0; i < hni; ++i) {
                        vectorField[i][j][k] = B( (i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2) ).y();  // x component
                        n++;
                        }
                }
        }
        H5title = "./h5data/By_" +  std::to_string( (time.curTime()-time.dt()) ) + ".h5";
        writeVectorFieldToHDF5(parallelCom::world(), H5title, vectorField);

        //write Bz
        n=0;
        for (int k = 0; k < hnk; ++k) {
                for (int j = 0; j < hnj; ++j) {
                    for (int i = 0; i < hni; ++i) {
                        vectorField[i][j][k] = B( (i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2) ).z();  // x component
                        n++;
                        }
                }
        }
        H5title = "./h5data/Bz_" +  std::to_string( (time.curTime()-time.dt()) ) + ".h5";
        writeVectorFieldToHDF5(parallelCom::world(), H5title, vectorField);

	//write Jx
        n=0;
        for (int k = 0; k < hnk; ++k) {
                for (int j = 0; j < hnj; ++j) {
                    for (int i = 0; i < hni; ++i) {
                        vectorField[i][j][k] = J( (i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2) ).x();  // x component
                        n++;
                        }   
                }   
        }   
        H5title = "./h5data/Jx_" +  std::to_string( (time.curTime()-time.dt()) ) + ".h5";
        writeVectorFieldToHDF5(parallelCom::world(), H5title, vectorField);

        //write Jy
        n=0;
        for (int k = 0; k < hnk; ++k) {
                for (int j = 0; j < hnj; ++j) {
                    for (int i = 0; i < hni; ++i) {
                        vectorField[i][j][k] = J( (i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2) ).y();  // x component
                        n++;
                        }
                }
        }
        H5title = "./h5data/Jy_" +  std::to_string( (time.curTime()-time.dt()) ) + ".h5";
        writeVectorFieldToHDF5(parallelCom::world(), H5title, vectorField);

        //write Jz
        n=0;
        for (int k = 0; k < hnk; ++k) {
                for (int j = 0; j < hnj; ++j) {
                    for (int i = 0; i < hni; ++i) {
                        vectorField[i][j][k] = J( (i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2) ).z();  // x component
                        n++;
                        }
                }
        }
        H5title = "./h5data/Jz_" +  std::to_string( (time.curTime()-time.dt()) ) + ".h5";
        writeVectorFieldToHDF5(parallelCom::world(), H5title, vectorField);





	if( parallelCom::master() )
	  {
          std::cout << "Vector field written to vector_field.h5" << std::endl;
	  }
	}


        scalar Ek=0.0;
        scalar Em=0.0;
        scalar Jmax=0.0;
	scalar Umax=0.0;
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
		    Umax = std::max( Umax, sqrt(U(i, j, k).x() * U(i, j, k).x() + U(i, j, k).y() * U(i, j, k).y() + U(i, j, k).z() * U(i, j, k).z()));
		    omega += sqrt(omegav(i, j, k).x() * omegav(i, j, k).x() + omegav(i, j, k).y() * omegav(i, j, k).y() + omegav(i, j, k).z() * omegav(i, j, k).z());
                    Jav += sqrt(J(i, j, k).x() * J(i, j, k).x() + J(i, j, k).y() * J(i, j, k).y() + J(i, j, k).z() * J(i, j, k).z());
		    n++;
                }
            }
        }
    
        reduce( Ek, plusOp<scalar>() );
	reduce( Jmax, maxOp<scalar>() );
	reduce( Umax, maxOp<scalar>() );
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
            data<<time.curTime()<<" "<<std::setprecision(15)<<Ek<<" "<<std::setprecision(15)<<Em<<" "<<std::setprecision(15)<<Jav<<" "<<std::setprecision(15)<<Jmax<<" "<<std::setprecision(15)<<omega<<" "<<std::setprecision(15)<<epsilon<<" "<<std::setprecision(15)<<-(Ek-EkOld)/time.dt()<<std::setprecision(15)<<" "<<(Em-EmOld)/time.dt()<<std::setprecision(15)<<" "<<Umax<<std::endl;
	}

        EkOld = Ek;
	EmOld = Em;
 
    }


    std::cout<< timer.elapsed().wall / 1e9 <<std::endl;
    
    mesh.write(settings::zoneName()+"/data");

    return 0;
}
