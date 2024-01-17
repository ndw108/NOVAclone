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
#include "parallelCom/parallelCom.h"
#include <sstream>
#include <fstream>

#ifdef HAVE_MPI
boost::mpi::environment parallelCom::env_;
#ifdef HAVE_MUI
boost::mpi::communicator parallelCom::world_(mui::mpi_split_by_app(), boost::mpi::comm_take_ownership);
#else
boost::mpi::communicator parallelCom::world_;
#endif
int parallelCom::myRank_= world_.rank();
int parallelCom::worldSize_= world_.size();
int parallelCom::ni_= 1;
int parallelCom::nj_= 1;
int parallelCom::nk_= world_.size();
int parallelCom::i_= 0;
int parallelCom::j_= 0;
int parallelCom::k_= myRank_;

#else
int parallelCom::myRank_= 0;
int parallelCom::worldSize_= 1;
int parallelCom::ni_= 1;
int parallelCom::nj_= 1;
int parallelCom::nk_= 1;
int parallelCom::i_= 0;
int parallelCom::j_= 0;
int parallelCom::k_= 0;
#endif

#ifdef HAVE_MUI
std::vector<mui::uniface<mui::default_config>*> parallelCom::muiInterfaces_;
std::array<std::vector<std::unique_ptr<shape> >, 10> parallelCom::muiShapesSend_;
std::array<std::vector<std::unique_ptr<shape> >, 10> parallelCom::muiShapesRecv_;
std::array<vector, 10> parallelCom::muiTrans_;
#endif



void parallelCom::decompose(int ni, int nj, int nk)
{
    int np=worldSize_;
    
    double lsq=1e10;

    for( int ii=1; ii<=np; ii++ )
    {
        for( int jj=1; jj<=np/ii; jj++ )
        {
            if( std::fabs( double(np)/double(ii)/double(jj) - np/ii/jj ) < 1e-5  ) 
            {
                int kk=np/ii/jj;

                double lsqc = std::pow( double(ni)/double(ii), 2 ) + std::pow( double(nj)/double(jj), 2 ) + std::pow( double(nk)/double(kk), 2 );

                if( lsqc < lsq )
                {
                    lsq=lsqc;
                    ni_=ii;
                    nj_=jj;
                    nk_=kk;
                }
            }
        }
    }
            
    if( master() ) 
    {
        std::cout<< "\n Decomposed as: ( " << ni_ << " " << nj_ << " " << nk_ << " )" <<std::endl;
    }


    int testRank=0;
    bool found=false;

    for( int i=0; i<ni_ && !found; i++ )
    {
        for( int j=0; j<nj_ && !found; j++ )
        {
            for( int k=0; k<nk_ && !found; k++ )
            {
                if( myProcNo() == testRank && !found )
                {
                    i_=i;
                    j_=j;
                    k_=k;
                    found = true;
                }
                testRank++;
            }
        }
    }
}

void parallelCom::decompose( std::string fileName )
{
    std::string line;
    std::ifstream fin( fileName );

    int ni, nj, nk;
    
    while(std::getline(fin, line))
    {
        std::stringstream strProcessor(line);
        std::string buf;

        strProcessor>>buf;

        if( buf == std::string( "N:" ) )
        {
            strProcessor>>ni>>nj>>nk;
        }
    }

    parallelCom::decompose( ni, nj, nk );
}


#ifdef HAVE_MUI
void parallelCom::initMUI( std::string fileName, Box bb )
{
    std::ifstream fin( fileName );
    std::string line;

    if( !fin.is_open() )
    {
        return;
    }

    std::vector<std::string> interfaces;
    bool halos = false;

    while(std::getline(fin, line))
    {
        std::stringstream strProcessor(line);
        std::string buf;

        strProcessor>>buf;

        if( buf == std::string( "interface:" ) )
        {
            strProcessor>>buf;
            interfaces.push_back( buf );

            std::getline(fin, line);
            std::stringstream strProcessor2(line);

            strProcessor2>>buf;

            assert( buf == std::string( "vector:" ) );
            scalar x, y, z;
            strProcessor2>>buf>>x>>y>>z;
            muiTrans(interfaces.size()-1) = vector( x, y, z );

            std::getline(fin, line);
           
            parallelCom::muiShapesRecv( interfaces.size()-1 ).push_back
            (
                std::unique_ptr<shape>
                (
                    new Box( line )
                )
            ); 

            std::getline(fin, line);
           
            parallelCom::muiShapesSend( interfaces.size()-1 ).push_back
            (
                std::unique_ptr<shape>
                (
                    new Box( line )
                )
            ); 

        }
    }

    if( interfaces.size() == 0 )
    {
        return;
    }

    muiInterfaces_ = mui::create_uniface<mui::default_config>( settings::zoneName(), interfaces );

    for( int i=0; i<muiInterfaces_.size(); i++ )
    {
        Box bboxSend( (* dynamic_cast<Box*>(parallelCom::muiShapesSend(i)[0].get())) && bb );  
        Box bboxRecv = ( muiShapesRecv(i).size()==1 ) ? ( * dynamic_cast<Box*>(muiShapesRecv(i)[0].get()) && bb ) : Box( vector( 1e10, 1e10, 1e10 ), vector( 1e10, 1e10, 1e10 ) );  

        bboxSend.min()+=muiTrans(i);
        bboxSend.max()+=muiTrans(i);
        bboxRecv.min()+=muiTrans(i);
        bboxRecv.max()+=muiTrans(i);


        muiInterfaces_[i]->announce_recv_span
        (
            -1,
            1e6,
            mui::geometry::box<mui::default_config>
            (
                mui::point3d( bboxRecv.min().x(), bboxRecv.min().y(), bboxRecv.min().z() ), 
                mui::point3d( bboxRecv.max().x(), bboxRecv.max().y(), bboxRecv.max().z() ) 
            )
        );

        muiInterfaces_[i]->announce_send_span
        (
            -1,
            1e6,
            mui::geometry::box<mui::default_config>
            (
                mui::point3d( bboxSend.min().x(), bboxSend.min().y(), bboxSend.min().z() ), 
                mui::point3d( bboxSend.max().x(), bboxSend.max().y(), bboxSend.max().z() ) 
            )
        );

    }
    
}
#endif


