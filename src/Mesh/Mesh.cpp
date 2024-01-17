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
#include "Mesh/Mesh.h"
#include <sstream>


Mesh::Mesh( std::string fileName, Time& runTime )
:
    runTime_( runTime ),
    bbox_( vector( -1e10, -1e10, -1e10 ), vector( 1e10, 1e10, 1e10 ) ),
    bboxHalo_( vector( -1e10, -1e10, -1e10 ), vector( 1e10, 1e10, 1e10 ) ),
    hyperbollic_(false)
{
    std::string line;
    std::ifstream fin;
    fin.open( fileName );

    int ni, nj, nk;
    scalar lx, ly, lz, ox, oy, oz, s1, s2;

    while(std::getline(fin, line))
    {
        std::stringstream strProcessor(line);
        std::string buf;
    
        strProcessor>>buf;

        if( buf == std::string( "origin:" ) )
        {
            strProcessor>>ox>>oy>>oz;
        }
        else if( buf == std::string( "N:" ) )
        {
            strProcessor>>ni>>nj>>nk;
        }
        else if( buf == std::string( "L:" ) )
        {
            strProcessor>>lx>>ly>>lz;
        }
        else if( buf == std::string( "hypz:" ) )
        {
            hyperbollic_ = true;
            strProcessor>>s1>>s2;
        }
    }

    glob_n_[0] = ni;
    glob_n_[1] = nj;
    glob_n_[2] = nk;

    dx_=lx/ni;
    dy_=ly/nj;
    dz_=lz/nk;

    lx_=lx;
    ly_=ly;
    lz_=lz;
    
//for( int i=settings::m()/2; i<ni-settings::m()/2; i++ )
//    {   
//        for( int j=settings::m()/2; j<nj-settings::m()/2; j++ )
//        {   
//            for( int k=settings::m()/2; k<nk()-settings::m()/2; k++ )
//            {   
//                cx = (i-settings::m()/2)*mesh.dx +mesh.origin().x();
//                cy = (j-settings::m()/2)*mesh.dy +mesh.origin().y();
//                cz = (k-settings::m()/2)*mesh.dz +mesh.origin().z();
//	    }
//	}
//    }

    ni_= ni/parallelCom::ni() + 1; 
    nj_= nj/parallelCom::nj() + 1; 
    nk_= nk/parallelCom::nk() + 1; 

    if( (parallelCom::i() == parallelCom::ni()-1) && parallelCom::ni()>1 )
    {
        ni_ = ni + 1 - (ni_-1)*(parallelCom::ni()-1);
    }
    if( (parallelCom::j() == parallelCom::nj()-1) && parallelCom::nj()>1 )
    {
        nj_ = nj + 1 - (nj_-1)*(parallelCom::nj()-1);
    }
    if( (parallelCom::k() == parallelCom::nk()-1) && parallelCom::nk()>1 )
    {
        nk_ = nk + 1 - (nk_-1)*(parallelCom::nk()-1);
    }

    origin_ = vector( (ni_-1)*dx_*parallelCom::i()+ox, (nj_-1)*dy_*parallelCom::j()+oy, (nk_-1)*dz_*parallelCom::k()+oz );
#ifdef HAVE_FFTMPI
    localStart_[0] = origin_[0];
    localStart_[1] = origin_[1];
    localStart_[2] = origin_[2];
#endif
    int ni_copy = ni_;
    int nj_copy = nj_;
    int nk_copy = nk_;
    reduce( ni_copy, minOp<int>() );
    reduce( nj_copy, minOp<int>() );
    reduce( nk_copy, minOp<int>() );

    if( parallelCom::master() )
    {
        std::cout << "Minimum number of cells per direction " << ni_copy << " " << nj_copy << " " << nk_copy << "." << std::endl;
    }

    ni_copy = ni_;
    nj_copy = nj_;
    nk_copy = nk_;
    reduce( ni_copy, maxOp<int>() );
    reduce( nj_copy, maxOp<int>() );
    reduce( nk_copy, maxOp<int>() );

    if( parallelCom::master() )
    {
        std::cout << "Maximum number of cells per direction " << ni_copy << " " << nj_copy << " " << nk_copy << "." << std::endl;
    }

    ni_+= settings::m();
    nj_+= settings::m();
    nk_+= settings::m();
    
    s1_=s1;
    s2_=s2;

    bbox_ = Box( this->loc(settings::m()/2,settings::m()/2,settings::m()/2), this->loc(ni_-settings::m()/2,nj_-settings::m()/2,nk_-settings::m()/2) );
    bboxHalo_ = Box( this->loc(0,0,0), this->loc(ni_,nj_,nk_) );

    //non-uniform spacing in z:
    dzdxi_ = new scalar[nk_];

    for( int k=settings::m()/2; k<nk_-settings::m()/2; k++ )
    {
        dzdxi_[k]=0;
        for( int l=0; l<settings::m()/2; l++ )
        {
            dzdxi_[k] +=
            (
                settings::coef()[l] / 2.0 / (l+1) / dz_ * z( k+1+l )
               -settings::coef()[l] / 2.0 / (l+1) / dz_ * z( k-1-l )
            );
        }
    }

    //non-uniform spacing in z:
    d2zdxi2_ = new scalar[nk_];

    for( int k=settings::m()/2; k<nk_-settings::m()/2; k++ )
    {
        d2zdxi2_[k]=0;
        for( int l=0; l<settings::m()/2; l++ )
        {
            d2zdxi2_[k] +=
            (
                settings::coef()[l] / pow( (l+1) * dz_, 2 ) * z(k+1+l) 
               +settings::coef()[l] / pow( (l+1) * dz_, 2 ) * z(k-1-l) 
               +(
                   -2.0 * settings::coef()[l] / pow( (l+1) * dz_, 2 )
                ) * z(k) //P
            );
        }
    }

    //output mesh details
    scalar max=-1e99;
    for( int k=settings::m()/2; k<nk_-settings::m()/2; k++ )
    {
        if( dzdxi_[k] > max )
        {
            max = dzdxi_[k];
        }
    }

    reduce( max, maxOp<scalar>() );
  
    scalar r1 = -1e99;
    if( parallelCom::k()==0 )
    {
        r1 = dzdxi_[settings::m()/2+1] / dzdxi_[settings::m()/2];
    }

    reduce( r1, maxOp<scalar>() );
  
    scalar r2 = -1e99;
    if( parallelCom::k()==parallelCom::nk()-1 )
    {
        r2 = dzdxi_[nk_-settings::m()/2-2] / dzdxi_[nk_-settings::m()/2-1];
    }

    reduce( r2, maxOp<scalar>() );
 

    if( parallelCom::master() )
    { 
        std::cout<<"Maximum z mesh spacing: " << max * dz_<<std::endl;
        std::cout<<"Growth ratio 1: " << r1 <<std::endl;
        std::cout<<"Growth ratio 2: " << r2 <<std::endl;
    }

#ifdef HAVE_MUI
    parallelCom::initMUI( settings::zoneName() + "/muiInterfaces", bboxHalo_ );
    this->sendMUIPoints();
#endif
}


void Mesh::write( std::string path )
{
    if( !runTime_.write() )
    {
        return;
    }

    std::stringstream fileName;
    fileName<< path << "/geom.p" << boost::format("%|04|.geo") % parallelCom::myProcNo();

	std::ofstream fout( fileName.str(), std::ios::out|std::ios::binary ) ;

    std::ostringstream buf;
	buf << "C Binary";
    buf.str().resize(80);
	fout.write( buf.str().c_str(), 80 ) ;
	
    buf.str("");
    buf.clear();
    buf << "desc1";
    buf.str().resize(80);
	fout.write( buf.str().c_str(), 80 ) ;
		
    buf.str("");
    buf.clear();
    buf << "desc2";
    buf.str().resize(80);
	fout.write( buf.str().c_str(), 80 ) ;
   
    buf.str("");
    buf.clear();
    buf << "node id assign";
    buf.str().resize(80);
	fout.write( buf.str().c_str(), 80 ) ;
	
    buf.str("");
    buf.clear();
    buf << "element id assign";
    buf.str().resize(80);
	fout.write( buf.str().c_str(), 80 ) ;
	
    buf.str("");
    buf.clear();
    buf << "part";
    buf.str().resize(80);
	fout.write( buf.str().c_str(), 80 ) ;
	int partNum = 1 ;
	fout.write( reinterpret_cast<char*>( &partNum ), sizeof( int ) ) ;
	buf.str("");
    buf.clear();
    buf << "part desc";
    buf.str().resize(80);
	fout.write( buf.str().c_str(), 80 ) ;

    buf.str("");
    buf.clear();
    buf << "block rectilinear";
    buf.str().resize(80);
	fout.write( buf.str().c_str(), 80 ) ;

    int nip = ni_-settings::m(), njp=nj_-settings::m(), nkp=nk_-settings::m(); 
	fout.write( reinterpret_cast<char*>( &nip ), sizeof( int ) ) ;
	fout.write( reinterpret_cast<char*>( &njp ), sizeof( int ) ) ;
	fout.write( reinterpret_cast<char*>( &nkp ), sizeof( int ) ) ;

    for( int i=0; i<nip; i++ )
    {
        float x = origin_.x() + dx_*i;
    	fout.write( reinterpret_cast<char*>( &x ), sizeof( float ) ) ;
    }

    for( int j=0; j<njp; j++ )
    {
        float y = origin_.y() + dy_*j;
    	fout.write( reinterpret_cast<char*>( &y ), sizeof( float ) ) ;
    }

    for( int k=0; k<nkp; k++ )
    {
        float z = this->z(k+settings::m()/2);
    	fout.write( reinterpret_cast<char*>( &z ), sizeof( float ) ) ;
    }

    fout.close();

    this->writeCase( path );
}

void Mesh::writeCase( std::string path )
{ 
    std::ofstream fout;
    std::stringstream fileName;
    fileName<< path << "/fd.p" << boost::format("%|04|.case") % parallelCom::myProcNo();

	fout.open( fileName.str() ) ;

	fout << "FORMAT" << std::endl ;
	fout << "type: ensight gold" <<std::endl ;

    fout << "GEOMETRY" << std::endl ;
    fout << "model: geom.p" << boost::format("%|04|.geo") % parallelCom::myProcNo() <<std::endl ;

    fout << "VARIABLE" << std::endl;
    fout << "vector per node: velocity " << boost::format("U.%|04|.*****.dat") % parallelCom::myProcNo() << std::endl ;
    fout << "scalar per node: p  " << boost::format("p.%|04|.*****.dat") % parallelCom::myProcNo() << std::endl ;
    fout << "scalar per node: T  " << boost::format("T.%|04|.*****.dat") % parallelCom::myProcNo() << std::endl ;

	fout << "TIME" << std::endl ;
	fout << "time set: 1" << std::endl ;
	fout << "number of steps: " << runTime_.steps().size() << std::endl ;	
	fout << "filename start number: 00001" << std::endl ;
	fout << "filename increment: 1" << std::endl ;
	fout << "time values: " ;

	for( unsigned i=0 ; i<runTime_.steps().size() ; i++ ) 
    {
		fout << runTime_.steps()[i] << std::endl ;
	}

	fout << std::endl ;

	fout.close() ;

    this->writeRestartInfo( path );
}

void Mesh::writeRestartInfo( std::string path )
{
    std::ofstream fout;
    std::stringstream fileName;
    fileName<< path << "/restart.p" << boost::format(".%|04|.case") % parallelCom::myProcNo();
    fout.open( fileName.str() ) ;

    fout<<runTime_.steps().size()<<std::endl;

    for( unsigned i=0 ; i<runTime_.steps().size() ; i++ )
    {   
        fout << runTime_.steps()[i] << std::endl ;
    }

    fout.close();
}


scalar Mesh::z( int k ) const
{
    scalar z;
    if( hyperbollic_ )
    {
        scalar A = std::pow( s1_/s2_, 0.5 );
#ifdef HAVE_FFTMPI
        int globK = k-settings::m()/2 + localStart_[2];
#else 
        int globK = k-settings::m()/2 + (parallelCom::k())*(nk_-settings::m()-1);
#endif
        int globN = glob_n_[2]+1;
        scalar R = scalar(globK)/scalar(globN-1)-0.5;
        static scalar b = 1.0;
        static bool setb = false;
        
        if( !setb )
        {
            setb = true;
            for( int i=1; i<10000; i++ )
            {
                b = std::asinh( b / ((globN-1)*std::pow(s1_*s2_, 0.5)));
            }
        }

        scalar U = 1.0 + std::tanh( b*R ) / std::tanh( b/2.0 );

        z=U / (2.0*A+(1-A)*U)*lz_;
    }
    else
    {
        z = (k-settings::m()/2)*dz_ + origin_.z();
    }
    return z;
}


#ifdef HAVE_MUI
void Mesh::sendMUIPoints()
{
    for( int inf=0; inf<parallelCom::numInterfaces(); inf++ )
    {
        for( int i=0; i<this->ni(); i++ )
        {
            for( int j=0; j<this->nj(); j++ )
            {
                for( int k=0; k<this->nk(); k++ )
                {
                    vector loc( this->loc(i, j, k) );                       
                    vector locTrans( loc+parallelCom::muiTrans(inf) );
 
                    bool inside=false;
                    for( unsigned s=0; s<parallelCom::muiShapesSend(inf).size(); s++ )
                    {
                        if( parallelCom::muiShapesSend(inf)[s]->inside( loc ) )
                        {
                            inside = true;
                            break;
                        }
                    }
    
                    if( inside )
                    {
                        parallelCom::muiInterface(inf)->push
                        ( 
                            mui::point3d( locTrans.x(), locTrans.y(), locTrans.z() ) 
                        );
                    }
                }
            }
        }
    }
}
#endif
