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
#include <cmath>
#include <iterator>
template <class T>
Field<T>::Field
( 
    Mesh& m
)
: 
m_(m),
curTime_(-1)
{
    void* mem = boost::alignment::aligned_alloc(64, sizeof(T)*ni()*nj()*nk());

    for(std::size_t i=0; i!=unsigned(ni()*nj()*nk()); ++i)
    {
        new (static_cast<T*>(mem) + i ) T;
    }

    v_=static_cast<T*>(mem); 
}    

template <class T>
Field<T>::Field
( 
    Mesh& m,
    const T& val
)
:
Field<T>( m )
{
    BOOST_ALIGN_ASSUME_ALIGNED(v_, 64); 

    #pragma omp parallel for simd
    for( int i=0; i<ni()*nj()*nk(); i++ )
    {
        v_[i] = val;
    } 
}

template <class T>
Field<T>::Field
( 
    Mesh& m,
    std::string FieldName
)
: 
m_(m),
FieldName_(FieldName),
curTime_(-1)
{
    void* mem = boost::alignment::aligned_alloc(64, sizeof(T)*ni()*nj()*nk());

    for(std::size_t i=0; i!=unsigned(ni()*nj()*nk()); ++i)
    {
        new (static_cast<T*>(mem) + i ) T;
    }

    v_=static_cast<T*>(mem); 

    if( !FieldName.empty() ) 
    {
        this->setBC( FieldName + "BC" );
    }
}

template <class T>
Field<T>::Field
( 
    Mesh& m,
    const T& val,
    std::string FieldName
)
:
Field<T>( m, FieldName )
{
    BOOST_ALIGN_ASSUME_ALIGNED(v_, 64); 

    #pragma omp parallel for simd
    for( int i=0; i<ni()*nj()*nk(); i++ )
    {
        v_[i] = val;
    } 
}

template <class T>
Field<T>::Field
( 
    Field<T>& fld
)
:
Field<T>( fld.m_, fld.FieldName_ )
{
    #pragma omp parallel for collapse(3)
    for( int i=0; i<ni(); i++ )
    {
        for( int j=0; j<nj(); j++ )
        {
            for( int k=0; k<nk(); k++ )
            {
                v_[k+j*nk()+i*(nk()*nj())] = fld(i, j, k);
            }
        }
    }
}


template <class T>
Field<T>::Field
( 
    std::shared_ptr<Field<T> > tf
)
:
m_(tf->mesh())
{
    *this = (*tf);
}



//implementation
//=====

template <class T>
void Field<T>::operator=(const Field<T>& fld)
{
    int nn = ni()*nj()*nk();
    T* ptr = fld.v_;

    BOOST_ALIGN_ASSUME_ALIGNED(ptr, 64); 

    #pragma omp parallel for simd
    for( int i=0; i<nn; ++i )
    {
        v_[i] = ptr[i];
    }
}


template <class T>
void Field<T>::operator=(const scalar val)
{
    BOOST_ALIGN_ASSUME_ALIGNED(v_, 64); 

    #pragma omp parallel for simd
    for( int i=0; i<ni()*nj()*nk(); ++i )
    {
        v_[i] = val;
    }
}

template <class T>
void Field<T>::operator=(const std::shared_ptr<Field<T> > fld)
{
    *this = (*fld);
}

template <class T>
void Field<T>::operator+=( const Field<T>& fld)
{ 
    T* ptr = fld.v_;
    
    BOOST_ALIGN_ASSUME_ALIGNED(v_, 64); 

    #pragma omp parallel for simd 
    for( int i=0; i<ni()*nj()*nk(); i++ )
    {
        v_[i] += ptr[i];
    }
}


template <class T>
void Field<T>::operator+=( const scalar val)
{ 
    BOOST_ALIGN_ASSUME_ALIGNED(v_, 64); 

    #pragma omp parallel for simd
    for( int i=0; i<ni()*nj()*nk(); i++ )
    {
        v_[i] += val;
    }
}

template <class T>
void Field<T>::operator+=( const std::shared_ptr<Field<T> > fld)
{ 
    T* ptr = (*fld).v_;
    BOOST_ALIGN_ASSUME_ALIGNED(v_, 64); 

    #pragma omp parallel for simd
    for( int i=0; i<ni()*nj()*nk(); i++ )
    {
        v_[i] += ptr[i];
    }
}

template <class T>
void Field<T>::operator-=( const Field<T>& fld)
{
    T* ptr = fld.v_;
    BOOST_ALIGN_ASSUME_ALIGNED(v_, 64); 

    
    #pragma omp parallel for simd
    for( int i=0; i<ni()*nj()*nk(); i++ )
    {
        v_[i] -= ptr[i];
    }
}

template <class T>
void Field<T>::operator-=( const std::shared_ptr<Field<T> > fld)
{
    T* ptr = (*fld).v_;
    BOOST_ALIGN_ASSUME_ALIGNED(v_, 64); 

    #pragma omp parallel for simd 
    for( int i=0; i<ni()*nj()*nk(); i++ )
    {
        v_[i] -= ptr[i];
    }
}

template <class T>
void Field<T>::operator-=( const scalar val)
{
    BOOST_ALIGN_ASSUME_ALIGNED(v_, 64); 

    #pragma omp parallel for simd
    for( int i=0; i<ni()*nj()*nk(); i++ )
    {
        v_[i] -= val;
    }
}

template <class T>
void Field<T>::operator*=( const Field<T>& fld)
{
    T* ptr = fld.v_;

    BOOST_ALIGN_ASSUME_ALIGNED(v_, 64); 

    #pragma omp parallel for simd
    for( int i=0; i<ni()*nj()*nk(); i++ )
    {
        v_[i] *= ptr[i];
    }
}

template <class T>
void Field<T>::operator*=( const scalar v )
{
    BOOST_ALIGN_ASSUME_ALIGNED(v_, 64); 

    #pragma omp parallel for simd
    for( int i=0; i<ni()*nj()*nk(); i++ )
    {   
        v_[i] *= v;
    }
}

template <class T>
void Field<T>::operator*=( const std::shared_ptr<Field<T> > fld)
{
    T* ptr = (*fld).v_;

    BOOST_ALIGN_ASSUME_ALIGNED(v_, 64); 
  
    #pragma omp parallel for simd
    for( int i=0; i<ni()*nj()*nk(); i++ )
    {
        v_[i] *= ptr[i];
    }
}


template <class T>
void Field<T>::operator/=( const scalar v )
{
    BOOST_ALIGN_ASSUME_ALIGNED(v_, 64); 

    #pragma omp parallel for simd
    for( int i=0; i<ni()*nj()*nk(); i++ )
    {   
        v_[i] /= v;
    }
}




template <class T>
void Field<T>::correctBoundaryConditions()
{
    this->updateMUI();

    //all non mpi/mui BCs first
    for( int d=0; d<6*pTraits<T>::nComp; d++ )
    {
        if( bc_[d] && dynamic_cast<muiBC<T>*>( bc_[d].get() ) == NULL && dynamic_cast<mpiBC<T>*>( bc_[d].get() ) == NULL )
        {
            bc_[d]->update(this);
        }
    }

    for( int d=0; d<6*pTraits<T>::nComp; d++ )
    {
        if(bc_[d]->dir() == east || bc_[d]->dir() == west )
        {
            if( bc_[d] && dynamic_cast<mpiBC<T>*>( bc_[d].get() ) != NULL )
            {
                bc_[d]->initUpdate(this);
            }
        }
    }

    for( int d=0; d<6*pTraits<T>::nComp; d++ )
    {
        if(bc_[d]->dir() == east || bc_[d]->dir() == west )
        {
            if( bc_[d] && dynamic_cast<mpiBC<T>*>( bc_[d].get() ) != NULL )
            {
                bc_[d]->update(this);
            }
        }
    }


    for( int d=0; d<6*pTraits<T>::nComp; d++ )
    {
        if(bc_[d]->dir() == north || bc_[d]->dir() == south )
        {
            if( bc_[d] && dynamic_cast<mpiBC<T>*>( bc_[d].get() ) != NULL )
            {
                bc_[d]->initUpdate(this);
            }
        }
    }

    for( int d=0; d<6*pTraits<T>::nComp; d++ )
    {
        if(bc_[d]->dir() == south || bc_[d]->dir() == north )
        {
            if( bc_[d] && dynamic_cast<mpiBC<T>*>( bc_[d].get() ) != NULL )
            {
                bc_[d]->update(this);
            }
        }
    }


    for( int d=0; d<6*pTraits<T>::nComp; d++ )
    {
        if(bc_[d]->dir() == top || bc_[d]->dir() == bottom )
        {
            if( bc_[d] && dynamic_cast<mpiBC<T>*>( bc_[d].get() ) != NULL )
            {
                bc_[d]->initUpdate(this);
            }
        }
    }

    for( int d=0; d<6*pTraits<T>::nComp; d++ )
    {
        if(bc_[d]->dir() == top || bc_[d]->dir() == bottom )
        {
            if( bc_[d] && dynamic_cast<mpiBC<T>*>( bc_[d].get() ) != NULL )
            {
                bc_[d]->update(this);
            }
        }
    }


    for( int d=0; d<6*pTraits<T>::nComp; d++ )
    {
        if( bc_[d] && dynamic_cast<muiBC<T>*>( bc_[d].get() ) != NULL )
        {
            bc_[d]->update(this);
        }
    }

}

//edit to hdf5
template <class T>
void Field<T>::write( const std::string path, const std::string name )
{
    if( !this->mesh().time().write() )
    {
        return;
    }

    std::stringstream fileName;
    fileName<< path << "/" << name << boost::format(".%|04|.%|05|.dat") % parallelCom::myProcNo() % this->mesh().time().fileNum();

    std::ofstream fout;
    fout.open( fileName.str(), std::ios::out|std::ios::binary );

    std::ostringstream buf;
    buf << "description";
    buf.str().resize(80);
    fout.write( buf.str().c_str(), 80 );

    //part header
    buf.str("");
    buf.clear();
    buf << "part";
    buf.str().resize(80);
    fout.write( buf.str().c_str(), 80 );
    int partNum = 1;
    fout.write( reinterpret_cast<char*>( &partNum ), sizeof( int ) );
    buf.str("");
    buf.clear();
    buf << "point";
    buf.str().resize(80);
    fout.write( buf.str().c_str(), 80 );


    //data
    for( int c=0; c<pTraits<T>::nComp; c++ )
    {
        for( int k=settings::m()/2; k<this->nk()-settings::m()/2; k++ )
        {
            for( int j=settings::m()/2; j<this->nj()-settings::m()/2; j++ )
            {
                for( int i=settings::m()/2; i<this->ni()-settings::m()/2; i++ ) 
                {
                    float v = component(this->operator()(i, j, k), c);
                    fout.write( reinterpret_cast<char*>( &v ), sizeof( float ) );
                }
            }
        }
    }

    fout.close();
}


template <class T>
void Field<T>::read( std::string path, std::string name )
{
    std::stringstream fileName;
    fileName<< path << "/" << name << boost::format(".%|04|.%|05|.dat") % parallelCom::myProcNo() % this->mesh().time().fileNum();

    std::ifstream fin;
    fin.open( fileName.str(), std::ios::in|std::ios::binary ) ;

    fin.seekg(240+sizeof(int),std::ios::beg);

    for( int c=0; c<pTraits<T>::nComp; c++ )
    {
        for( int k=settings::m()/2; k<this->nk()-settings::m()/2; k++ )
        {
            for( int j=settings::m()/2; j<this->nj()-settings::m()/2; j++ )
            {
                for( int i=settings::m()/2; i<this->ni()-settings::m()/2; i++ ) 
                {
                    float v;
                    fin.read( reinterpret_cast<char*>(&v), sizeof( float ) );
                    component(this->operator()(i, j, k), c) = scalar(v);
                }
            }
        }
    }

    fin.close() ;
}

template <class T>
void Field<T>::read( std::string path, std::string name, int num )
{
    std::stringstream fileName;
    fileName<< path<<"/"<<name<<boost::format(".%|04|.%|05|.dat") % name % parallelCom::myProcNo() % num;

    std::ifstream fin;
    fin.open( fileName.str(), std::ios::in|std::ios::binary ) ;

    fin.seekg(240+sizeof(int),std::ios::beg);

    for( int c=0; c<pTraits<T>::nComp; c++ )
    {
        for( int k=settings::m()/2; k<this->nk()-settings::m()/2; k++ )
        {
            for( int j=settings::m()/2; j<this->nj()-settings::m()/2; j++ )
            {
                for( int i=settings::m()/2; i<this->ni()-settings::m()/2; i++ ) 
                {
                    float v;
                    fin.read( reinterpret_cast<char*>(&v), sizeof( float ) );
                    component(this->operator()(i, j, k), c) = scalar(v);
                }
            }
        }
    }

    fin.close() ;
}


template <class T>
void Field<T>::setBC( std::string file )
{
    std::string line;
    std::ifstream fin;
    fin.open( settings::zoneName() + "/" + file );

    if( !fin.is_open() || file == std::string( "" ) )
    {
        return;
    }

    while(std::getline(fin, line))
    {
        std::stringstream strProcessor(line);
        std::string dir, type;
        int comp;
        scalar val;

        strProcessor>>dir>>type>>val>>comp;

        if( parallelCom::internal( BCdirmap[dir], BCcodemap[type] ) )
        {
            type = "mpi";
        } 
        
        bc_[BCdirmap[dir] + comp*6] = BCmaker<T>::makeBC[BCcodemap[type]]( line );        
    }
}


template <class T>
void Field<T>::bound
(
    T min,
    T max
)
{
    for( int c=0; c<pTraits<T>::nComp; c++ )
    {
        for( int k=settings::m()/2; k<this->nk()-settings::m()/2; k++ )
        {
            for( int j=settings::m()/2; j<this->nj()-settings::m()/2; j++ )
            {
                for( int i=settings::m()/2; i<this->ni()-settings::m()/2; i++ )
                {
                    component(this->operator()(i, j, k), c) = 
                    std::min( component(this->operator()(i, j, k), c), component( max, c ) );
                    
                    component(this->operator()(i, j, k), c) = 
                    std::max( component(this->operator()(i, j, k), c), component( min, c ) );
                }
            }
        }
    }
}


template <class T>
T Field<T>::weightedAverage
(
)
{
    T avg(0);
    int n=0;

    for( int c=0; c<pTraits<T>::nComp; c++ )
    {
        for( int k=settings::m()/2; k<this->nk()-settings::m()/2; k++ )
        {
            for( int j=settings::m()/2; j<this->nj()-settings::m()/2; j++ )
            {
                for( int i=settings::m()/2; i<this->ni()-settings::m()/2; i++ )
                {
                    component( avg, c ) += component(this->operator()(i, j, k), c);
                    n++;
                }
            }
        }
    }

    all_reduce( avg, plusOp<T>() );
    all_reduce( n, plusOp<int>() );

    return avg/n;
}
         

//destructor
//=====
template <class T>
Field<T>::~Field()
{
    boost::alignment::aligned_free(v_);
}



//global operations
template <class T>
std::shared_ptr<Field<T> > operator+
( 
    const Field<T>& f1,
    const Field<T>& f2
)
{ 
    reuseTmp<T> res( f1.mesh() );
    T* ptr = (*res()).ptr();
    T* ptr1 = f1.ptr();
    T* ptr2 = f2.ptr();
    int nn = f1.ni()*f1.nj()*f1.nk();

    #pragma omp parallel for simd 
    for( int i=0; i<nn; ++i )
    {
        ptr[i] = ptr1[i] + ptr2[i];
    }
    
    return res();
}

template <class T>
std::shared_ptr<Field<T> > operator+
( 
    const Field<T>& f1,
    const std::shared_ptr<Field<T> > f2
)
{ 
    reuseTmp<T> res( f1.mesh() );
    T* ptr = (*res()).ptr();
    T* ptr1 = f1.ptr();
    T* ptr2 = f2->ptr();
    int nn = f1.ni()*f1.nj()*f1.nk();

    #pragma omp parallel for simd 
    for( int i=0; i<nn; ++i )
    {
        ptr[i] = ptr1[i] + ptr2[i];
    }
    
    return res();
}

template <class T>
std::shared_ptr<Field<T> > operator-
( 
    const std::shared_ptr<Field<T> > f
)
{
    reuseTmp<T> res( f->mesh() );
    T* ptr = res()->ptr();
    T* ptrf = f->ptr();
    int nn = f->ni()*f->nj()*f->nk();

    BOOST_ALIGN_ASSUME_ALIGNED(ptr, 64);
    BOOST_ALIGN_ASSUME_ALIGNED(ptrf, 64); 

    #pragma omp parallel for simd 
    for( int i=0; i<nn; ++i )
    {
        ptr[i] = -ptrf[i];
    }
    
    return res();
}

template <class T>
std::shared_ptr<Field<T> > operator-
( 
    const Field<T>& f1,
    const Field<T>& f2
)
{
    reuseTmp<T> res( f1.mesh() );
    T* ptr = res()->ptr();
    T* ptr1 = f1.ptr();
    T* ptr2 = f2.ptr();
    int nn = f1.ni()*f1.nj()*f1.nk();

    BOOST_ALIGN_ASSUME_ALIGNED(ptr, 64);
    BOOST_ALIGN_ASSUME_ALIGNED(ptr1, 64); 
    BOOST_ALIGN_ASSUME_ALIGNED(ptr2, 64); 

    #pragma omp parallel for simd 
    for( int i=0; i<nn; ++i )
    {   
        ptr[i] = ptr1[i] - ptr2[i];
    }

    return res();
}

template <class T>
std::shared_ptr<Field<T> > operator-
( 
    const Field<T>& f,
    const T& s
)
{
    reuseTmp<T> res( f.mesh() );
    T* ptr = res()->ptr();
    T* ptrf = f.ptr();

    int nn = f.ni()*f.nj()*f.nk();

    #pragma omp parallel for simd 
    for( int i=0; i<nn; ++i )
    {   
        ptr[i] = ptrf[i] - s;
    }

    return res();
}

template <class T>
std::shared_ptr<Field<T> > operator-
(
    const std::shared_ptr<Field<T> > f,
    const T& s
)
{
    return (*f)-s;
}

template <class T>
std::shared_ptr<Field<T> > operator-
( 
    const Field<T>& f1,
    const std::shared_ptr<Field<T> > f2
)
{
    return f1-(*f2);
}

template <class T>
std::shared_ptr<Field<T> > operator-
( 
    const std::shared_ptr<Field<T> > f1,
    const Field<T>& f2
)
{
    return (*f1)-f2;
}

template <class T>
std::shared_ptr<Field<T> > operator-
( 
    const std::shared_ptr<Field<T> > f1,
    const std::shared_ptr<Field<T> > f2
)
{
    return (*f1)-(*f2);
}

template <class T>
std::shared_ptr<Field<T> > operator+
( 
    const std::shared_ptr<Field<T> > f1,
    const std::shared_ptr<Field<T> > f2
)
{
    reuseTmp<T> res( f1->mesh() );
    T* ptr = (*res()).ptr();
    T* ptr1 = f1->ptr();
    T* ptr2 = f2->ptr();
    int nn = f1->ni()*f1->nj()*f1->nk();

    BOOST_ALIGN_ASSUME_ALIGNED(ptr, 64); 

    #pragma omp parallel for simd 
    for( int i=0; i<nn; ++i )
    {
        ptr[i] = ptr1[i] + ptr2[i];
    }
    
    return res();
}

template <class T>
std::shared_ptr<Field<T> > operator*
( 
    scalar d,
    const Field<T>& f
)
{
    reuseTmp<T> res( f.mesh() );
    T* ptr = (*res()).ptr();
    T* ptrf = f.ptr();
    int nn = f.ni()*f.nj()*f.nk();

    #pragma omp parallel for simd 
    for( int i=0; i<nn; ++i )
    {
        ptr[i] = ptrf[i] * d;
    }
    
    return res();
}


template <class T1, class T2>
std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > operator*
( 
    const Field<T1>& f1,
    const Field<T2>& f2
)
{
    reuseTmp<typename outerProductType<T1, T2>::type> rest( f1.mesh() );
    std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > res( rest() );

    typename outerProductType<T1, T2>::type* ptr = res->ptr();
    T1* ptr1 = f1.ptr();
    T2* ptr2 = f2.ptr();

    int nn = f1.ni()*f1.nj()*f1.nk();

    #pragma omp parallel for simd 
    for( int i=0; i<nn; i++ )
    {
        ptr[i] = ptr1[i] * ptr2[i];
    }

    return res;
}


template <class T1, class T2>
std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > operator*
( 
    const std::shared_ptr<Field<T1> >& f1,
    const Field<T2>& f2
)
{
    reuseTmp<typename outerProductType<T1, T2>::type> rest( f1->mesh() );
    std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > res( rest() );

    typename outerProductType<T1, T2>::type* ptr = res->ptr();
    T1* ptr1 = f1->ptr();
    T2* ptr2 = f2.ptr();

    int nn = f1->ni()*f1->nj()*f1->nk();

    #pragma omp parallel for simd 
    for( int i=0; i<nn; i++ )
    {
        ptr[i] = ptr1[i] * ptr2[i];
    }

    return res;
}

template <class T1, class T2>
std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > operator*
( 
    const Field<T1>& f1,
    const std::shared_ptr<Field<T2> >& f2
)
{
    return f2*f1;
}

template <class T1, class T2>
std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > operator*
( 
    const Field<T1>& f,
    const T2& v
)
{
    reuseTmp<typename outerProductType<T1, T2>::type> rest( f.mesh() );
    std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > res( rest() );

    typename outerProductType<T1, T2>::type* ptr = res->ptr();
    T1* ptr1 = f.ptr();

    int nn = f.ni()*f.nj()*f.nk();

    #pragma omp parallel for simd 
    for( int i=0; i<nn; i++ )
    {   
        ptr[i] = ptr1[i] * v;
    }

    return res;
}


template <class T1, class T2>
std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > operator*
( 
    const std::shared_ptr<Field<T1> >& f,
    const T2& v
)
{
    return (*f) * v;
}


template <class T1, class T2>
std::shared_ptr<Field<typename innerProductType<T1, T2>::type> > operator&
( 
    const Field<T1>& f1,
    const Field<T2>& f2
)
{
    reuseTmp<typename innerProductType<T1, T2>::type> rest( f1.mesh() );
    std::shared_ptr<Field<typename innerProductType<T1, T2>::type> > res( rest() );

    typename innerProductType<T1, T2>::type* ptr = res->ptr();
    T1* ptr1 = f1.ptr();
    T2* ptr2 = f2.ptr();

    int nn = f1.ni()*f1.nj()*f1.nk();

    #pragma omp parallel for simd 
    for( int i=0; i<nn; i++ )
    {   
        ptr[i] = ptr1[i] & ptr2[i];
    }

    return res;
}

template <class T1, class T2>
std::shared_ptr<Field<typename innerProductType<T1, T2>::type> > operator&
( 
    const Field<T1>& f1,
    const std::shared_ptr<Field<T2> >& f2
)
{
    return f1 & (*f2);
}

template <class T1, class T2>
std::shared_ptr<Field<typename innerProductType<T1, T2>::type> > operator&
( 
    const std::shared_ptr<Field<T1> >& f1,
    const std::shared_ptr<Field<T2> >& f2
)
{
    return (*f1) & (*f2);
}

template <class T1, class T2>
std::shared_ptr<Field<typename innerProductType<T1, T2>::type> > operator&
( 
    const T1& val,
    const Field<T2>& f
)
{
    reuseTmp<typename innerProductType<T1, T2>::type> rest( f.mesh() );
    std::shared_ptr<Field<typename innerProductType<T1, T2>::type> > res( rest() );

    typename innerProductType<T1, T2>::type* ptr = res->ptr();
    T2* ptr2 = f.ptr();

    int nn = f.ni()*f.nj()*f.nk();

    #pragma omp parallel for simd 
    for( int i=0; i<nn; i++ )
    {   
        ptr[i] = val & ptr2[i];
    }

    return res;
}


template <class T1, class T2>
std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > operator/
( 
    const Field<T1>& f1,
    const Field<T2>& f2
)
{
    reuseTmp<typename outerProductType<T1, T2>::type> rest( f1.mesh() );
    std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > res( rest() );

    typename outerProductType<T1, T2>::type* ptr = res->ptr();
    T1* ptr1 = f1.ptr();
    T2* ptr2 = f2.ptr();

    int nn = f1.ni()*f1.nj()*f1.nk();

    #pragma omp parallel for simd 
    for( int i=0; i<nn; i++ )
    {
        ptr[i] = ptr1[i] / ptr2[i];
    }

    return res;
}

template <class T1, class T2>
std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > operator/
( 
    const Field<T1>& f,
    const T2& v
)
{
    reuseTmp<typename outerProductType<T1, T2>::type> rest( f.mesh() );
    std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > res( rest() );

    typename outerProductType<T1, T2>::type* ptr = res->ptr();
    T1* ptr1 = f.ptr();

    int nn = f.ni()*f.nj()*f.nk();

    #pragma omp parallel for simd 
    for( int i=0; i<nn; i++ )
    {
        ptr[i] = ptr1[i] / v;
    }

    return res;
}

template <class T1, class T2>
std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > operator/
( 
    const std::shared_ptr<Field<T1> >& f,
    const T2& v
)
{
    return (*f)/v;
}

template <class T1, class T2>
std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > operator/
( 
    const std::shared_ptr<Field<T1> > f1,
    const Field<T2>& f2
)
{
    return std::shared_ptr<Field<typename outerProductType<T1, T2>::type> >( *f1 / f2 );
} 

template <class T1, class T2>
std::shared_ptr<Field<typename outerProductType<T1, T2>::type> > operator/
( 
    const std::shared_ptr<Field<T1> > f1,
    const std::shared_ptr<Field<T2> > f2
)
{
    return std::shared_ptr<Field<typename outerProductType<T1, T2>::type> >( (*f1) / (*f2) );
}

template <class T>
std::shared_ptr<Field<T> > operator*
( 
    const scalar d,
    const std::shared_ptr<Field<T> >& f
)
{
    reuseTmp<T> res( f->mesh() );
    T* ptr = (*res()).ptr();
    T* ptrf = f->ptr();
    int nn = f->ni()*f->nj()*f->nk();

    BOOST_ALIGN_ASSUME_ALIGNED(ptr, 64); 

    #pragma omp parallel for simd 
    for( int i=0; i<nn; ++i )
    {
        ptr[i] = ptrf[i] * d;
    }

    return res();
}

template <class T>
std::shared_ptr<Field<T> > operator/
( 
    const Field<T>& f,
    scalar d
)
{
    reuseTmp<T> res( f.mesh() );
    T* ptr = (*res()).ptr();
    T* ptrf = f.ptr();
    int nn = f.ni()*f.nj()*f.nk();

    #pragma omp parallel for simd 
    for( int i=0; i<nn; ++i )
    {
        ptr[i] = ptrf[i] / d;
    }

    return res();
}


template <class T>
std::shared_ptr<Field<T> > transpose
(
    std::shared_ptr<Field<T> > fld
)
{
    reuseTmp<T> res( fld->mesh() );
    T* ptr = (*res()).ptr();
    T* f_ptr = fld->ptr();

    int nn = fld->ni()*fld->nj()*fld->nk();

    #pragma omp parallel for simd 
    for( int i=0; i<nn; i++ )
    {
        ptr[i] = transpose( f_ptr[i] );
    }
   
    return res();
}

template <class T>
std::shared_ptr<Field<T> > dev2
(
    std::shared_ptr<Field<T> > fld
)
{
    reuseTmp<T> res( fld->mesh() );
    T* ptr = (*res()).ptr();
    T* f_ptr = fld->ptr();

    int nn = fld->ni()*fld->nj()*fld->nk();

    #pragma omp parallel for simd 
    for( int i=0; i<nn; i++ )
    {
       ptr[i] = dev2( f_ptr[i] );
    }
    
    return res();
}


