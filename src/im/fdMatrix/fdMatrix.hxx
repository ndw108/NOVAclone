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
template <class T>
fdMatrix<T>::fdMatrix
( 
    std::shared_ptr<Field<T> > transportedF
)
:
    Field<T>::Field( transportedF->mesh() ),
    transportedFld_( transportedF )
{
    for( int i=0; i<6; i++ )
    {
        reuseTmp<T> alt( transportedF->mesh() );
        auto al( alt() );
        al_.push_back( al ) ;
    }

    reuseTmp<T> apt( transportedF->mesh() );
    ap_ = apt();

    auto& m = transportedF->mesh();

    for( int n=0; n<m.ni()*m.nj()*m.nk(); n++ )
    {
        this->operator()(n) = T(0.0);
    }
}

template <class T>
fdMatrix<T>::fdMatrix
(
    std::shared_ptr<fdMatrix<T> > m
)
{
    *this = (*m);
}
 


//====

template <class T>
void fdMatrix<T>::solve(scalar tol, int maxIts)
{
    int its=0;
    T res;
    while( its < maxIts )
    {
        res = this->residual();
        bool converged = true;
        
        for( int c=0; c< pTraits<T>::nComp; c++ )
        {
            if( res > tol )
            {
                converged = false;
            }
        }

        if( converged ) break;

        its++;


        //todo: outer loop for components of T if this would ever be used for vectors, etc..
        //
        std::vector<scalar> a, b, c, d, x;
        //j sweeps
        a.resize(this->nj());
        b.resize(this->nj());
        c.resize(this->nj());
        d.resize(this->nj());
        x.resize(this->nj());

        for( int i=settings::m()/2; i<this->ni()-settings::m()/2; i++ )
        {
            for( int k=settings::m()/2; k<this->nk()-settings::m()/2; k++ )
            {
                for( int j=0; j<this->nj(); j++ )
                {
                    bool boundary( j<settings::m()/2 || j>=this->nj()-settings::m()/2-1 );

                    a[j] = boundary ? 0.0 : this->al(1)(i,j,k);
                    b[j] = boundary ? 1.0 : this->ap()(i,j,k);
                    c[j] = boundary ? 0.0 : this->al(0)(i,j,k);
                    d[j] = boundary ? transportedFld_->operator()(i,j,k) :
                        this->operator()(i, j, k)
                       -this->al(2)(i, j, k) * transportedFld_->operator()(i+1, j, k)
                       -this->al(3)(i, j, k) * transportedFld_->operator()(i-1, j, k)
                       -this->al(4)(i, j, k) * transportedFld_->operator()(i, j, k+1)
                       -this->al(5)(i, j, k) * transportedFld_->operator()(i, j, k-1);
                }
                
                tools::tdma( a, b, c, d, x );

                for( int j=settings::m()/2; j<this->nj()-settings::m()/2; j++ )
                {
                    transportedFld_->operator()(i, j, k) = x[j];
                }
            }
        }

        //i sweeps
        a.resize(this->ni());
        b.resize(this->ni());
        c.resize(this->ni());
        d.resize(this->ni());
        x.resize(this->ni());

        for( int j=settings::m()/2; j<this->nj()-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<this->nk()-settings::m()/2; k++ )
            {
                for( int i=0; i<this->ni(); i++ )
                {
                    bool boundary( i<settings::m()/2 || i>=this->ni()-settings::m()/2-1 );

                    a[i] = boundary ? T(0.0) : this->al(3)(i,j,k);
                    b[i] = boundary ? T(1.0) : this->ap()(i,j,k);
                    c[i] = boundary ? T(0.0) : this->al(2)(i,j,k);
                    d[i] = boundary ? transportedFld_->operator()(i,j,k) :
                        this->operator()(i, j, k)
                       -this->al(0)(i, j, k) * transportedFld_->operator()(i, j+1, k)
                       -this->al(1)(i, j, k) * transportedFld_->operator()(i, j-1, k)
                       -this->al(4)(i, j, k) * transportedFld_->operator()(i, j, k+1)
                       -this->al(5)(i, j, k) * transportedFld_->operator()(i, j, k-1);
                }
                
                tools::tdma( a, b, c, d, x );

                for( int i=settings::m()/2; i<this->ni()-settings::m()/2; i++ )
                {
                    transportedFld_->operator()(i, j, k) = x[i];
                }
            }
        }

        //k sweeps
        a.resize(this->nk());
        b.resize(this->nk());
        c.resize(this->nk());
        d.resize(this->nk());
        x.resize(this->nk());

        for( int i=settings::m()/2; i<this->ni()-settings::m()/2; i++ )
        {
            for( int j=settings::m()/2; j<this->nj()-settings::m()/2; j++ )
            {
                for( int k=0; k<this->nk(); k++ )
                {
                    bool boundary( k<settings::m()/2 || k>=this->nk()-settings::m()/2-1 );

                    a[k] = boundary ? 0.0 : this->al(5)(i,j,k);
                    b[k] = boundary ? 1.0 : this->ap()(i,j,k);
                    c[k] = boundary ? 0.0 : this->al(4)(i,j,k);
                    d[k] = boundary ? transportedFld_->operator()(i,j,k) :
                        this->operator()(i, j, k)
                       -this->al(2)(i, j, k) * transportedFld_->operator()(i+1, j, k)
                       -this->al(3)(i, j, k) * transportedFld_->operator()(i-1, j, k)
                       -this->al(0)(i, j, k) * transportedFld_->operator()(i, j+1, k)
                       -this->al(1)(i, j, k) * transportedFld_->operator()(i, j-1, k);
                }
                
                tools::tdma( a, b, c, d, x );

                for( int k=settings::m()/2; k<this->nk()-settings::m()/2; k++ )
                {
                    transportedFld_->operator()(i, j, k) = x[k];
                }
            }
        }
        
        transportedFld_->correctBoundaryConditions();
    }

    if( parallelCom::master() && this->mesh().time().writelog() )
    {
        std::cout<<"Poisson Equation converged in "<<its<<" iterations."<<std::endl;
        std::cout<<"Residual: "<<res<<"."<<std::endl;
    }
}


template <class T>
void fdMatrix<T>::deferredCorrection
(
    const std::shared_ptr<Field<T> > hoa
)
{
    for( int i=settings::m()/2; i<this->ni()-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<this->nj()-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<this->nk()-settings::m()/2; k++ )
            {
                this->operator()(i, j, k) += //rhs
               -hoa->operator()(i, j, k)
               +this->al(0)(i, j, k) * transportedFld_->operator()(i, j+1, k)
               +this->al(1)(i, j, k) * transportedFld_->operator()(i, j-1, k)
               +this->al(2)(i, j, k) * transportedFld_->operator()(i+1, j, k)
               +this->al(3)(i, j, k) * transportedFld_->operator()(i-1, j, k)
               +this->al(4)(i, j, k) * transportedFld_->operator()(i, j, k+1)
               +this->al(5)(i, j, k) * transportedFld_->operator()(i, j, k-1)
               +this->ap()(i, j, k) * transportedFld_->operator()(i, j, k);
            }
        }
    } 
}

template <class T>
void fdMatrix<T>::relax( const scalar& urf )
{
    (*dynamic_cast<Field<T>* >(this)) += (1.0-urf)/urf * ( this->ap() * transportedFld_ );
    this->ap() /= urf;
}


template <class T>
void fdMatrix<T>::setReference
(
    int i,
    int j,
    int k,
    scalar val
)
{
    this->operator()(i, j, k) += 1e15 * val;
    this->ap()(i, j, k) *= 1e15;
}
    

template <class T>
T fdMatrix<T>::residual()
{
    T res(0.0);

    for( int i=settings::m()/2; i<this->ni()-settings::m()/2; i++ )
    {   
        for( int j=settings::m()/2; j<this->nj()-settings::m()/2; j++ )
        {   
            for( int k=settings::m()/2; k<this->nk()-settings::m()/2; k++ )
            {   
                T resi =
                (
                    this->operator()(i, j, k)
                   -this->al(0)(i, j, k) * transportedFld_->operator()(i, j+1, k)
                   -this->al(1)(i, j, k) * transportedFld_->operator()(i, j-1, k)
                   -this->al(2)(i, j, k) * transportedFld_->operator()(i+1, j, k)
                   -this->al(3)(i, j, k) * transportedFld_->operator()(i-1, j, k)
                   -this->al(4)(i, j, k) * transportedFld_->operator()(i, j, k+1)
                   -this->al(5)(i, j, k) * transportedFld_->operator()(i, j, k-1)
                   -this->ap()(i, j, k) * transportedFld_->operator()(i, j, k)
                );

                for( int c=0; c<pTraits<T>::nComp; c++ )
                {
                    component(res, c) += pow( component( resi, c ), 2 );
                }
            }
        }
    }

    int n = (this->ni()-settings::m())*(this->nj()-settings::m())*(this->nk()-settings::m());
    all_reduce( res, plusOp<T>() );
    all_reduce( n, plusOp<int>() );

    for( int c=0; c< pTraits<T>::nComp; c++ )
    {
        component(res,c) = pow( component(res,c)/( n ), 0.5 );
    }

    return res;
}

//====



template <class T>
fdMatrix<T>& fdMatrix<T>::operator=( const fdMatrix<T>& r )
{
    Field<T>::operator=( dynamic_cast<const Field<T>& >( r ) );

    for( int l=0; l<6; l++ )
    {
        this->al(l) = r.al(l);
    }

    this->ap() = r.ap(); 
    
    return *this;
}

template <class T>
fdMatrix<T>& fdMatrix<T>::operator+=( const fdMatrix<T>& r )
{
    Field<T>::operator+=( dynamic_cast<const Field<T>& >( r ) );

    for( int l=0; l<6; l++ )
    {
        this->al(l) += r.al(l);
    }

    this->ap() += r.ap();

    return *this;
}

template <class T>
fdMatrix<T>& fdMatrix<T>::operator-=( const fdMatrix<T>& r )
{
    Field<T>::operator-=( dynamic_cast<const Field<T>& >( r ) );

    for( int l=0; l<6; l++ )
    {
        this->al(l) -= r.al(l);
    }

    this->ap() -= r.ap();

    return *this;
}

template <class T>
fdMatrix<T>& fdMatrix<T>::operator*=( const fdMatrix<T>& r )
{
    Field<T>::operator*=( dynamic_cast<const Field<T>& >( r ) );

    for( int l=0; l<6; l++ )
    {
        this->al(l) *= r.al(l);
    }

    this->ap() *= r.ap();

    return *this;
}

//-------
//global operators

template <class T> 
std::shared_ptr<fdMatrix<T> > operator+
(
    const std::shared_ptr<fdMatrix<T> > t1,
    const std::shared_ptr<fdMatrix<T> > t2
)
{
    std::shared_ptr<fdMatrix<T> > ptr( t1 );
    *ptr += *t2;
    return ptr;
}

template <class T> 
std::shared_ptr<fdMatrix<T> > operator==
(
    const std::shared_ptr<fdMatrix<T> > t1,
    const std::shared_ptr<Field<T> > t2
)
{
    std::shared_ptr<fdMatrix<T> > ptr( t1 );
    (*dynamic_cast<Field<T>* >(ptr.get())) += *t2;
    return ptr;
}

template <class T> 
std::shared_ptr<fdMatrix<T> > operator-
(
    const std::shared_ptr<fdMatrix<T> > t1,
    const std::shared_ptr<fdMatrix<T> > t2
)
{
    std::shared_ptr<fdMatrix<T> > ptr( t1 );
    *ptr -= *t2;
    return ptr;
}

template <class T> 
std::shared_ptr<fdMatrix<T> > operator-
(
    const std::shared_ptr<fdMatrix<T> > t
)
{
    std::shared_ptr<fdMatrix<T> > ptr( t );

    for( int l=0; l<6; l++ )
    {   
        ptr->al(l) *= -1;
    }

    ptr->ap() *= -1;

    (*dynamic_cast<Field<T>* >(ptr.get())) *= -1;

    return ptr;
}

