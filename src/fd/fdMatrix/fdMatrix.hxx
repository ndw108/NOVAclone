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
    Field<T>& transportedF
)
:
    Field<T>::Field( transportedF.mesh(), 0.0 ),
    al_(6, Field<double>( transportedF.mesh(), 0.0 )),
    ap_( transportedF.mesh(), 0.0 ),
    transportedFld_( transportedF )
{
}

template <class T>
fdMatrix<T>::fdMatrix
( 
    const fdMatrix<T>& mat 
)
:
    fdMatrix<T>::fdMatrix( mat.transportedFld_ )
{
    (*this) = mat;
}

template <class T>
fdMatrix<T>::fdMatrix
(
    const std::shared_ptr<fdMatrix<T> >& m
)
:
    fdMatrix<T>::fdMatrix( (*m) )
{
}
 


//====

template <class T>
void fdMatrix<T>::solve()
{
    int dim = (this->mesh().ni()-settings::m())*(this->mesh().nj()-settings::m())*(this->mesh().nk()-settings::m());

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x0(dim,pTraits<T>::nComp);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b(dim,pTraits<T>::nComp);
    std::vector<Eigen::Triplet<double> > coefs;
 
    coefs.reserve(7*dim);

    for( int k=settings::m()/2; k<this->mesh().nk()-settings::m()/2; k++ )
    {
        for( int j=settings::m()/2; j<this->mesh().nj()-settings::m()/2; j++ )
        {
            for( int i=settings::m()/2; i<this->mesh().ni()-settings::m()/2; i++ )
            {
                int p = (k-settings::m()/2)*(this->mesh().nj()-settings::m())*(this->mesh().ni()-settings::m())+(j-settings::m()/2)*(this->mesh().ni()-settings::m())+(i-settings::m()/2);

                coefs.push_back( Eigen::Triplet<double>( p, p, this->ap()()[i][j][k] ) );
                
                for( int c=0; c< pTraits<T>::nComp; c++ )
                {
                    b(p,c) = component(this->operator()()[i][j][k], c);
                    x0(p,c) = component(transportedFld_()[i][j][k], c);
                }


                int eas = p + 1;
                int wes = p - 1;
                int top = p + (this->mesh().nj()-settings::m())*(this->mesh().ni()-settings::m());
                int bot = p - (this->mesh().nj()-settings::m())*(this->mesh().ni()-settings::m());
                int nor = p + (this->mesh().ni()-settings::m());
                int sou = p - (this->mesh().ni()-settings::m());


                if( i>settings::m()/2 )
                {
                    coefs.push_back( Eigen::Triplet<double>( p, wes, this->al()[3]()[i][j][k] ) );
                }
                else
                {
                    coefs.push_back( Eigen::Triplet<double>( p, p + (this->ni()-settings::m())-1, this->al()[3]()[i][j][k] ) );
//                    for( int c=0; c< pTraits<T>::nComp; c++ )
//                    {
//                        b(p,c) -= this->al()[3]()[i][j][k] * component(transportedFld_()[i-1][j][k], c);
//                    }
                }


                if( i<this->mesh().ni()-settings::m()/2-1 )
                {
                    coefs.push_back( Eigen::Triplet<double>( p, eas, this->al()[2]()[i][j][k] ) );
                }
                else
                {
                    coefs.push_back( Eigen::Triplet<double>( p, p - (this->ni()-settings::m())+1, this->al()[2]()[i][j][k] ) );
//                    for( int c=0; c< pTraits<T>::nComp; c++ )
//                    {
//                        b(p,c) -= this->al()[2]()[i][j][k] * component(transportedFld_()[i+1][j][k], c);
//                    }
                }


                if( j>settings::m()/2 )
                {
                    coefs.push_back( Eigen::Triplet<double>( p, sou, this->al()[1]()[i][j][k] ) );
                }
                else
                {
                    coefs.push_back( Eigen::Triplet<double>( p, p+(this->mesh().ni()-settings::m())*(this->nj()-settings::m()-1), this->al()[1]()[i][j][k] ) );

//                    for( int c=0; c< pTraits<T>::nComp; c++ )
//                    {
//                        b(p,c) -= this->al()[1]()[i][j][k] * component(transportedFld_()[i][j-1][k], c);
//                    }
                }


                if( j<this->mesh().nj()-settings::m()/2-1 )
                {
                    coefs.push_back( Eigen::Triplet<double>( p, nor, this->al()[0]()[i][j][k] ) );
                }
                else
                {
                    coefs.push_back( Eigen::Triplet<double>( p, p-(this->mesh().ni()-settings::m())*(this->nj()-settings::m()-1), this->al()[0]()[i][j][k] ) );
//                    for( int c=0; c< pTraits<T>::nComp; c++ )
//                    {
//                        b(p,c) -= this->al()[0]()[i][j][k] * component(transportedFld_()[i][j+1][k], c);
//                    }
                }


                if( k>settings::m()/2 )
                {
                    coefs.push_back( Eigen::Triplet<double>( p, bot, this->al()[5]()[i][j][k] ) );
                }
                else
                {
                    coefs.push_back( Eigen::Triplet<double>( p, p+ (this->mesh().nj()-settings::m())*(this->mesh().ni()-settings::m())*(this->nk()-settings::m()-1), this->al()[5]()[i][j][k] ) );
//                    for( int c=0; c< pTraits<T>::nComp; c++ )
//                    {
//                        b(p,c) -= this->al()[5]()[i][j][k] * component(transportedFld_()[i][j][k-1], c);
//                    }
                }


                if( k<this->mesh().nk()-settings::m()/2-1 )
                {
                    coefs.push_back( Eigen::Triplet<double>( p, top, this->al()[4]()[i][j][k] ) );
                }
                else
                {
                    coefs.push_back( Eigen::Triplet<double>( p, p-(this->mesh().nj()-settings::m())*(this->mesh().ni()-settings::m())*(this->nk()-settings::m()-1), this->al()[4]()[i][j][k] ) );
//                    for( int c=0; c< pTraits<T>::nComp; c++ )
//                    {
//                        b(p,c) -= this->al()[4]()[i][j][k] * component(transportedFld_()[i][j][k+1], c);
//                    }
                }

            }
        }
    }
    
    Eigen::SparseMatrix<double> A(dim,dim);
    A.reserve(coefs.size());
    A.setFromTriplets(coefs.begin(), coefs.end());

//    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver(A);
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver(A);
//    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::IncompleteCholesky<double, Eigen::Lower> > solver(A);
//    Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver(A);

    solver.setTolerance(1e-8);
    solver.setMaxIterations(300);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x(dim,pTraits<T>::nComp);
    x = solver.solveWithGuess(b, x0);



    std::cout << "#iterations:     " << solver.iterations() << ". Error: "<< solver.error()<< std::endl;

    for( int i=settings::m()/2; i<this->mesh().ni()-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<this->mesh().nj()-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<this->mesh().nk()-settings::m()/2; k++ )
            { 
                int p = (k-settings::m()/2)*(this->mesh().nj()-settings::m())*(this->mesh().ni()-settings::m())+(j-settings::m()/2)*(this->mesh().ni()-settings::m())+(i-settings::m()/2);
                for( int c=0; c< pTraits<T>::nComp; c++ )
                {
                    component( transportedFld_()[i][j][k], c ) = x(p,c);
                }
            }
        }
    }

}


template <class T>
void fdMatrix<T>::deferredCorrection
(
    const std::shared_ptr<Field<T> >& hoa
)
{
    for( int i=settings::m(); i<this->ni()-settings::m(); i++ )
    {
        for( int j=settings::m(); j<this->nj()-settings::m(); j++ )
        {
            for( int k=settings::m(); k<this->nk()-settings::m(); k++ )
            {
                this->operator()()[i][j][k] = //rhs
               -hoa->operator()()[i][j][k]
               +this->al()[0]()[i][j][k] * transportedFld_()[i][j+1][k]
               +this->al()[1]()[i][j][k] * transportedFld_()[i][j-1][k]
               +this->al()[2]()[i][j][k] * transportedFld_()[i+1][j][k]
               +this->al()[3]()[i][j][k] * transportedFld_()[i-1][j][k]
               +this->al()[4]()[i][j][k] * transportedFld_()[i][j][k+1]
               +this->al()[5]()[i][j][k] * transportedFld_()[i][j][k-1]
               +this->ap()()[i][j][k] * transportedFld_()[i][j][k];
            }
        }
    } 
}

template <class T>
void fdMatrix<T>::relax( const double& urf )
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
    double val
)
{
    this->operator()()[i][j][k] += 1e15 * val;
    this->ap()()[i][j][k] *= 1e15;
}
    

template <class T>
T fdMatrix<T>::residual()
{
    T res(0.0);

    for( int i=3; i<this->ni()-3; i++ )
    {   
        for( int j=3; j<this->nj()-3; j++ )
        {   
            for( int k=3; k<this->nk()-3; k++ )
            {   
                for( int c=0; c<pTraits<T>::nComp; c++ )
                {
                    component(res, c) += 
                    pow
                    (
                        component
                        (
                            this->operator()()[i][j][k]
                           -this->al()[0]()[i][j][k] * transportedFld_()[i][j+1][k]
                           -this->al()[1]()[i][j][k] * transportedFld_()[i][j-1][k]
                           -this->al()[2]()[i][j][k] * transportedFld_()[i+1][j][k]
                           -this->al()[3]()[i][j][k] * transportedFld_()[i-1][j][k]
                           -this->al()[4]()[i][j][k] * transportedFld_()[i][j][k+1]
                           -this->al()[5]()[i][j][k] * transportedFld_()[i][j][k-1]
                           -this->ap()()[i][j][k] * transportedFld_()[i][j][k],
                            c
                        ),
                        2
                    );
                }
            }
        }
    }

    for( int c=0; c< pTraits<T>::nComp; c++ )
    {
        component(res,c) = pow( component(res,c)/( (this->ni()-settings::m())*(this->nj()-settings::m())*(this->nk()-settings::m()) ), 0.5 );
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
        this->al()[l] = r.al()[l];
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
        this->al()[l] += r.al()[l];
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
        this->al()[l] -= r.al()[l];
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
        this->al()[l] *= r.al()[l];
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
        ptr->al()[l] *= -1;
    }

    ptr->ap() *= -1;

    (*dynamic_cast<Field<T>* >(ptr.get())) *= -1;

    return ptr;
}

