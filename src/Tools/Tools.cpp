#include "Tools/Tools.h"

namespace tools
{

const scalar pi = std::numbers::pi;
const scalar eps = 1.0e-16;

void CFL( const Field<vector>& U, const Mesh& mesh )
{
    scalar cfl_x = 0.0;
    scalar cfl_y = 0.0;
    scalar cfl_z = 0.0;

    for( int i=settings::m()/2; i<mesh.ni()-settings::m()/2; i++ )
    {
        for( int j=settings::m()/2; j<mesh.nj()-settings::m()/2; j++ )
        {
            for( int k=settings::m()/2; k<mesh.nk()-settings::m()/2; k++ )
            {
                cfl_x = std::max( cfl_x, U(i, j, k).x() * mesh.time().dt() / mesh.dx());
                cfl_y = std::max( cfl_x, U(i, j, k).y() * mesh.time().dt() / mesh.dy());
                cfl_z = std::max( cfl_x, U(i, j, k).z() * mesh.time().dt() / mesh.dz());

            }
        }
    }

    reduce( cfl_x, maxOp<scalar>() );
    reduce( cfl_y, maxOp<scalar>() );
    reduce( cfl_z, maxOp<scalar>() );

    if( parallelCom::master() )
    {
        std::cout << "================================" << std::endl;
        std::cout << "  CFL x: " << cfl_x << ".\n";
        std::cout << "  CFL y: " << cfl_y << ".\n";
        std::cout << "  CFL z: " << cfl_z << ".\n";
        std::cout << "================================" << std::endl;
    }
}

void tdma( std::vector<scalar>& a, std::vector<scalar>& b, std::vector<scalar>& c, std::vector<scalar>& d, std::vector<scalar>& x )
{
    std::vector<scalar> cp, dp;
    cp.resize( a.size() );
    dp.resize( a.size() );

    cp[0] = c[0]/b[0];
    dp[0] = d[0]/b[0];
    for( int i=1; i<a.size(); i++ )
    {
        cp[i] = c[i]/(b[i]-a[i]*cp[i-1]);
        dp[i] = (d[i] - a[i]*dp[i-1])/(b[i]-a[i]*cp[i-1]);
    }

    x[x.size()-1] = dp[x.size()-1];
    for( int i=x.size()-2; i>=0; i-- )
    {
        x[i] = dp[i] - cp[i]*x[i+1];
    }
}

}
