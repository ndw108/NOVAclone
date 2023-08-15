#include "Tools/Tools.h"

namespace tools
{

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

}
