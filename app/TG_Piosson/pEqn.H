{    
    auto m = (im::laplacian(p_ptr) == (ex::div(U) / time.dt() ) );

    m->relax(0.1);

    m->solve( 1e-15, 50 );

    U -= time.dt() * ex::grad(p);
    U.correctBoundaryConditions();
}
