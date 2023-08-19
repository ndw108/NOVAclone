{    
    auto m = (im::laplacian(p_ptr) == (ex::div(U) / time.dt() ) );

    m->relax(0.3);

    m->solve( 1e-10, 100 );

    U -= time.dt() * ex::grad(p);
}
