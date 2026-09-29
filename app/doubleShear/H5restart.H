{
//TODO restart function is unfinished and untested
//restart from HDF5 file at given arg -t
const H5std_string FILE_NAME_RESTART("./h5data/" + std::to_string(time.curTime()) + ".h5");
    if (parallelCom::master())
        std::cout << "Reading Restart data from file " << FILE_NAME_RESTART
                  << " at current time = " << time.curTime() << std::endl;

    // --- Open file with MPI I/O
    hid_t plist_fa = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_fa, parallelCom::world(), MPI_INFO_NULL);
    H5::H5File file(FILE_NAME_RESTART, H5F_ACC_RDONLY, H5P_DEFAULT, plist_fa);
    H5Pclose(plist_fa);

    // --- Create dataset transfer property list for collective reads
    hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);

    // --- Per-rank local block size
    hsize_t ini = mesh.glob_n()[0] / parallelCom::ni();
    hsize_t inj = mesh.glob_n()[1] / parallelCom::nj();
    hsize_t ink = mesh.glob_n()[2] / parallelCom::nk();

    // --- Compute hyperslab start/count/stride/block
    hsize_t start[3]  = { ini * parallelCom::i(),
                           inj * parallelCom::j(),
                           ink * parallelCom::k() };
    hsize_t count[3]  = { ini, inj, ink };
    hsize_t stride[3] = { 1, 1, 1 };
    hsize_t block[3]  = { 1, 1, 1 };


    // --- Create memory dataspace for this rank (shape = local count)
    hsize_t mem_dims[3] = { count[0], count[1], count[2] };
    hid_t memspace_id = H5Screate_simple(3, mem_dims, nullptr);

    // --- List of datasets and buffers 
    std::vector<std::string> dataset_names = { "Ux.h5", "Uy.h5", "Uz.h5",
                                               "Bx.h5", "By.h5", "Bz.h5",
                                               "p.h5",  "pB.h5" };
    std::vector<std::vector<double>> buffers(dataset_names.size());

    size_t nloc = static_cast<size_t>(count[0]) *
                  static_cast<size_t>(count[1]) *
                  static_cast<size_t>(count[2]);

    // Allocate all buffers
    for (auto& buf : buffers)
        buf.resize(nloc);


    // --- Loop over all datasets
    for (size_t idx_ds = 0; idx_ds < dataset_names.size(); ++idx_ds)
    {
        H5::DataSet dataset = file.openDataSet(dataset_names[idx_ds]);
        hid_t dataset_id = dataset.getId();

        // File dataspace
        hid_t filespace_id = H5Dget_space(dataset_id);
        H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, start, stride, count, block);

        // Read into local buffer
        herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE,
                                memspace_id, filespace_id,
                                xfer_plist, buffers[idx_ds].data());
        H5Sclose(filespace_id);

        if (status < 0)
        {
            std::cerr << "H5Dread failed for dataset " << dataset_names[idx_ds] << std::endl;
            MPI_Abort(parallelCom::world(), 1);
        }
    }

    // --- Map buffers into solver arrays
    size_t idx = 0;
    for (hsize_t i = 0; i < count[0]; ++i)
        for (hsize_t j = 0; j < count[1]; ++j)
            for (hsize_t k = 0; k < count[2]; ++k)
            {
                idx = i * count[1] * count[2] + j * count[2] + k;

                U((i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2)).x() = buffers[0][idx];
                U((i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2)).y() = buffers[1][idx];
                U((i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2)).z() = buffers[2][idx];

                B((i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2)).x() = buffers[3][idx];
                B((i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2)).y() = buffers[4][idx];
                B((i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2)).z() = buffers[5][idx];

                p((i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2))  = buffers[6][idx];
                pB((i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2)) = buffers[7][idx];
            }

    // --- Cleanup
    H5Sclose(memspace_id);
    H5Pclose(xfer_plist);
}

