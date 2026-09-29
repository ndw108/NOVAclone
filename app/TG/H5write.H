{
        std::string H5title = "./h5data/" + std::to_string( (time.curTime()-time.dt()) ) + ".h5";

        //Create properties list & MPI access
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, parallelCom::world(), MPI_INFO_NULL);

        //collective metadata operations
        H5Pset_coll_metadata_write(plist_id, 1);
        H5Pset_all_coll_metadata_ops(plist_id, 1);

	//data striping options, adjust as needed
        H5Pset_alignment(plist_id, 0, 4*1024*1024);  // 4 MiB

        //create a file
        hid_t file_id = H5Fcreate(H5title.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

        hsize_t hni = mesh.glob_n()[0]/parallelCom::ni();
        hsize_t hnj = mesh.glob_n()[1]/parallelCom::nj();
        hsize_t hnk = mesh.glob_n()[2]/parallelCom::nk();

        //create a data space (full domain)
        hsize_t dims[] = {   static_cast<hsize_t>((hni)*parallelCom::ni() ),
                             static_cast<hsize_t>((hnj)*parallelCom::nj() ),
                             static_cast<hsize_t>((hnk)*parallelCom::nk() )  };

        //create the local dataspace for each MPI rank
        hsize_t local_dims[3] = {hni, hnj, hnk};
        hsize_t start[3] = {static_cast<hsize_t>((hni*parallelCom::i())) , static_cast<hsize_t>((hnj*parallelCom::j())) , static_cast<hsize_t>((hnk*parallelCom::k())) };
        hsize_t count[3] = {hni, hnj, hnk};

        //Define dataset creation property list
        hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_NEVER);

        //MPI options
	hid_t xf_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(xf_id, H5FD_MPIO_COLLECTIVE);

        //create a dataspace
        hid_t dataspace_id = H5Screate_simple(3, dims, NULL);

        //size of memory
        hsize_t mem_dims[] = { hnk, hnj, hni};
        hid_t memspace_id = H5Screate_simple(3, mem_dims, NULL);

        //select area to write to
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL);

        //1d vector buffer
        std::unique_ptr<std::vector<double>> up_vectorField;
        std::vector<double>* vectorField;

        up_vectorField = std::make_unique<std::vector<double>>(hni * hnj * hnk);
        vectorField = up_vectorField.get();

	//Add or remove fields as desired
        //Ux
        int n=0;
        for (int i = 0; i < hni; ++i) {
            for (int j = 0; j < hnj; ++j) {
                for (int k = 0; k < hnk; ++k) {
                vectorField[0][n] = U( (i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2) ).x();
                n++;
                }
            }
        }
        std::string dataSetTitle = "Ux.h5";

        hid_t dataset_id = H5Dcreate2(file_id, dataSetTitle.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, xf_id, vectorField->data());
        H5Dclose(dataset_id);
        vectorField->clear();


        //Uy
        n=0;
        for (int i = 0; i < hni; ++i) {
            for (int j = 0; j < hnj; ++j) {
                for (int k = 0; k < hnk; ++k) {
                vectorField[0][n] = U( (i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2) ).y();
                n++;
                }
           }
        }
        dataSetTitle = "Uy.h5";
        dataset_id = H5Dcreate2(file_id, dataSetTitle.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, xf_id, vectorField->data());
        H5Dclose(dataset_id);
        vectorField->clear();


        //Uz
        n=0;
        for (int i = 0; i < hni; ++i) {
            for (int j = 0; j < hnj; ++j) {
                for (int k = 0; k < hnk; ++k) {
                vectorField[0][n] = U( (i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2) ).z();
                n++;
                }
           }
        }
        dataSetTitle = "Uz.h5";
        dataset_id = H5Dcreate2(file_id, dataSetTitle.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, xf_id, vectorField->data());
        H5Dclose(dataset_id);
        vectorField->clear();


        //p
        n=0;
        for (int i = 0; i < hni; ++i) {
            for (int j = 0; j < hnj; ++j) {
                for (int k = 0; k < hnk; ++k) {
                vectorField[0][n] = p( (i+settings::m()/2), (j+settings::m()/2), (k+settings::m()/2) );
                n++;
                }
           }
        }
        dataSetTitle = "p.h5";
        dataset_id = H5Dcreate2(file_id, dataSetTitle.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, xf_id, vectorField->data());
        H5Dclose(dataset_id);
        vectorField->clear();


	//cleanup
        H5Sclose(dataspace_id);
        H5Sclose(memspace_id);
        H5Pclose(dcpl_id);
        H5Pclose(xf_id);
        H5Pclose(plist_id);
        H5Fclose(file_id);

        if( parallelCom::master() )
          {
          std::cout << "data for t = "<< time.curTime()-time.dt() <<" written to ./h5data" << std::endl;
          }

}
