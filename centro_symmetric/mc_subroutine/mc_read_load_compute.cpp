//
// Created by polya on 11/30/24.
//
#include "mc_read_load_compute.hpp"

///
/// @param x
/// @param leftEnd
/// @param rightEnd
/// @param eps
/// @return return a value within distance eps from x, on the open interval (leftEnd, rightEnd)
   double mc_computation::generate_uni_open_interval(const double &x, const double &leftEnd, const double &rightEnd, const double &eps)
   {

       double xMinusEps=x-eps;
       double xPlusEps=x+eps;

       double unif_left_end=xMinusEps<leftEnd?leftEnd:xMinusEps;
       double unif_right_end=xPlusEps>rightEnd?rightEnd:xPlusEps;

       //    std::random_device rd;
       //    std::ranlux24_base e2(rd());

       double unif_left_end_double_on_the_right=std::nextafter(unif_left_end, std::numeric_limits<double>::infinity());



       std::uniform_real_distribution<> distUnif(unif_left_end_double_on_the_right,unif_right_end); //[unif_left_end_double_on_the_right, unif_right_end)

       double xNext=distUnif(e2);
       return xNext;

   }



///
/// @param x proposed value
/// @param y current value
/// @param a left end of interval
/// @param b right end of interval
/// @param epsilon half length
/// @return proposal probability S(x|y)
double mc_computation::S_uni(const double &x, const double &y,const double &a, const double &b, const double &epsilon){

       if (a<y and y<a+epsilon){
           return 1.0/(y-a+epsilon);
       } else if( a+epsilon<=y and y<b+epsilon){
           return 1.0/(2.0*epsilon);
       }else if(b-epsilon<=y and y<b){
           return 1.0/(b-y+epsilon);
       } else{

           std::cerr<<"value out of range."<<std::endl;
           std::exit(10);


       }


   }//end S_uni




///
/// @param lat
/// @param pos position (n0*N+n1) of theta, to be updated
/// @return proposed new value of theta
double mc_computation::propose_uni(const lattice_centro_symmetric& lat, const int &pos)
   {
       double theta_old=lat.theta_mat[pos];

       double theta_new=generate_uni_open_interval(theta_old,lat.lower_bound_theta,lat.upper_bound_theta,h);

       return theta_new;
   }//end propose_uni




///
/// @param UCurr
/// @param UNext
/// @param theta_Curr
/// @param theta_Next
/// @return acceptantce ratio
double mc_computation::acceptance_ratio_uni_theta_update(const double &interval_lower_bound, const double& interval_upper_bound,const double& UCurr, const double &UNext,const double&theta_Curr, const double &theta_Next)
   {
       double numerator = -this->beta*UNext;
       double denominator=-this->beta*UCurr;

       double R=std::exp(numerator - denominator);

       double S_Curr_Next=S_uni(theta_Curr,theta_Next,interval_lower_bound,interval_upper_bound,h);

       double S_Next_Curr=S_uni(theta_Next,theta_Curr,interval_lower_bound,interval_upper_bound,h);

       double ratio=S_Curr_Next/S_Next_Curr;
       if (std::fetestexcept(FE_DIVBYZERO)) {
           std::cout << "Division by zero exception caught." << std::endl;
           std::exit(15);
       }
       if (std::isnan(ratio)) {
           std::cout << "The result is NaN." << std::endl;
           std::exit(15);
       }

       R*=ratio;

       return std::min(1.0,R);



   }//end acceptantce_ratio_uni_theta_update





void mc_computation::execute_mc_one_sweep(lattice_centro_symmetric & lat,  double& UCurr)
   {
       double UNext;
//       U_time=0;
//       proposal_time=0;
//       rand_time=0;
//       acc_reject_time=0;

//       std::chrono::duration<double> total_U_Elapsed{0};
//       std::chrono::duration<double> total_proposal_Elapsed{0};
//       std::chrono::duration<double> total_rand_Elapsed{0};
//       std::chrono::duration<double> total_acc_reject_Elapsed{0};
       double theta_Curr,theta_Next;
       for(int n0=0;n0<N;n0++)
       {
           for(int n1=0;n1<N;n1++)
           {
               int theta_ind=lat.twoD_ind_2_flat_ind(n0,n1);
               //before update
               theta_Curr=lat.theta_mat[theta_ind];
               UCurr=lat.Energy_to_update(n0,n1);
               //propose new theta value
               theta_Next=this->propose_uni(lat,theta_ind);

               //modify lat
               lat.theta_mat[theta_ind]=theta_Next;
               UNext=lat.Energy_to_update(n0,n1);

               double r=acceptance_ratio_uni_theta_update(lat.lower_bound_theta,lat.upper_bound_theta,UCurr,UNext,theta_Curr,theta_Next);
               double u = distUnif01(e2);
               if(u<=r)
               {
                   UCurr=UNext;

               }//accept
               else{
                   lat.theta_mat[theta_ind]=theta_Curr;
               }//reject

           }//end for n1
       }//end for n0



   }//end execute_mc_one_sweep



void mc_computation::execute_mc( const int & flushNum){

    double UCurr=0;
    int flushThisFileStart=this->flushLastFile+1;
//    std::cout<<"entering mc"<<std::endl;
    for (int fls = 0; fls < flushNum; fls++){
        const auto tMCStart{std::chrono::steady_clock::now()};
        for (int swp = 0; swp < sweepToWrite*sweep_multiple; swp++){
            const auto t_swp_Start{std::chrono::steady_clock::now()};
            this->execute_mc_one_sweep(this->lat,UCurr);
            if(swp%sweep_multiple==0){
                int swp_out=static_cast<int>(swp/sweep_multiple);
                U_data_ptr[swp_out]=UCurr;
                std::memcpy(theta_data_ptr.get()+swp_out*N*N,lat.theta_mat.get(),N*N*sizeof (double));
                //write U
            }//end write
            const auto t_swp_End{std::chrono::steady_clock::now()};
            const std::chrono::duration<double> swp_elapsed_secondsAll{t_swp_End - t_swp_Start};
//            std::cout << "sweep " + std::to_string(swp)  + ": "
//                      << swp_elapsed_secondsAll.count()  << " s" << std::endl;
        }//end sweep for
        int flushEnd=flushThisFileStart+fls;
        std::string fileNameMiddle =  "flushEnd" + std::to_string(flushEnd);
        std::string out_U_PickleFileName = out_U_path+"/" + fileNameMiddle + ".U.pkl";
        std::string out_theta_PickleFileName=out_theta_path+"/"+fileNameMiddle+".theta.pkl";

        save_array_to_pickle(U_data_ptr,sweepToWrite,out_U_PickleFileName);
        save_array_to_pickle(theta_data_ptr,sweepToWrite*N*N,out_theta_PickleFileName);
        const auto tMCEnd{std::chrono::steady_clock::now()};
        const std::chrono::duration<double> elapsed_secondsAll{tMCEnd - tMCStart};
        std::cout << "flush " + std::to_string(flushEnd)  + ": "
                  << elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;
    }//end flush for
    std::cout << "mc executed for " << flushNum << " flushes." << std::endl;
    std::cout<<"U_data_ptr[8]="<<U_data_ptr[8]<<std::endl;
}//end execute_mc


void mc_computation::save_array_to_pickle(const std::shared_ptr<double[]> &ptr, const int& size, const std::string& filename) {
    using namespace boost::python;
    namespace np = boost::python::numpy;
    std::cout << "Entering save, fileName=" << filename << ", size=" << size << std::endl;

    // Initialize Python interpreter and NumPy
    if (!Py_IsInitialized()) {
        Py_Initialize();
        if (!Py_IsInitialized()) {
            throw std::runtime_error("Failed to initialize Python interpreter");
        }
        np::initialize();  // Initialize NumPy after Python is initialized
    }

    try {
        std::cout << "Entering serialization" << std::endl;

        // Import the pickle module
        object pickle = import("pickle");
        object pickle_dumps = pickle.attr("dumps");

        // Convert C++ array to NumPy array using shared_ptr
        np::ndarray numpy_array = np::from_data(
                ptr.get(),                               // Use shared_ptr's raw pointer
                np::dtype::get_builtin<double>(),        // NumPy data type (double)
                boost::python::make_tuple(size),         // Shape of the array (1D array)
                boost::python::make_tuple(sizeof(double)),  // Strides
                object()                                 // Optional base object
        );

        // Serialize the NumPy array using pickle.dumps
        object serialized_array = pickle_dumps(numpy_array);

        // Access the raw bytes from the Python bytes object
        PyObject* py_bytes = serialized_array.ptr();
        if (PyBytes_Check(py_bytes)) {
            char* buffer;
            Py_ssize_t length;

            if (PyBytes_AsStringAndSize(py_bytes, &buffer, &length) == -1) {
                throw std::runtime_error("Failed to extract bytes data from serialized array");
            }

            // Write the raw bytes to a file
            std::ofstream file(filename, std::ios::binary);
            if (!file) {
                throw std::runtime_error("Failed to open file for writing");
            }
            file.write(buffer, length);
            file.close();

            std::cout << "Array serialized and written to file successfully." << std::endl;
        } else {
            throw std::runtime_error("Serialized data is not a bytes object");
        }

    } catch (const error_already_set&) {
        PyErr_Print();
        std::cerr << "Boost.Python error occurred." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
}


void mc_computation::load_pickle_data(const std::string& filename, std::shared_ptr<double[]>& data_ptr, std::size_t size)
{

    // Initialize Python and NumPy
    Py_Initialize();
    np::initialize();

    try {
        // Use Python's 'io' module to open the file directly in binary mode
        py::object io_module = py::import("io");
        py::object file = io_module.attr("open")(filename, "rb");  // Open file in binary mode

        // Import the 'pickle' module
        py::object pickle_module = py::import("pickle");

        // Use pickle.load to deserialize from the Python file object
        py::object loaded_data = pickle_module.attr("load")(file);

        // Close the file
        file.attr("close")();

        // Check if the loaded object is a NumPy array
        if (py::extract<np::ndarray>(loaded_data).check()) {
            np::ndarray np_array = py::extract<np::ndarray>(loaded_data);

            // Convert the NumPy array to a Python list using tolist()
            py::object py_list = np_array.attr("tolist")();

            // Ensure the list size matches the expected size
            ssize_t list_size = py::len(py_list);
            if (static_cast<std::size_t>(list_size) > size) {
                throw std::runtime_error("The provided shared_ptr array size is smaller than the list size.");
            }

            // Copy the data from the Python list to the shared_ptr array
            for (ssize_t i = 0; i < list_size; ++i) {
                data_ptr[i] = py::extract<double>(py_list[i]);
            }
        } else {
            throw std::runtime_error("Loaded data is not a NumPy array.");
        }
    }
    catch (py::error_already_set&) {
        PyErr_Print();
        throw std::runtime_error("Python error occurred.");
    }

}


void mc_computation::init_ptr_theta_mat(){

    std::string name;
    std::string theta_inFileName;
    if(this->flushLastFile==-1){
        name="init";
        theta_inFileName=out_theta_path+"/theta_"+name+".pkl";
        this->load_pickle_data(theta_inFileName,ptr_theta_mat,N*N);
    }//end -1 case
    else{
        name="flushEnd"+std::to_string(this->flushLastFile);
        theta_inFileName=out_theta_path+"/"+name+".theta.pkl";
        this->load_pickle_data(theta_inFileName,theta_data_ptr,sweepToWrite*N*N);
        std::memcpy(ptr_theta_mat.get(),theta_data_ptr.get()+N*N*(sweepToWrite-1),N*N*sizeof(double));
    }
}


void mc_computation::init_and_run(){

this->execute_mc(newFlushNum);

}