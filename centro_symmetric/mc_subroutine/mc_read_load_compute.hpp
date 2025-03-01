//
// Created by polya on 11/30/24.
//

#ifndef MC_READ_LOAD_COMPUTE_HPP
#define MC_READ_LOAD_COMPUTE_HPP
#include "../lattice/lattice.hpp"
#include <boost/filesystem.hpp>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <cfenv>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#endif //MC_READ_LOAD_COMPUTE_HPP
namespace fs = boost::filesystem;
namespace py = boost::python;
namespace np = boost::python::numpy;

class mc_computation
{
    //set seed
private:
    int seed=11;
public:
    mc_computation(const std::string &cppInParamsFileName):e2(seed),distUnif01(0.0, 1.0)
    {

        std::ifstream file(cppInParamsFileName);
        if (!file.is_open()) {
            std::cerr << "Failed to open the file." << std::endl;
            std::exit(20);
        }

        std::string line;
        int paramCounter = 0;

        while (std::getline(file, line))
        {

            // Check if the line is empty
            if (line.empty()) {
                continue; // Skip empty lines
            }
            std::istringstream iss(line);
            //read T
            if (paramCounter == 0)
            {
                iss >> T;
                if (T <= 0) {
                    std::cerr << "T must be >0" << std::endl;
                    std::exit(1);
                }//end if

                this->beta = 1 / T;

                paramCounter++;
                continue;


            }//end reading T

            //read a
            if (paramCounter==1)
            {
                iss>>a;
                if (a<=0)
                {
                    std::cerr << "a must be >0" << std::endl;
                    std::exit(1);
                }// end if
                paramCounter++;
                continue;
            }// end reading a

            //read J
            if(paramCounter==2)
            {
                iss>>J;
                paramCounter++;
                continue;
            }//end reading J

            //read N
            if(paramCounter==3)
            {
                iss>>N;
                if(N<=0)
                {
                    {
                        std::cerr << "N must be >0" << std::endl;
                        std::exit(1);
                    }// end if
                }

                paramCounter++;
                continue;
            }//end reading N


            //read sweepToWrite
            if(paramCounter==4)
            {
                iss>>sweepToWrite;
                paramCounter++;
                continue;
            }//end reading sweepToWrite

            //read newFlushNum
            if(paramCounter==5)
            {
                iss>>newFlushNum;
                paramCounter++;
                continue;
            }//end reading newFlushNum

            //read flushLastFile
            if(paramCounter==6){
                //if flushLastFileStr is "-1"
                //flushLastFile+1 will be 0
                iss>>flushLastFile;
                paramCounter++;
                continue;
            }//end reading flushLastFile
            //read TDirRoot
            if (paramCounter==7){
                iss>>TDirRoot;
                paramCounter++;
                continue;
            }//end reading TDirRoot

            //read U_dipole_dataDir
            if(paramCounter==8){
                iss>>U_dipole_dataDir;
                paramCounter++;
                continue;
            }//end reading U_dipole_dataDir


            //read h
            if(paramCounter==9){
                iss>>h;
                paramCounter++;



                continue;
            }// end h

            //read sweep_multiple
            if(paramCounter==10){
                iss>>sweep_multiple;
                paramCounter++;


                continue;
            }//end sweep_multiple


        }//end while

        try {
            this->ptr_theta_mat= std::shared_ptr<double[]>(new double[N*N],
                                                        std::default_delete<double[]>());
            this->ptr_dist_inv_5_mat= std::shared_ptr<double[]>(new double[N*N],
                                                        std::default_delete<double[]>());
            this->ptr_dist_inv_3_mat= std::shared_ptr<double[]>(new double[N*N],
                                                        std::default_delete<double[]>());

            this->U_data_ptr=std::shared_ptr<double[]>(new double[sweepToWrite ],
                                                       std::default_delete<double[]>());

            this->theta_data_ptr=std::shared_ptr<double[]>(new double[sweepToWrite*N*N ],
                                                       std::default_delete<double[]>());

        }
        catch (const std::bad_alloc &e) {
            std::cerr << "Memory allocation error: " << e.what() << std::endl;
            std::exit(2);
        } catch (const std::exception &e) {
            std::cerr << "Exception: " << e.what() << std::endl;
            std::exit(2);
        }
        this->out_theta_path=U_dipole_dataDir+"/theta/";
        this->out_U_path=U_dipole_dataDir+"/U/";
        this->init_ptr_theta_mat();

        this->lat=lattice_centro_symmetric(J,N,a,ptr_theta_mat
         ,   ptr_dist_inv_5_mat, ptr_dist_inv_3_mat);
//        lat.print_array(lat.theta_mat,N);
        std::cout << "T=" << T << std::endl;
        std::cout<<"beta="<<beta<<std::endl;
        std::cout << "a=" << a << std::endl;
        std::cout << "N=" << N << std::endl;
        std::cout << "J=" << J << std::endl;

        std::cout<<"sweepToWrite="<<sweepToWrite<<std::endl;
        std::cout<<"newFlushNum="<<newFlushNum<<std::endl;
        std::cout<<"flushLastFile="<<flushLastFile<<std::endl;
        std::cout<<"TDirRoot="<<TDirRoot<<std::endl;
        std::cout<<"U_dipole_dataDir="<<U_dipole_dataDir<<std::endl;
        std::cout << "h=" << h << std::endl;
        std::cout << "sweep_multiple=" << sweep_multiple << std::endl;
    }//end constructor



public:
    void init_and_run();
    void execute_mc( const int & flushNum);


    void execute_mc_one_sweep(lattice_centro_symmetric & lat,  double& UCurr);
    ///
    /// @param UCurr
    /// @param UNext
    /// @param theta_Curr
    /// @param theta_Next
    /// @return acceptantce ratio
    double acceptance_ratio_uni_theta_update(const double &interval_lower_bound, const double& interval_upper_bound,const double& UCurr, const double &UNext,const double&theta_Curr, const double &theta_Next);
    ///
    /// @param lat
    /// @param pos position (n0*N+n1) of theta, to be updated
    /// @return proposed new value of theta
    double propose_uni(const lattice_centro_symmetric& lat, const int &pos);
    ///
    /// @param x proposed value
    /// @param y current value
    /// @param a left end of interval
    /// @param b right end of interval
    /// @param epsilon half length
    /// @return proposal probability S(x|y)
    double S_uni(const double &x, const double &y,const double &a, const double &b, const double &epsilon);
    ///
    /// @param x
    /// @param leftEnd
    /// @param rightEnd
    /// @param eps
    /// @return return a value within distance eps from x, on the open interval (leftEnd, rightEnd)
    double generate_uni_open_interval(const double &x, const double &leftEnd, const double &rightEnd, const double &eps);

    void save_array_to_pickle(const std::shared_ptr<double[]> &ptr,const int& size,const std::string& filename);
    void load_pickle_data(const std::string& filename, std::shared_ptr<double[]>& data_ptr, std::size_t size);
    void init_ptr_theta_mat();

public:
    double T;// temperature
    double beta;
    double a;
    double J;
    int N;
    double h;// step size
    int sweepToWrite;
    int newFlushNum;
    int flushLastFile;
    std::string TDirRoot;
    std::string U_dipole_dataDir;
    std::ranlux24_base e2;
    std::uniform_real_distribution<> distUnif01;
    int sweep_multiple;
    std::string out_U_path;
    std::string out_theta_path;
    lattice_centro_symmetric lat;

    std::shared_ptr<double[]> ptr_theta_mat;//array containing values of angles
    std::shared_ptr<double[]> ptr_dist_inv_5_mat;//-3*J*a^{-3}*r^{-5},  V part
    std::shared_ptr<double[]> ptr_dist_inv_3_mat;//J*a^{-3}*r^{-3},  U part

    //data
    std::shared_ptr<double[]> U_data_ptr;//all U data
    std::shared_ptr<double[]>theta_data_ptr;//all theta_mat data

};