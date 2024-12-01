//
// Created by polya on 11/30/24.
//

#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <cmath>
#include <iostream>

#include <random>



const auto PI=M_PI;

class lattice_centro_symmetric
{
public:
    lattice_centro_symmetric(const double& JVal, const int& NVal,const double &aVal):e2(std::random_device{}()),uni_0_2pi(0.0,2*PI)
    {
        this->J=JVal;
        this->N=NVal;
        this->a=aVal;

        this->theta_mat=new double[N*N];
        this->dist_inv_5_mat=new double[N*N];
        this->dist_inv_3_mat=new double[N*N];



        this->compute_dist_arrays();
        this->init_angles();


    }
    ~lattice_centro_symmetric()
    {
        delete [] theta_mat;
        delete [] dist_inv_5_mat;
        delete [] dist_inv_3_mat;

    }

public:
    ///
    /// @param n0
    /// @param n1
    /// @param m0
    /// @param m1
    /// @return E_{n0n1m0m1}
    double E(const int& n0, const int &n1,const int &m0,const int & m1);
    ///
    /// @param n0
    /// @param n1
    /// @param m0
    /// @param m1
    /// @return V_{n0n1m0m1}
    double V(const int& n0, const int &n1,const int &m0,const int & m1);
    ///
    /// @param n0
    /// @param n1
    /// @param m0
    /// @param m1
    /// @return U_{n0n1m0m1}
    double U(const int& n0, const int &n1,const int &m0,const int & m1);
    ///
    /// @param theta angle
    /// @return dipole
    void p_val(const double& theta,double &p0,double &p1);



    /// compute dist_inv_3_mat and dist_inv_5_mat
    void compute_dist_arrays();

    /// initialize theta_mat
    void init_angles();

    ///
    /// @param n0
    /// @param n1
    /// @return flat index
    int twoD_ind_2_flat_ind(const int &n0,const int &n1);
    void print_array(double* array, int N) {
        for (int m0 = 0; m0 < N; ++m0) {
            for (int m1 = 0; m1 < N; ++m1) {
                std::cout << array[m0 * N + m1] << " ";
            }
            std::cout << std::endl;
        }
    }

public:
    double J;// coupling constant
    int N;//lattice length
    double a;//length of a unit cell
    double * theta_mat;//array containing values of angles
    double * dist_inv_5_mat;//-3*J*a^{-3}*r^{-5},  V part
    double * dist_inv_3_mat;//J*a^{-3}*r^{-3},  U part


    std::ranlux24_base e2;
    std::uniform_real_distribution<> uni_0_2pi;
};



#endif //LATTICE_HPP
