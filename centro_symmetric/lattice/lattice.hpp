//
// Created by polya on 11/30/24.
//

#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <cmath>
#include <iostream>
#include <memory>
#include <random>



const auto PI=M_PI;

class lattice_centro_symmetric
{
public:
    lattice_centro_symmetric() : e2(std::random_device{}()), uni_0_2pi(0.0, 2 * PI) {
        // Initialize default values for members
    }
    lattice_centro_symmetric(const double& JVal, const int& NVal,const double &aVal, const std::shared_ptr<double[]>& ptr_theta_mat
        ,  const std::shared_ptr<double[]>& ptr_dist_inv_5_mat, const std::shared_ptr<double[]>& dist_inv_3_mat):e2(std::random_device{}()),uni_0_2pi(0.0,2*PI)
    {
        this->J=JVal;
        this->N=NVal;
        this->a=aVal;

        this->theta_mat=ptr_theta_mat;
        this->dist_inv_5_mat=ptr_dist_inv_5_mat;
        this->dist_inv_3_mat=dist_inv_3_mat;

    this->lower_bound_theta=-0.5*PI;
        this->upper_bound_theta=2.5*PI;

        this->compute_dist_arrays();
        this->init_angles();

        // int m0=2,m1=1;
        // int n0=m0,n1=m1;
        // int ind=m0*N+m1;
        // std::cout<<"dist_inv_3_mat[ind]="<<dist_inv_3_mat[ind]<<std::endl;
        // std::cout<<"dist_inv_5_mat[ind]="<<dist_inv_5_mat[ind]<<std::endl;
        // std::cout<<"E(n0,n1,m0,m1)="<<E(n0,n1,m0,m1)<<std::endl;

        // print_array(theta_mat,N);


    }


public:
    ///
    /// @param n0
    /// @param n1
    /// @return interaction energy between the dipole on posistion [n0,n1] with other dipoles
    double Energy_to_update(const int& n0, const int &n1);
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
    void print_array(const std::shared_ptr<double[]>& array, int N) {
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
    std::shared_ptr<double[]> theta_mat;//array containing values of angles
    std::shared_ptr<double[]> dist_inv_5_mat;//-3*J*a^{-3}*r^{-5},  V part
    std::shared_ptr<double[]> dist_inv_3_mat;//J*a^{-3}*r^{-3},  U part

    double lower_bound_theta,upper_bound_theta;

    std::ranlux24_base e2;
    std::uniform_real_distribution<> uni_0_2pi;
};



#endif //LATTICE_HPP
