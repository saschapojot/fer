//
// Created by polya on 11/30/24.
//

#include "lattice.hpp"

/// compute dist_inv_3_mat and dist_inv_5_mat
void lattice_centro_symmetric::compute_dist_arrays()
{
    for(int m0=0;m0<N;m0++)
    {
        for(int m1=0;m1<N;m1++)
        {
            int diff_tmp=std::pow(m0,2)+std::pow(m1,2);
            double diff_tmp_double=static_cast<double>(diff_tmp);
            int ind_tmp=m0*N+m1;
            if(diff_tmp==0)
            {
                this->dist_inv_3_mat[ind_tmp]=0.0;
                this->dist_inv_5_mat[ind_tmp]=0.0;
            }else
            {
                //U part
                this->dist_inv_3_mat[ind_tmp]=J*std::pow(a,-3.0)* std::pow(diff_tmp_double,-1.5);
                //V part
                this->dist_inv_5_mat[ind_tmp]=-3.0*J*std::pow(a,-3.0)* std::pow(diff_tmp_double,-2.5);
            }//end if-else diff_tmp

        }//end for m1
    }//end for m0
}//end compute_dist_arrays()



/// initialize theta_mat
void lattice_centro_symmetric::init_angles()
{
    for(int m0=0;m0<N;m0++)
    {
        for(int m1=0;m1<N;m1++)
        {
            int ind_tmp=twoD_ind_2_flat_ind(m0,m1);
            double thetaTmp=uni_0_2pi(e2);
            this->theta_mat[ind_tmp]=thetaTmp;
        }//end for m1
    }//end for m0
}


///
/// @param theta angle
/// @return dipole
void lattice_centro_symmetric::p_val(const double& theta,double &p0,double &p1)
{
    p0=std::cos(theta);
    p1=std::sin(theta);
}


///
/// @param n0
/// @param n1
/// @param m0
/// @param m1
/// @return U_{n0n1m0m1}
double lattice_centro_symmetric::U(const int& n0, const int &n1,const int &m0,const int & m1)
{
double p_n0n1_0,p_n0n1_1;
    double p_m0m1_0, p_m0m1_1;

    int abs_diff0=std::abs(m0-n0);
    int abs_diff1=std::abs(m1-n1);
    int flat_ind=twoD_ind_2_flat_ind(abs_diff0,abs_diff1);
    double dist_part=this->dist_inv_3_mat[flat_ind];

    double theta_n0n1,theta_m0m1;

    int ind_n0n1_flat=this->twoD_ind_2_flat_ind(n0,n1);
    int ind_m0m1_flat=this->twoD_ind_2_flat_ind(m0,m1);

    theta_n0n1=this->theta_mat[ind_n0n1_flat];
    theta_m0m1=this->theta_mat[ind_m0m1_flat];

    p_val(theta_n0n1,p_n0n1_0,p_n0n1_1);
    p_val(theta_m0m1,p_m0m1_0,p_m0m1_1);

    double U_val_tmp=dist_part*(p_n0n1_0*p_m0m1_0+p_n0n1_1*p_m0m1_1);

    return U_val_tmp;


}

///
/// @param n0
/// @param n1
/// @return flat index
int lattice_centro_symmetric::twoD_ind_2_flat_ind(const int &n0,const int &n1)
{

    return n0*N+n1;
}

///
/// @param n0
/// @param n1
/// @param m0
/// @param m1
/// @return V_{n0n1m0m1}
double lattice_centro_symmetric::V(const int& n0, const int &n1,const int &m0,const int & m1)
{

    double p_n0n1_0,p_n0n1_1;
    double p_m0m1_0, p_m0m1_1;

    double diff0=static_cast<double> (m0-n0);
    double diff1=static_cast<double> (m1-n1);

    int abs_diff0=std::abs(m0-n0);
    int abs_diff1=std::abs(m1-n1);

    int flat_ind=twoD_ind_2_flat_ind(abs_diff0,abs_diff1);

    double dist_part=this->dist_inv_5_mat[flat_ind];
    double theta_n0n1,theta_m0m1;

    int ind_n0n1_flat=this->twoD_ind_2_flat_ind(n0,n1);
    int ind_m0m1_flat=this->twoD_ind_2_flat_ind(m0,m1);

    theta_n0n1=this->theta_mat[ind_n0n1_flat];
    theta_m0m1=this->theta_mat[ind_m0m1_flat];

    p_val(theta_n0n1,p_n0n1_0,p_n0n1_1);
    p_val(theta_m0m1,p_m0m1_0,p_m0m1_1);

    double y1=dist_part*std::pow(diff0,2.0);

    double y2=dist_part*diff0*diff1;

    double y3=y2;

    double y4=dist_part*std::pow(diff1,2.0);

    double V_val=y1*p_n0n1_0*p_m0m1_0
                +y2*p_n0n1_0*p_m0m1_1
                +y3*p_n0n1_1*p_m0m1_0
                +y4*p_n0n1_1*p_m0m1_1;


    return V_val;






}

///
/// @param n0
/// @param n1
/// @param m0
/// @param m1
/// @return E_{n0n1m0m1}
double lattice_centro_symmetric::E(const int& n0, const int &n1,const int &m0,const int & m1)
{

    return U(n0,n1,m0,m1)+V(n0,n1,m0,m1);
}