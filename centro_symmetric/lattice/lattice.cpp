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
            int ind_tmp=m0*N+m1;
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