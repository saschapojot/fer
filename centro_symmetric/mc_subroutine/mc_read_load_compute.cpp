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





void mc_computation::execute_mc_one_sweep(lattice_centro_symmetric & lat,  double& UCurr,
         double &U_time, double& proposal_time,double &rand_time,double &acc_reject_time)
   {
       double UNext;
       U_time=0;
       proposal_time=0;
       rand_time=0;
       acc_reject_time=0;

       std::chrono::duration<double> total_U_Elapsed{0};
       std::chrono::duration<double> total_proposal_Elapsed{0};
       std::chrono::duration<double> total_rand_Elapsed{0};
       std::chrono::duration<double> total_acc_reject_Elapsed{0};
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

               }

           }//end for n1
       }//end for n0



   }