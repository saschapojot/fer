#include "lattice/lattice.hpp"
#include "mc_subroutine/mc_read_load_compute.hpp"


int main(int argc, char *argv[])
{
    double J=2.0;
    int N=5;
    double a=1.1;
    lattice_centro_symmetric lat(J,N,a);
    lat.print_array(lat.dist_inv_3_mat,N);
    std::cout<<"======================"<<std::endl;
    lat.print_array(lat.dist_inv_5_mat,N);



}