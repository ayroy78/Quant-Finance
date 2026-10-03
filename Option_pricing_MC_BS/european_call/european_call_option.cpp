
// finite-difference scheme for black-scholes

#include <iostream>
#include <cmath>
#include <fstream>

//this is a function to create linearly separated points. 
void linspace(double start, double end, double numpoints, double* points) {
	//double * points  = new double [numpoints];
	double dx = (end - start) / (numpoints - 1);
	for (int i = 0; i < numpoints; i++) {
		points[i] = (start + i * dx);
	}
}


int main(){

    // risk free rate
    double r = 0.05 ;

    //risk
    double sig = 0.2 ;

    //strike price
    double K = 1. ;

    //time to maturity
    double T = 1. ;

    //width of transformed S coordinate
    double w = 1.2 ; 

    //time integration array
    int Ntau = 1001; 
    double* tau = new double[Ntau];
    double tau0 = 0;
    double tauf = T;
    linspace(tau0, tauf, Ntau, tau);
    double dtau = T/(Ntau - 1);

    //grid points array
    int Nx = 101;
    double* x = new double[Nx];
    double x0 = -w;
    double xf = w;
    linspace(x0, xf, Nx, x);
    double dx = 2*w/(Nx - 1);

    //call option data storage array
    double ** f = new double*[Ntau];
    for (int i = 0; i<Ntau; i++){
        f[i] = new double[Nx];
    }

    //Set initial data
    for (int j = 0; j<Nx; j++){
        if ( K*(exp(x[j]) - 1. ) > 0 ) {
            f[0][j] = K*(exp(x[j]) - 1. );
        }
        else {
            f[0][j] = 0.;
        }
    }

    //Time Evolve with FTCS using BCs
    for (int i = 0; i<Ntau - 1; i++){
        //computational domain
        for (int j = 1; j<Nx -1; j++){
            f[i+1][j] = f[i][j] + dtau*( 0.5*pow(sig, 2.)*(f[i][j+1] - 2*f[i][j] + f[i][j-1])/pow(dx,2.) + (r - 0.5*pow(sig, 2.))*(f[i][j+1] - f[i][j-1])/(2*dx) - r*f[i][j] );
        }
        //boundary conditions applied to grid edges
        f[i+1][0] = 0.;
        f[i+1][Nx-1] = K*(exp(1.2) - exp(-r*tau[i+1]));
    }

    //save call data to output file
    std::ofstream myFile("european_call_option_data.txt");
    myFile.precision(17);   

    for (int i=0; i<Ntau; i++){
        myFile << tau[i] << " ";
        for (int j=0; j<Nx; j++){
            myFile << f[i][j] << " ";
        }
        myFile << "\n";
    }
   
    myFile.close();

    //delete pointers
    delete [] tau;

    delete [] x;

    for(int i = 0; i < Ntau; i++) {
        delete[] f[i]; // Free each row
    }
    delete[] f; // Free the pointer array

    return 0;
}