#!/usr/local/bin/rdmd
import std.stdio;
import myfft;
import std.math;
immutable uint N = 128;
immutable real dx = 0.01;
immutable real var_epsilon = 8.85418782 * 10.0^^(-12.0); // (s^4 A^2) / (m^3 kg)
immutable real Rqe = 1.60217657 * 10.0^^(-19.0); // C

void outputCSV(T)(T[] value, string filename, string spliter = ",", bool newfile = true) {
    if(newfile) {
        auto output = File(filename, "w");
        foreach(i, elem; value) {
            output.writeln(i, spliter, elem);
        }
    } else {
        auto output = File(filename, "a");
        output.writeln("\n\n## hogehoge");
        foreach(i, elem; value) {
            output.writeln(i, spliter, elem);
        }
    }
}

real[] solvePoissonNeumann1D(real[] arg_rho, uint n, real dx){
    real[] rho = arg_rho.dup;
    
    /* for neumann boundary condition */
    rho[0] *= 2.0; rho[n-1] *= 2.0;

    real[] fourier_phi = cosfft(rho, n);

    for(int i = 0; i < n; i++) {
        int ng = i - n/2;
        real k = PI * (i - 1.0) / n;
        real b = dx * (-2.0 * cos(k) + PI);
        
        if( b != 0.0 ){
            fourier_phi[i] /= (var_epsilon * b^^2.0);
        }
    }

    return icosfft(fourier_phi, n);
}

real[] solvePoissonDirichlet1D(real[] arg_rho, uint n, real dx){
    real[] rho = arg_rho.dup;

    real[] fourier_phi = simple_sint(rho, n);

    for(int i = 0; i < n; i++) {
        int ng = i - n/2;
        real k = PI * (i - 1.0) / n;
        real b = dx * (-2.0 * (cos(k) - 2.0));
        
        if( b != 0.0 ){
            fourier_phi[i] /= ( var_epsilon * b^^2.0 );
        }
    }

    return simple_isint(fourier_phi, n);
}

real[] solvePoisson1D(real[] rho, uint n, real dx, bool is_dirichlet = true){
    if(is_dirichlet) {
        return solvePoissonDirichlet1D(rho, n, dx);
    } else {
        return solvePoissonNeumann1D(rho, n, dx);
    }
}

real sinc(real x){
    if(x == 0.0) {
        return 1.0;
    } else {
        return sin(x)/x;
    }
}

void main() {
    real[] rho = new real[N];
    foreach(uint i, ref elem; rho){
        if(i == N/2){
            elem = 1.0 * Rqe;
        /*} else if(i == 3 * N/4){
            elem = -1.0 * Rqe;
        */} else {
            elem = 0.0;
        }
    }

    real[] ideal_phi = new real[N];
    for(int i = 0; i < N ; i++) {
        real r = dx * cast(real)abs(i - cast(int)N / 2);
        ideal_phi[i] = Rqe / ( var_epsilon * 4.0 * PI * r^^2.0);
    }
    /*for(int i = 0; i < N ; i++) {
        real r = dx * cast(real)abs(i - 3.0 * cast(int)N / 4);
        ideal_phi[i] += - Rqe / (var_epsilon * 4.0 * PI * r^^2.0);
    }*/

    real[] phi = solvePoisson1D(rho, N, dx, true);
    outputCSV!(real)(phi, "test-phi.csv", " ", true);
    outputCSV!(real)(ideal_phi, "test-phi.csv", " ", false);
}
