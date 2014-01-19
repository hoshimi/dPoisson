#!/usr/local/bin/rdmd
import std.stdio;
import myfft;
import std.math;
immutable uint N = 128;
immutable real dx = 1.0;
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

real[] solvePoisson1D(real[] arg_rho, uint n, real dx){
    real[] rho = arg_rho.dup;
    
    /* for neumann boundary condition */
    rho[0] *= 2.0; rho[n-1] *= 2.0;

    real[] fourier_phi = cosfft(rho, n, dx);

    for(int i = 1; i < n; i++) {
        int ng = i - n/2;
        real k = PI * i / n;
        real kdx_per_2 = k * dx / 2.0;
        real sinckdx_per_2 = sinc(kdx_per_2);
        real K2 = k^^2 * (sinckdx_per_2)^^2;
        fourier_phi[i] /= (K2);
    }
    
    foreach(ref elem;fourier_phi){
        if(elem == real.infinity){
            elem = real.max;
        }
        writeln(elem);
    }

    return icosfft(fourier_phi, n, 1.0 / (n * dx));
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
        if(i == N/4){
            elem = 1.0;
        } else if(i == 3 * N/4){
            elem = -1.0;
        } else {
            elem = 0.0;
        }
    }

    real[] ideal_phi = new real[N];
    for(int i = 0; i < N ; i++) {
        ideal_phi[i] = 1.0 / (4.0 * PI * (abs(i - cast(int)N/4)^^2.0));
    }
    for(int i = 0; i < N ; i++) {
        ideal_phi[i] += -1.0 / (4.0 * PI * (abs(i - 3.0 * cast(int)N/4)^^2.0));
    }

    real[] phi = solvePoisson1D(rho, N, dx);

    outputCSV!(real)(phi, "test-phi.csv", " ", true);
    outputCSV!(real)(ideal_phi, "test-phi.csv", " ", false);
}
