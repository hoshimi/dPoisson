#!/usr/local/bin/rdmd
import std.stdio;
import myfft;
import myUtils;
import std.math;
immutable uint N = 32;
immutable real dx = 0.01;
immutable real dy = 0.01;
immutable real var_epsilon = 8.85418782 * 10.0^^(-12.0); // (s^4 A^2) / (m^3 kg)
immutable real Rqe = 1.60217657 * 10.0^^(-19.0); // C

real[] solvePoissonNeumann1D(real[] arg_rho, uint n, real dx){
    real[] rho = arg_rho.dup;

    /* for neumann boundary condition */
    rho[0] *= 2.0; rho[n-1] *= 2.0;

    real[] fourier_phi = cosfft(rho, n);

    for(int i = 0; i < n; i++) {
        int ng = i - n/2;
        real k = PI * (i - 1.0) / n;
        real b = dx * (-2.0 * cos(k) + 3.5);

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
        real b = dx * (-2.0 * cos(k) + 3.5);

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

real[N][N] solvePoissonNeumann2D(real[N][N] arg_rho, uint n, real dx, real dy, bool is_dirichlet = true){
    real[N][N] rho = arg_rho;
    real[] temp = new real[n];

    /* for neumann boundary condition */
    rho[0][] *= 2.0; rho[n-1][] *= 2.0;
    for(uint i = 0; i<n; i++) {
        rho[i][0] *= 2.0;
        rho[i][n-1] *= 2.0;
    }

    /* -- x-axis sine transform -- */
    for(uint i = 0; i<n; i++) {
        temp[] = rho[i].dup;
        rho[i] = cosfft(temp, n);
    }

    /* -- y-axis sine transform -- */
    for(uint i = 0; i<n; i++) {
        for(uint j = 0; j<n; j++) {
            temp[j] = rho[j][i];
        }
        
        temp = cosfft(temp, n);

        for(uint j = 0; j<n; j++) {
            rho[j][i] = temp[j];
        }
    }

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            real kx = PI * (i - 1.0) / n;
            real ky = PI * (j - 1.0) / n;
            /* Why this minus 2.4 is need??????? */
            real b = dx * (-2.0 * (cos(kx) + cos(ky) - 2.34));

            if( b != 0.0 ){
                rho[i][j] /= ( var_epsilon * b^^2.0 );
            }
        }
    }

    /* -- y-axis inverse sine transform -- */
    for(uint i = 0; i<n; i++) {
        for(uint j = 0; j<n; j++) {
            temp[j] = rho[j][i];
        }
        
        temp = icosfft(temp, n);

        for(uint j = 0; j<n; j++) {
            rho[j][i] = temp[j];
        }
    }

    /* -- x-axis inverse sine transform -- */
    for(uint i = 0; i<n; i++) {
        temp[] = rho[i].dup;
        rho[i] = icosfft(temp, n);
    }

    return rho;
}

real[N][N] solvePoissonDirichlet2D(real[N][N] arg_rho, uint n, real dx, real dy, bool is_dirichlet = true){
    real[N][N] rho = arg_rho;
    real[] temp = new real[n];

    /* -- x-axis sine transform -- */
    for(uint i = 0; i<n; i++) {
        temp[] = rho[i].dup;
        rho[i] = simple_sint(temp, n);
    }

    /* -- y-axis sine transform -- */
    for(uint i = 0; i<n; i++) {
        for(uint j = 0; j<n; j++) {
            temp[j] = rho[j][i];
        }
        
        temp = simple_sint(temp, n);

        for(uint j = 0; j<n; j++) {
            rho[j][i] = temp[j];
        }
    }

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            real kx = PI * (i - 1.0) / n;
            real ky = PI * (j - 1.0) / n;
            /* Why this minus 2.4 is need??????? */
            real b = dx * (-2.0 * (cos(kx) + cos(ky) - 2.4));

            if( b != 0.0 ){
                rho[i][j] /= ( var_epsilon * b^^2.0 );
            }
        }
    }

    /* -- y-axis inverse sine transform -- */
    for(uint i = 0; i<n; i++) {
        for(uint j = 0; j<n; j++) {
            temp[j] = rho[j][i];
        }
        
        temp = simple_isint(temp, n);

        for(uint j = 0; j<n; j++) {
            rho[j][i] = temp[j];
        }
    }

    /* -- x-axis inverse sine transform -- */
    for(uint i = 0; i<n; i++) {
        temp[] = rho[i].dup;
        rho[i] = simple_isint(temp, n);
    }

    return rho;
}

real[N][N] solvePoisson2D(real[N][N] rho, uint n, real dx, real dy, bool is_dirichlet = true){
    if(is_dirichlet) {
        return solvePoissonDirichlet2D(rho, n, dx, dy);
    } else {
        return solvePoissonNeumann2D(rho, n, dx, dy);
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
    real[N][N] rho = 0.0;
    real[] rho_1d = new real[N];
    for(uint i; i < N; i++) {
        for(uint j; j < N; j++) {
            if(i == N/2 && j == N/2){
                rho[i][j] = 1.0 * Rqe;
            /*} else if(i == 3 * N/4 && j == N/2){
                rho[i][j] = -1.0 * Rqe;
            */} else {
                rho[i][j] = 0.0;
            }
        }
    }
    for(uint i; i < N; i++) {
        if(i == N/2){
            rho_1d[i] = 1.0 * Rqe;
        } else {
            rho_1d[i] = 0.0;
        }
    }

    real[N][N] ideal_phi = 0.0;
    real[] ideal_phi_1d = new real[N];
    for(int i = 0; i < N ; i++) {
        for(int j = 0; j < N ; j++) {
            real rx = dx * cast(real)abs(i - cast(int)N / 2);
            real ry = dy * cast(real)abs(j - cast(int)N / 2);
            ideal_phi[i][j] = Rqe / ( var_epsilon * 4.0 * PI * (rx^^2.0 + ry^^2.0));
            ideal_phi_1d[i] = Rqe / ( var_epsilon * 4.0 * PI * rx^^2.0);
        }
    }

    /*
    for(int i = 0; i < N ; i++) {
        for(int j = 0; j < N ; j++) {
            real rx = dx * cast(real)abs(i - 3.0 * cast(int)N / 4);
            real ry = dy * cast(real)abs(j - cast(int)N / 2);
            ideal_phi[i][j] += -1.0 * Rqe / ( var_epsilon * 4.0 * PI * (rx^^2.0 + ry^^2.0));
        }
    }*/

    // real[N][N] phi = solvePoisson2D(rho, N, dx, dy, false);
    // outputCSV!(real[N][N])(phi, "test-phi.csv", " ", true);
    // outputCSV!(real[N][N])(ideal_phi, "test-phi.csv", " ", false);
    // real[] phi_1d = solvePoisson1D(rho_1d, N, dx, false);
    // outputCSV1d!(real[])(phi_1d, "test-phi-1d.csv", " ", true);
    // outputCSV1d!(real[])(ideal_phi_1d, "test-phi-1d.csv", " ", false);
}
