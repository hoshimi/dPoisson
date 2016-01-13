#!/usr/local/bin/rdmd
import std.file, std.stdio;
import std.math;
import std.complex;
import std.numeric : Fft;
import std.conv;
import iteration;
import myUtils;

void main() {
    // Define const
    const int nx = 32; const int ny = nx; const int nz = nx;
    const double dx = 0.7; const double dy = dx; const double dz = dx;
    const double dt = dx*dx/6.0;
    const double eps0 = 8.85418782 * pow(10.0, -12);
    const double omega = 2.0/(1.0 + PI/nx);
    const double e = 1.60217657 * pow(10.0, -19);

    // Define field value
    double[nx][ny][nz] rho = 0.0;
    double[nx][ny][nz] phi = 0.0;
    double phi_b = 0.0;

    // Initial Charge
    for(int j = 0; j < ny; j++){
        if(j % 2 == 0){
            rho[nx/2][j][nz/2] = 1.0 * e/eps0;
        }else{
            rho[nx/2][j][nz/2] = -1.0 * e/eps0;
        }
    }
    rho[nx/2][ny/2][nz/2] = 1.0 * e/eps0;

    // Add boundary processing
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            // for x
            phi[0][i][j] = phi_b; phi[nx-1][i][j] = phi_b;

            // for y
            phi[i][0][j] = phi_b; phi[i][ny-1][j] = phi_b;

            // for z
            phi[i][j][0] = phi_b; phi[i][j][nz-1] = phi_b;
        }
        //phi[0][0][i] = sin(2.0 * 2.0 * PI * cast(double)i * dt);
    }

    //const double df = (1.0/(dt * cast(double)nx));

    //auto fft = new Fft(nx);
    //Complex!(double)[nx][ny][nz] c_phi;

    //for(int i = 0; i < nx; i++){
    //    for(int j = 0; j < ny; j++){
    //        for(int k = 0; k < nz; k++){
    //            c_phi[i][j] = fft.fft(phi[i][j][]);
    //        }
    //    }
    //}

    //writeln();
    //for(int i = 0; i < nx/2; i++){
    //    writef("%f: %f\n", i*df, abs(c_phi[0][0][i]));
    //}


    double limit = eps0*eps0;
    writeln(limit);
    double diff = 10000.0;
    double pre = 0.0;
    int n = 0;

    while(n < 1000 || fabs(diff) > limit){
        n++;
        //phi = GaussSeidel3D!(double[nx][ny][nz])(phi, nx, ny, nz, dt, rho);
        phi = Sor3D!(double[nx][ny][nz])(phi, nx, ny, nz, dt, rho, omega);
        diff = calcResidualAverage3D!(double[nx][ny][nz])(phi, nx, ny, nz, dx, rho);
        writefln("%d: %3e", n, diff);

        if(n > 10000) break;
    }

    outputCSV3D!(double[nx][ny][nz])(phi, "phi-sor.csv", " ", true);
}
