module iteration;
import std.math;

T Jacobi3D(T)(T data, const int nx, const int ny, const int nz, double h, T rho){
    T t_data = data;
    for(int i = 1; i < nx-1; i++){
        for(int j = 1; j < ny-1; j++){
            for(int k = 1; k < nz-1; k++){
                t_data[i][j][k] = (data[i+1][j][k] + data[i-1][j][k] + data[i][j+1][k] + data[i][j-1][k] + data[i][j][k+1] + data[i][j][k-1])/6.0 + h * rho[i][j][k];
            }
        }
    }

    return t_data;
}

T GaussSeidel3D(T)(T data, const int nx, const int ny, const int nz, double h, T rho){
    T t_data = data;
    for(int i = 1; i < nx-1; i++){
        for(int j = 1; j < ny-1; j++){
            for(int k = 1; k < nz-1; k++){
                t_data[i][j][k] = (data[i+1][j][k] + t_data[i-1][j][k] + data[i][j+1][k] + t_data[i][j-1][k] + data[i][j][k+1] + t_data[i][j][k-1])/6.0 + h * rho[i][j][k];
            }
        }
    }

    return t_data;
}

T Sor3D(T)(T data, const int nx, const int ny, const int nz, double h, T rho, double omega){
    T t_data = data;
    for(int i = 1; i < nx-1; i++){
        for(int j = 1; j < ny-1; j++){
            for(int k = 1; k < nz-1; k++){
                t_data[i][j][k] = (1-omega) * data[i][j][k] + omega*(data[i+1][j][k] + t_data[i-1][j][k] + data[i][j+1][k] + t_data[i][j-1][k] + data[i][j][k+1] + t_data[i][j][k-1] + 6.0*h*rho[i][j][k])/6.0;
            }
        }
    }

    return t_data;
}

T ChebyshevSor3D(T)(T data, const int nx, const int ny, const int nz, double h, T rho, double omega){
    T t_data = data;
    for(int i = 1; i < nx-1; i++){
        for(int j = 1; j < ny-1; j++){
            for(int k = 1; k < nz-1; k++){
                t_data[i][j][k] = (1-omega) * data[i][j][k] + omega*(data[i+1][j][k] + t_data[i-1][j][k] + data[i][j+1][k] + t_data[i][j-1][k] + data[i][j][k+1] + t_data[i][j][k-1] + 6.0*h*rho[i][j][k])/6.0;
            }
        }
    }

    return t_data;
}

T calcResidual3D(T)(T data, const int nx, const int ny, const int nz, double h, T rho){
    T res = data;
    for(int i = 1; i < nx-1; i++){
        for(int j = 1; j < ny-1; j++){
            for(int k = 1; k < nz-1; k++){
                res[i][j][k] = (data[i+1][j][k] + data[i][j][k] + data[i][j+1][k] + data[i][j-1][k] + data[i][j][k+1] + data[i][j][k-1] - 6.0*data[i][j][k])/6.0 + rho[i][j][k];
            }
        }
    }

    return res;
}

double calcAverage3D(T)(T data, const int nx, const int ny, const int nz){
    double ave = 0.0; int cnt = 0;

    for(int i = 1; i < nx-1; i++){
        for(int j = 1; j < ny-1; j++){
            for(int k = 1; k < nz-1; k++){
                if(data[i][j][k] != 0.0){
                    ave += data[i][j][k];
                    cnt++;
                }
            }
        }
    }

    return ave/(cast(double)cnt);
}

double calcResidualAverage3D(T)(T data, const int nx, const int ny, const int nz, double h, T rho){
    return calcAverage3D!(T)(calcResidual3D!(T)(data, nx, ny, nz, h, rho), nx, ny, nz);
}
