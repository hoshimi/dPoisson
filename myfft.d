#!/usr/local/bin/rdmd
module myfft;
import std.math;

// Simple Cosine Transform
// Type:DCT-II 
// Calculate F_i = dx * \sum_{j=0}^N f_j cos(PI * i * (j + 1/2)) / N)
real[] simple_cost(real[] value, immutable uint size) {
    real[] fourier_value = new real[size];
    fourier_value[] = 0.0;
    
    foreach(i;0..size) {
        foreach(j;0..size) {
            fourier_value[i] += value[j] * cos(cast(real)PI * i * (j + 0.5) / size) / (sqrt(cast(real)size));
        }
    }

    return fourier_value;
}

// Simple Sine Transform
// Type:DST-II 
// Calculate F_i = dx * \sum_{j=0}^N f_j sin(PI * i * (j - 1/2)) / N)
real[] simple_sint(real[] value, immutable uint size) {
    real[] fourier_value = new real[size];
    fourier_value[] = 0.0;
    
    foreach(i;0..size) {
        foreach(j;0..size) {
            fourier_value[i] += value[j] * sin(cast(real)PI * (i + 1.0) * (j + 0.5) / size);
        }
    }
    fourier_value[] *= sqrt(1.0 / cast(real)size);

    return fourier_value;
}

// Simple Inverse Sine Transform
// Type:DST-III 
// Calculate F_i = dx * \sum_{j=0}^N f_j sin(PI * i * (j - 1/2)) / N)
real[] simple_isint(real[] value, immutable uint size) {
    real[] fourier_value = new real[size];
    fourier_value[] = 0.0;
    
    foreach(i;0..size) {
        foreach(j;0..size) {
            fourier_value[i] += value[j] * sin(cast(real)PI * (i + 0.5) * (j + 1.0) / size);
        }
    }
    fourier_value[] *= sqrt(2.0 / cast(real)size);

    return fourier_value;
}

// Bit reverse scrambler
real[] scramble(real[] arg_value, uint n){
    real[] value = arg_value.dup;
    uint[] ip = new uint[n];
    uint k = n;
    uint m = 1;
    while (2 * m < k) {
        k /= 2;
        for(uint j = 0; j < m ; j++){
            ip[m + j] = ip[j] + k;
        }
        m = m * 2;
    }

    if( m == k) {
        for(uint i = 1; i < m ; i++ ) {
            for(uint j = 0; j < i; j++ ){
                uint ji = j + ip[i];
                uint ij = i + ip[j];
                real temp = value[ji];
                value[ji] = value[ij];
                value[ij] = temp;
            }
        }
    } else {
        for(uint i = 1; i < m ; i++ ) {
            for(uint j = 0; j < i; j++ ){
                uint ji = j + ip[i];
                uint ij = i + ip[j];
                real temp = value[ji];
                value[ji] = value[ij];
                value[ij] = temp;

                temp = value[ji + m];
                value[ji + m] = value[ij + m];
                value[ij + m] = temp;
            }
        }
    }

    return value;
}

// Time 1/2 Shift Real Fast Fourier Transform
// Type:DCT-II + DST-II
real[] rfft(real[] arg_value, immutable uint n) {
    real[] value = arg_value.dup;
    uint m, mh, mq, j0, j1, jr, ji, kr, ki, irev;
    real wr, wi, xr, xi;
    value = scramble(value, n);
    real theta = - 2.0 * PI / cast(real)n;
    value[0] /= sqrt(2.0);

    // for mh is 1, 2, 4, .. , n/2
    // kr, ki is complex conjugate of jr and ji
    for(mh = 1; (m = mh << 1) <= n; mh = m) {
        // jr = 0, 2mh, 4mh, .. n-m
        // kr = mh, 3mh, 5mh, ..
        for(jr = 0; jr < n; jr += m) {
            kr = jr + mh;
            xr = value[jr] - value[kr];
            value[jr] += value[kr];
            value[kr] = xr;
        }

        if ( mh == 1 ) continue;

        // mq is 1, 2, 4, .. , n/4
        mq = mh >> 1;
        irev = mq;

        for(ji = mq; ji < n; ji += m) {
            wr = cos(0.5 * theta * irev);
            wi = sin(0.5 * theta * irev);
            for(uint k = n >> 2; k > (irev ^= k); k >>= 1){};
            // ki = 3mh/2, 7mh/2, ..
            ki = ji + mh;

            xr = value[ji];
            xi = value[ki];

            value[ji] = wr * xr + wi * xi;
            value[ki] = wr * xi - wi * xr;
        }

        if(mh == 2) continue;

        irev = mq;
        // i = 0, 2mh, 4mh, .. n
        for(uint i = 0; i < n; i += m) {
            wr = cos(theta * irev);
            wi = sin(theta * irev);
            for(uint k = n >> 2; k > (irev ^= k); k >>= 1){};
            // j = 1, 2, 3, .. , mh/2
            for(uint j = 1; j < mq; j++){
                jr = i + j;
                kr = jr + mh;
                ki = i + mh - j;
                ji = ki + mh;

                xr = value[jr] - value[kr];
                xi = value[ji] - value[ki];

                value[jr] += value[kr];
                value[ji] += value[ki];

                value[ki] = wr * xr + wi * xi;
                value[kr] = wr * xi - wi * xr;
            }
        }
    }

    value[] *= sqrt(2.0/n);
    return value;
}

// Fast Sine Transform
// Type:DST-II 
// Calculate F_i = sqrt(2/n) * \sum_{j=0}^N f_j sin(PI * i * (j + 1/2)) / N)
real[] sinfft(real[] arg_value, immutable uint n) {
    real[] value = arg_value.dup;
    uint m, mh, mq, j0, j1, jr, ji, kr, ki, irev;
    real wr, wi, xr, xi;
    value = scramble(value, n);
    real theta = - PI / cast(real)n;
    value[0] /= sqrt(2.0);

    // for mh is 2, 4, .. , n/2
    // kr, ki is complex conjugate of jr and ji
    for(mh = 2; (m = mh << 1) <= n; mh = m) {
        // mq is 1, 2, 4, .. , n/4
        mq = mh >> 1;
        irev = 0;

        for(jr = 0; jr < n; jr += m) {
            wr = cos(0.5 * theta * (irev + mq));
            wi = sin(0.5 * theta * (irev + mq));
            for(uint k = n >> 2; k > (irev ^= k); k >>= 1){};

            // ki = jr + 1, jr + 2, jr + 4, .. jr + n/4
            // kr = n - 1 - jr, n - 2 - jr, .. 3n/4 - jr
            ki = jr + mq;
            kr = n - ki;
            
            // ji = n - mh - jr
            ji = kr - mq;

            xr = value[jr] - value[kr];
            xi = -value[ji] + value[ki];

            value[jr] += value[kr];
            value[ji] += value[ki];

            value[ki] = wr * xr + wi * xi;
            value[kr] = wr * xi - wi * xr;
        }

        if(mh == 2) continue;

        irev = 0;
        for(uint i = 0; i < n; i += m) {
            wr = cos(theta * (irev + mq));
            wi = sin(theta * (irev + mq));
            for(uint k = n >> 2; k > (irev ^= k); k >>= 1){};
            for(uint j = 1; j < mq; j++){
                jr = i + j;
                ki = i + mh - j;

                ji = n - jr;
                kr = n - ki;

                xr = value[jr] - value[kr];
                xi = -value[ji] - value[ki];

                value[jr] += value[kr];
                value[ji] -= value[ki];

                value[ki] = wr * xr + wi * xi;
                value[kr] = wr * xi - wi * xr;
            }
        }
    }

    if(n > 1) {
        m = n >> 1;
        xr = value[0] - value[m];

        value[0] += value[m];
        value[m] = xr * cos(0.5 * theta * m);
    }

    value[] *= sqrt(1.0/n);
    return value;
}

// Fast Cosine Transform
// Type:DCT-II 
// Calculate F_i = sqrt(2/n) * \sum_{j=0}^N f_j cos(PI * i * (j + 1/2)) / N)
real[] cosfft(real[] arg_value, immutable uint n) {
    real[] value = arg_value.dup;
    uint m, mh, mq, j0, j1, jr, ji, kr, ki, irev;
    real wr, wi, xr, xi;
    value = scramble(value, n);
    real theta = - PI / cast(real)n;
    value[0] /= sqrt(2.0);

    // for mh is 2, 4, .. , n/2
    // kr, ki is complex conjugate of jr and ji
    for(mh = 2; (m = mh << 1) <= n; mh = m) {
        // mq is 1, 2, 4, .. , n/4
        mq = mh >> 1;
        irev = 0;

        for(jr = 0; jr < n; jr += m) {
            wr = cos(0.5 * theta * (irev + mq));
            wi = sin(0.5 * theta * (irev + mq));
            for(uint k = n >> 2; k > (irev ^= k); k >>= 1){};

            // ki = jr + 1, jr + 2, jr + 4, .. jr + n/4
            // kr = n - 1 - jr, n - 2 - jr, .. 3n/4 - jr
            ki = jr + mq;
            kr = n - ki;
            
            // ji = n - mh - jr
            ji = kr - mq;

            xr = value[jr] - value[kr];
            xi = -value[ji] + value[ki];

            value[jr] += value[kr];
            value[ji] += value[ki];

            value[ki] = wr * xr + wi * xi;
            value[kr] = wr * xi - wi * xr;
        }

        if(mh == 2) continue;

        irev = 0;
        for(uint i = 0; i < n; i += m) {
            wr = cos(theta * (irev + mq));
            wi = sin(theta * (irev + mq));
            for(uint k = n >> 2; k > (irev ^= k); k >>= 1){};
            for(uint j = 1; j < mq; j++){
                jr = i + j;
                ki = i + mh - j;

                ji = n - jr;
                kr = n - ki;

                xr = value[jr] - value[kr];
                xi = -value[ji] - value[ki];

                value[jr] += value[kr];
                value[ji] -= value[ki];

                value[ki] = wr * xr + wi * xi;
                value[kr] = wr * xi - wi * xr;
            }
        }
    }

    if(n > 1) {
        m = n >> 1;
        xr = value[0] - value[m];

        value[0] += value[m];
        value[m] = xr * cos(0.5 * theta * m);
    }

    value[] *= sqrt(1.0/n);
    return value;
}

// Fast Inverse Cosine Transform
// Type:DCT-III
// Calculate F_i = sqrt(2/n) * \sum_{j=0}^N f_j cos(PI * j * (i + 1/2)) / N)
real[] icosfft(real[] arg_value, immutable uint n) {
    real[] value = arg_value.dup;

    uint m, mh, mq, j0, j1, jr, ji, kr, ki, irev;
    real wr, wi, xr, xi;
    real theta = PI / cast(real)n;
    value[0] /= sqrt(2.0);

    if(n > 1) {
        m = n >> 1;
        xr = value[m] * cos(0.5 * theta * m);

        value[m] = value[0] - xr;
        value[0] += xr;
    }

    for(m = n; (mh = m >> 1) >= 2; m = mh) {
        mq = mh >> 1;
        irev = 0;
        for(jr = 0; jr < n; jr += m) {
            wr = cos(0.5 * theta * (irev + mq));
            wi = sin(0.5 * theta * (irev + mq));
            for(uint k = n >> 2; k > (irev ^= k); k >>= 1){};
            ki = jr + mq;
            kr = n - ki;

            ji = kr - mq;
            
            xr = wr * value[ki] + wi * value[kr];
            xi = wr * value[kr] - wi * value[ki];

            value[kr] = value[jr] - xr;
            value[ki] = value[ji] + xi;

            value[jr] += xr;
            value[ji] -= xi;
        }

        if(mh == 2) continue;

        irev = 0;
        for(uint i = 0; i < n; i += m) {
            wr = cos(theta * (irev + mq));
            wi = sin(theta * (irev + mq));
            for(uint k = n >> 2; k > (irev ^= k); k >>= 1){};

            for(uint j = 1; j < mq; j++) {
                jr = i + j;
                ki = i + mh - j;
                ji = n - jr;
                kr = n - ki;

                xr = wr * value[ki] + wi * value[kr];
                xi = wr * value[kr] - wi * value[ki];

                value[kr] = value[jr] - xr;
                value[ki] = -value[ji] - xi;

                value[jr] += xr;
                value[ji] -= xi;
            }
        }
    }

    value[] *= sqrt(2.0/n);
    return scramble(value, n);
}

real[] bg_cosfft(real[] value, immutable uint size) {
    uint m, mh, j0, j1, irev;
    value = scramble(value, size);

    /* --- butterflies ---*/
    real theta = - PI / (2.0 * cast(real)size);

    for(mh = 1; (m = mh << 1) <= size; mh = m) {
        irev = 0;

        for(uint i = 0; i < size ; i += m) {
            real c2 = 2.0 * cos(theta * (irev + mh));
            for(uint k = size>>1; k > (irev ^= k); k >>= 1){};

            for(uint j = 0; j < mh ; j++){
                j0 = i + j;
                j1 = size - mh - i + j;
                real temp = value[j0] - value[j1];
                value[j0] += value[j1];
                value[j1] = c2 * temp;
            }
        }
    }

    /* --- adds --- */
    for(m = size >> 1;(mh = m >> 1) >= 1; m = mh) {
        for(uint j = 0; j < mh; j++) {
            j0 = m + mh + j;
            value[j0] = - value[j0] - value[j0 - m];
        }
        
        for(uint i = (m << 1) + mh; i < size; i += (m << 1)) {
            for(uint j = 0; j < mh; j++) {
                j0 = i + j;
                j1 = j0 + m;
                value[j0] -= value[j0 - m];
                value[j1] = -value[j1] - value[j0];
            }
        }
    }

    /* --- scaling --- */
    for(uint j = 1; j < size; j++){
        value[j] *= 0.5;
    }

    return value;
}
