#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>

// for testing.
// prints array of length n as x-by-y tabbed block.
// assumes x*y = n else accesses memory outside array.
void print_array(int a[], int n, int x, int y){
    for (int i = 0; i < x; i++){
        for (int j = 0; j < y; j++){
            printf("%i\t", a[i*y + j]);
        }
        printf("\n");
    }
}

// schoolbook multiplication, degree n-1 polys
void SM(int a[], int b[], int c[], int n){
    for (int i = 0; i < 2*n; i++){
        c[i] = 0;
    }
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            c[i + j] += a[i]*b[j];
        }
    }
}

// karatsuba multiplication of integer polynomials.
// assumes n = degree + 1 is a power of 2.
// switches to schoolbook at degree = base - 1,
// which should also be a power of 2.
void KM(int a[], int b[], int c[], int n, int base){
    for (int i = 0; i < 2*n; i++){
        c[i] = 0;
    }
    if (n <= base){
        //c[0] = a[0]*b[0];
        SM(a,b,c,n);
    } else {
        int a0[n/2], a1[n/2], b0[n/2], b1[n/2], a0_a1[n/2], b0_b1[n/2], x[n], y[n], z[n];
        for (int i = 0; i < n/2; i++){
            a0[i] = a[i];
            a1[i] = a[i + n/2];
            b0[i] = b[i];
            b1[i] = b[i + n/2];
            a0_a1[i] = a0[i] + a1[i];
            b0_b1[i] = b0[i] + b1[i];
        }
        KM(a1, b1, x, n/2, base);
        KM(a0, b0, y, n/2, base);
        KM(a0_a1, b0_b1, z, n/2, base);
        for (int i = 0; i < n; i++){
            c[i] += y[i];
            c[i + n/2] += (z[i] - x[i] - y[i]);
            c[i + n] += x[i];
        }
    }
}

void random_poly(int deg, int bound, int p[]){
    for (int i = 0; i <= deg; i++){
        p[i] = (rand() % (2*bound)) - bound;
    }
}

int main() {
    srand(time(NULL));
    int bound = 100; // polynomial coefficient bounds
    int n = 16; // starting degree - 1
    int base = 16; // switch to schoolbook at degree base-1
    int trials = 10; // number of multiplications and number of random polys
    double kmtimer = 0; // times KM
    double smtimer = 0; // times SM
    clock_t start, end;
    printf("********************************************************************************************\n");
    printf("* Timing schoolbook vs. karatsuba for random integer polynomials with coefficients bounded by %i.\n", bound);
    printf("* %i multiplications each of %i polynomials for 2-power lengths 16 <= n <= 4096.\n", trials, trials);
    printf("* The Karatsuba switches to schoolbook at length n = %i.\n*\n", base);
    while (n < (4097)){
        kmtimer = 0;
        smtimer = 0;
        int a[n], b[n], c1[2*n], c2[2*n];
        for (int i = 0; i < trials; i++){

            random_poly(n-1, bound, a);
            random_poly(n-1, bound, b);

            start = clock();
            for (int j = 0; j < trials; j++){
                KM(a, b, c1, n, base);
            }
            end = clock();
            kmtimer += (double)(end - start);

            start = clock();
            for (int j = 0; j < trials; j++){
                SM(a, b, c2, n);
            }
            end = clock();
            smtimer += (double)(end - start);
        }
        printf("*\tratio smtimer/kmtimer at n = %i: %f\n", n, smtimer/kmtimer);
        n <<=1;
    }
    printf("*\n********************************************************************************************\n");
    return 0;
}
