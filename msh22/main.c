#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <fenv.h>

double somatorio(double* pi, double tolerance, unsigned int* n, double* aproximate_absolut_error, long long int* flops){
    double last_value = 0;
    long double fact_up_double = 1; 
    long double fact_down_double = 1; 
    unsigned long int fact_up_int = 1;
    unsigned long int fact_down_int = 1;
    double up_part = 0;

    do{
        if(*n == 0){
            *pi = 2.0;
        }
        else{
            if(*n < 10){
                up_part = ((1ULL << *n) * fact_up_int*fact_up_int);
                *pi += 2 * (up_part / fact_down_int);
                *flops += 2;
            }else{
                up_part = ((1ULL << *n) * fact_up_double*fact_up_double);
                *pi += 2 * (up_part / fact_down_double);
                *flops += 4;
            }
        }
        // não considerei essa operação como um flop, pois não faz parte diretamente da aproximação do pi
        *aproximate_absolut_error = fabs(*pi - last_value); 
        if(*aproximate_absolut_error < tolerance){
            *n += 1; // como começa em 0, para ficar mais preciso o número de iterações
            break;
        }
        *n += 1;
        last_value = *pi;
        if(*n < 10){
            fact_up_int = fact_up_int * (*n);
            fact_down_int = fact_down_int * ((2 * *n + 1) * (2 * *n));
            fact_up_double = (long double)fact_up_int;
            fact_down_double = (long double)fact_down_int;
        }else{
            fact_up_double = fact_up_double * (*n);
            fact_down_double = fact_down_double * ((2 * *n + 1) * (2 * *n));
            *flops += 2;
        }
    }while(1);
}

void pi_calc(double* pi, double tolerance, unsigned int* n, double* aproximate_absolut_error, long long int* flops){

    fesetround(FE_DOWNWARD);
    somatorio(&pi[0], tolerance, n, aproximate_absolut_error, flops);
    fesetround(FE_UPWARD);
    *aproximate_absolut_error = 0;
    *n = 0;
    *flops = 0;
    somatorio(&pi[1], tolerance, n, aproximate_absolut_error, flops);
}

int main(int argc, char *argv[]){
    double tolerance = 0;
    int n = 0; // número de iterações
    double aproximate_absolut_error = 0;
    double exact_absolut_error = 0;
    double pi[2] = {0, 0};
    int64_t *ptr_pi_down = NULL;
    int64_t *ptr_pi_up = NULL;
    int64_t *ptr_abe = NULL; // aproximate_absolut_error
    int64_t *ptr_ebe = NULL; // exact_absolut_error
    int ulps = 0;
    long long int flops = 0; 

    scanf("%lf", &tolerance);

    pi_calc(pi, tolerance, &n, &aproximate_absolut_error, &flops);
    exact_absolut_error = fabs(M_PI - pi[1]);

    // convertendo para int64_t para poder representar em hexadecimal e calcular ulp
    ptr_pi_down = (int64_t *) &pi[0];
    ptr_pi_up = (int64_t *) &pi[1];
    ptr_abe = (int64_t *) &aproximate_absolut_error;
    ptr_ebe = (int64_t *) &exact_absolut_error;
    ulps = abs(*ptr_pi_up - *ptr_pi_down) -1;

    printf("%d\n", n);
    printf("%.15e %llx\n", aproximate_absolut_error, *ptr_abe);

    printf("%.15e %llx\n", exact_absolut_error, *ptr_ebe);
    
    printf("%.15e %llx\n", pi[0], *ptr_pi_down); // pi_down

    printf("%.15e %llx\n", pi[1], *ptr_pi_up); // pi_up

    printf("%d\n", ulps);
    printf("%lld\n", flops);

    return 0;
}
