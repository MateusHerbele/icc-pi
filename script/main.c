#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fenv.h>


double factorial(int n, unsigned int* flops){
    if(n == 0)
        return 1;
    else{
        *flops += 1;
        return n * factorial(n - 1, flops);
    }
}

double exponentiation(int base, int exponent, unsigned int* flops){
    if(exponent == 1)
        return base;
    else{
        *flops += 1;
        return base * exponentiation(base, exponent - 1, flops);
    }
}
double somatorio(double* pi, double tolerance, unsigned int* n, double* aproximate_absolut_error, unsigned int* ulps, unsigned int* flops){
    double last_value = 0;
    double exp = 1;
    double fact_up = 1; // Fatorial de cima da fração
    double fact_down = 1; // Fatorial de baixo da fração
    double up_part = 0; // Parte de cima da fração
    
    // tentar fazer as variáveis já estarem todas alinhadas a partir da execução 1 (2)

    // inicial erro absoluto
    printf("ABS INICIAL: %.15e\n", *aproximate_absolut_error);
    do{
        if(*n == 0){
            *pi = 2.0;
        }
        else{
            
            up_part = (exp * (fact_up*fact_up));
            *pi += 2.0 * (up_part / fact_down);
            *flops += 3;
        }
        *aproximate_absolut_error = fabs(*pi - last_value); // talvez conte como 1 operação de ponto flutuante
        *n += 1;
        last_value = *pi;
        printf("Iteração: %d, PI: %.15e, ABS: %.15e\n", *n, *pi, *aproximate_absolut_error);
        exp *= 2;
        fact_up = fact_up * (*n);
        fact_down = fact_down * ((2 * *n + 1) * (2 * *n));
        *flops += 3; // considerando abe
    }while(*aproximate_absolut_error >= tolerance);
    
}

void pi_calc(double* pi, double tolerance, unsigned int* n, double* aproximate_absolut_error, unsigned int* ulps, unsigned int* flops){

    fesetround(FE_DOWNWARD);
    somatorio(&pi[0], tolerance, n, aproximate_absolut_error, ulps, flops);
    fesetround(FE_UPWARD);
    *aproximate_absolut_error = 0;
    *n = 0;
    somatorio(&pi[1], tolerance, n, aproximate_absolut_error, ulps, flops);
}

void print_hex(double* pi){

    unsigned long long *ptr = (unsigned long long*) pi;

    printf("%llx\n", *ptr); 

}


int main(int argc, char *argv[]){
/*
IMPRIMIR: 
Número de interações
Erro absoluto aproximado obtido
Erro absoluto exato
Apróximado de pi arredondodamento para baixo
Apróximado de pi arredondodamento para cima
Diferença de ULPs entre arredoamento p/baixo e p/cima
Número de operações de ponto flutuante
*/
    double tolerance = 0;
    unsigned int n = 0; 
    double aproximate_absolut_error = 0;
    double exact_absolut_error = 0;
    double pi[2] = {0, 0};
    unsigned long long *ptr_pi_down = NULL;
    unsigned long long *ptr_pi_up = NULL;
    unsigned int ulps = 0;
    unsigned int flops = 0; 

    if(argc > 1){
        printf("INPUT: %s!\n", argv[1]);
    }else{
        printf("No input detected!\n");
        return 1;
    }

    tolerance = strtod(argv[1], NULL);
    printf("Tolerance: %.15e\n", tolerance);
    pi_calc(pi, tolerance, &n, &aproximate_absolut_error, &ulps, &flops);
    exact_absolut_error = fabs(M_PI - pi[1]);


    ptr_pi_down = (unsigned long long*) &pi[0];
    ptr_pi_up = (unsigned long long*) &pi[1];
    ulps = *ptr_pi_up - *ptr_pi_down;

    printf("ITERAÇÕES: %u\n", n);
    printf("APROXIMATE ABSOLUT ERROR: %.15e ", aproximate_absolut_error);
    print_hex(&aproximate_absolut_error);
    printf("EXACT ABSOLUT ERRO:       %.15e ", exact_absolut_error);
    print_hex(&exact_absolut_error);
    printf("Pi down : %.15e ", pi[0]); // pi_down
    print_hex(&pi[0]);
    printf("Pi up:    %.15e ", pi[1]); // pi_up
    print_hex(&pi[1]);  
    printf("ULPS: %u\n", ulps);
    printf("FLOPS:    %u\n", flops);




    printf("%.15e\n", M_PI);

    return 0;
}
