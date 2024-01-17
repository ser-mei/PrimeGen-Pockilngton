#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <time.h>

//Función floor de logaritmo base 2
int floorlog(int num, int initial[1]);

//Test de Pocklington
int pocklingtonTest(mpz_t n, mpz_t p, mpz_t r, mpz_t base, mpz_t criterion, mpz_t mcd);

// Retorna largo en bits del número
void bitcount(mpz_t n);

// Retorna un número impar aleatorio de n-bits
void randomNBitOddNumber(mpz_t num, int nbits, gmp_randstate_t state);

//Test de Miller-Rabin
int millerrabintest(mpz_t n, int exp, mpz_t t, mpz_t nMinus1, mpz_t base);


int main()
{
    //Declaración de variables
    // n es el candidato
    // n = pr + 1
    // r = 2k
    // El resto de variables se utiliza dentro de los tests de primalidad

    mpz_t k, n, p, r, randNumb, millerrabin, mrbase, nmenos1, base, criterion, mcd, filter;

    //Declaración de estado para rng
    gmp_randstate_t state;

    //Variables de tiempo
    clock_t startTotal, endTotal, startTest, endTest;

    //Variables enteras
    //phi_n son las bases para demostrar que el primer candidato es primo con test de Miller-Rabin

    int i, nbits, aux, exp, m, proof = 0, mrfactor, errorcount = 0, initial[1], initialStart;
    int phi[4] = {2, 3, 5, 7};
    int j, numtests, bsize, sizediff = 0;
    double avgtime = 0;

    //Inicialización de variables de gmp
    mpz_init(k);
    mpz_init(n);
    mpz_init(p);
    mpz_init(r);
    mpz_init(millerrabin);
    mpz_init(mrbase);
    mpz_init(nmenos1);
    mpz_init(randNumb);
    mpz_init(base);
    mpz_init(criterion);
    mpz_init(mcd);
    mpz_init(filter);

    //Inicialización de estado para rng
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, time(NULL));

    //Input de número de bits y cantidad de primos para generar
    printf("Enter the number of bits: ");
    scanf("%d", &nbits);

    printf("Enter the number of tests:\n");
    scanf("%d", &numtests);

    startTotal= clock();
    
    for(j = 0; j < numtests; j++)
    {
        startTest = clock();

        aux = floorlog(nbits, initial);
        initialStart = initial[0]+1;

        proof = 0;

        //Generación de un número primo aleatorio de 32 bits 
        //Empezar con candidato de 2^31

        randomNBitOddNumber(p, initialStart, state);

        //El primer primo se demuestra con Miller-Rabin
        //El ciclo termina cuando se encuentra un primo demostrado

        while(!proof)
        {
                mpz_sub_ui(nmenos1, p, 1);
                mrfactor = mpz_scan1(nmenos1, 0);
                mpz_tdiv_q_2exp(millerrabin, p, mrfactor);
 
                for(i = 0; i < 4; i++)
                {
                    mpz_set_ui(mrbase, phi[i]);
                    if(millerrabintest(p, mrfactor, millerrabin, nmenos1, mrbase) != 1)
                    {
                        randomNBitOddNumber(p, initialStart, state);
                        break;
                    }
                }
                if(i == 4)
                    proof = 1;
        }

        //gmp_printf("Primer primo p demostrado = %Zd\n", p);

        //Ciclo para generar primos más grandes desde 2^32 hasta 2^nbits

        bsize = mpz_sizeinbase(p, 2);

        for(i = 1; i < aux; i++)
        {

            //Generación de candidato y preparación de parámetros del test de Pocklington
            do
            {
                do
                {
                    mpz_rrandomb(k, state, bsize-1);
                }while(mpz_cmp(k, p) > 0);

            mpz_mul_ui(r, k, 2);
            mpz_mul(n, r, p);
            mpz_add_ui(n, n, 1);
            }while(!pocklingtonTest(n, p, r, base, criterion, mcd));

            mpz_set(p, n);
            bsize = mpz_sizeinbase(p, 2);
            //printf("bit size: %d\n", bsize);
        }

        //Ajuste de bits
        m = nbits - bsize;

        //printf("m = %d\n", m);

        do
        {
            mpz_rrandomb(k, state, m-1);
            while(mpz_cmp(k, p) > 0)
                mpz_rrandomb(k, state, m-1);
            mpz_mul_ui(r, k, 2);
            mpz_mul(n, r, p);
            mpz_add_ui(n, n, 1);
        }while(!pocklingtonTest(n, p, r, base, criterion, mcd));

        //printf("prime bitsize = %d\n", mpz_sizeinbase(n, 2));

        //mpz_set(p, n);

        endTest = clock();

        avgtime += ((double)endTest - startTest) / CLOCKS_PER_SEC;

        mpz_mul(p,p,p);
        if(mpz_cmp(p,n) > 0)
            sizediff += 1;

        if(mpz_probab_prime_p(n, 15) == 0)
            errorcount += 1;

    
        //bitcount(p);

        //printf("Tiempo de búsqueda para número primo de %d bits: %f\n segundos", nbits, ((double)end - start) / CLOCKS_PER_SEC);
    }

    endTotal = clock();

    printf("Tiempo de búsqueda promedio para primo de %d bits: %f segundos\n", nbits, avgtime / numtests);
    printf("Tiempo de ejecución total para %d pruebas: %f segundos\n", numtests, ((double)endTotal - startTotal) / CLOCKS_PER_SEC);
    printf("Errores: %d\n", errorcount);
    printf("Primos que cumplen el criterio de tamaño: %d\n", sizediff);

    //Liberación de memoria
    mpz_clear(k);
    mpz_clear(n);
    mpz_clear(p);
    mpz_clear(r);
    mpz_clear(randNumb);
    mpz_clear(millerrabin);
    mpz_clear(mrbase);
    mpz_clear(nmenos1);
    mpz_clear(base);
    mpz_clear(criterion);
    mpz_clear(mcd);

    gmp_randclear(state);

    return 0;
}

int floorlog(int num, int initial[1])
{
    int i = 0;
    while(num > 32)
    {
        num = num / 2;
        initial[0] = num;
        i++;
    }
    return i;
}

// Test de Pocklington con criterio de tamaño cuadrático

int pocklingtonTest(mpz_t n, mpz_t p, mpz_t r, mpz_t base, mpz_t criterion, mpz_t mcd)
{
    int i;

    //Se intenta demostrar que se cumple el criterio de Pocklington con la base 2

    //Primera parte del criterio:
    //base^((n-1)/p)
    mpz_set_ui(base,2);
    mpz_powm(criterion, base, r, n);
    mpz_sub_ui(criterion, criterion, 1);
    mpz_gcd(mcd, criterion, n);
    if(mpz_cmp_ui(mcd, 1) != 0)
        return 0;

    //Segunda parte del criterio:
    //criterio^p = criterio^(n-1) -> Test de Fermat
    //Si el candidato lo pasa es primo,  no requiere verificación de tamaño por construcción del candidato
    //Si no lo pasa, se intenta con otra base o se rechaza el número
    mpz_add_ui(criterion, criterion, 1);
    mpz_powm(criterion, criterion, p, n);
    if(mpz_cmp_ui(criterion, 1) == 0)
        return 1;

   return 0;
}

void randomNBitOddNumber(mpz_t num, int nbits, gmp_randstate_t state)
{
    mpz_rrandomb(num, state, nbits);
    if((mpz_get_ui(num) & 1) == 0)
        mpz_add_ui(num, num, 1);
}


// Test de Miller-Rabin
int millerrabintest(mpz_t n, int exp, mpz_t t, mpz_t nMinus1, mpz_t base)
{
    int i;
    mpz_t criterion;

    mpz_init(criterion);

    mpz_powm(criterion, base, t, n);

    if(mpz_cmp_ui(criterion, 1) == 0 || mpz_cmp(criterion, nMinus1) == 0)
    {
        mpz_clear(criterion);
        return 1;
    }

    for(i = 0; i < exp; i++)
    {
        mpz_powm_ui(criterion, criterion, 2, n);
        if(mpz_cmp(criterion, nMinus1) == 0)
        {
            mpz_clear(criterion);
            return 1;
        }
    }
    mpz_clear(criterion);
    return 0;
}
