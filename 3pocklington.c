#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <time.h>

// Función floor de logaritmo base 3
int floorlog(int num);

//Test de Pocklington
int pocklingtonTest(mpz_t n, mpz_t p, mpz_t r, mpz_t base, mpz_t criterion, mpz_t mcd, mpz_t k, mpz_t s, mpz_t u);

// Retorna largo en bits del número
void bitcount(mpz_t n);

// Retorna número impar aleatorio de n-bits
void randomNBitOddNumber(mpz_t num, int nbits, gmp_randstate_t state);

//Test de Miller-Rabin
int millerrabintest(mpz_t n, int exp, mpz_t t, mpz_t nMinus1, mpz_t base);


int main()
{
    //Declaración de variables
    // n es el candidato
    // n = pr + 1
    // r = 2k

    mpz_t k, n, p, r, randNumb, millerrabin, mrbase, nmenos1, base, criterion, mcd, u, s;

    //Declaración de estado para rng
    gmp_randstate_t state;

    //Variables de tiempo
    clock_t startTotal, endTotal, startTest, endTest;

    //Variables enteras
    //phi_n son las bases para demostrar que el primer candidato es primo
    int i, nbits, aux, exp, m, proof = 0, mrfactor, errorcount = 0;
    int phi[4] = {2, 3, 5, 7};
    int j, numtests;
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
    mpz_init(u);
    mpz_init(s);


    //Inicialización de estado para rng
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, time(NULL));

    //Input de número de bits
    printf("Enter the number of bits: ");
    scanf("%d", &nbits);

    printf("Enter the number of tests:\n");
    scanf("%d", &numtests);

    startTotal= clock();

    for(j = 0; j < numtests; j++)
    {
        startTest = clock();
        proof = 0;

        //Generación de un número primo aleatorio de 32 bits 
        //Empezar con 2^31

        randomNBitOddNumber(p, 27, state);
        //gmp_printf("Candidate p = %Zd\n", p);

        //El primer primo debe ser primo demostrado
        //Se puede hacer con Miller-Rabin con 4 bases phi_4

        while(!proof)
        {
            //gmp_printf("Candidate p = %Zd\n", p);
            mpz_gcd_ui(randNumb, p, 105);
            if(mpz_cmp_ui(randNumb, 1) == 0)
            {
                mpz_sub_ui(nmenos1, p, 1);
                mrfactor = mpz_scan1(nmenos1, 0);
                mpz_tdiv_q_2exp(millerrabin, p, mrfactor);
 
                for(i=0; i<4; i++)
                {
                    mpz_set_ui(mrbase, phi[i]);
                    if(millerrabintest(p, mrfactor, millerrabin, nmenos1, mrbase) != 1)
                    {
                        randomNBitOddNumber(p, 27, state);
                        break;
                    }
                }
                if(i == 4)
                    proof = 1;
            }
            else
                randomNBitOddNumber(p, 27, state);
        }

        //gmp_printf("Primer primo p demostrado = %Zd\n", p);

        //Ciclo para generar primos más grandes acotado por el log base 3 de nbits
        aux = floorlog(nbits);
        //printf("logfloor base 3 de nbits = %d\n", aux);

        //exponente base 3^3 y multiplicado por 2 para ajustar el tamaño del segundo primo generado
        exp = 9*2;

        for(i = 3; i < aux-1; i++)
        {
            exp = exp*3;

            //printf("exponente = 3^%d = %d\n", i, exp);

            do
            {
                randomNBitOddNumber(k, exp, state);
                mpz_gcd(r, k, p);
            }while(mpz_cmp_ui(r, 1) != 0);
        
            //Revisa criterio de tamaño k < p^2 + 1
            //mpz_mul(r, p, p);
            //mpz_mul_ui(r, r, 2);
            //mpz_add_ui(r,r,2);
            //if(mpz_cmp(k, r) < 0)
            //    printf("LEGAL\n");
            //else
            //    printf("ILLEGAL------------------------------------------------------------------\n");

            mpz_mul_ui(r, k, 2);
            mpz_mul(n, r, p);
            mpz_add_ui(n, n, 1);


            while(!pocklingtonTest(n, p, r, base, criterion, mcd, k, s, u))
            {
                //printf("Cant find prime\n");
                //mpz_urandomb(k, state, exp);

                do
                {
                    randomNBitOddNumber(k, exp, state);
                    mpz_gcd(r, k, p);
                }while(mpz_cmp_ui(r, 1) != 0);

                mpz_mul_ui(r, k, 2);
                mpz_mul(n, r, p);
                mpz_add_ui(n, n, 1);
            }

            //mpz_mul(r, p, p);
            //mpz_mul_ui(r, r, 2);
            //mpz_add_ui(r,r,2);
            //if(mpz_cmp(k, r) < 0)
            //    printf("LEGAL\n");
            //else
            //    printf("ILLEGAL------------------------------------------------------------------\n");


            //gmp_printf("prime n = %Zd = %Zd * %Zd + 1\n", n, p, r);
            mpz_set(p, n);
            //bitcount(n);
        }

        // Ajuste de bits para el último primo
        m = exp*3/2;
        //printf("logfloor de nbits = %d\n", aux);
        //printf("m = %d\n", m);
        m = nbits - m;
        //printf("m = %d\n", m);

        do
        {
            do
            {
                randomNBitOddNumber(k, m, state);
                mpz_gcd(r, k, p);
            }while(mpz_cmp_ui(r, 1) != 0);
            mpz_mul_ui(r, k, 2);
            mpz_mul(n, r, p);
            mpz_add_ui(n, n, 1);
        }while(!pocklingtonTest(n, p, r, base, criterion, mcd, k, s, u));

        //mpz_mul(r, p, p);
        //mpz_mul_ui(r, r, 2);
        //mpz_add_ui(r,r,2);
        //if(mpz_cmp(k, r) < 0)
        //    printf("LEGAL\n");
        //else
        //    printf("ILLEGAL------------------------------------------------------------------\n");

        mpz_set(p, n);



        //gmp_printf("Número primo n = %Zd \n", p);

        endTest = clock();

        avgtime += ((double)endTest - startTest) / CLOCKS_PER_SEC;

        if(mpz_probab_prime_p(p, 15) == 0)
            errorcount += 1;

        //bitcount(p);
        //printf("Tiempo de búsqueda para número primo de %d bits: %f\n segundos", nbits, ((double)end - start) / CLOCKS_PER_SEC);
    }

    endTotal = clock();


    printf("Tiempo de búsqueda promedio para primo de %d bits: %f segundos\n", nbits, avgtime / numtests);
    printf("Tiempo de ejecución total para %d pruebas: %f segundos\n", numtests, ((double)endTotal - startTotal) / CLOCKS_PER_SEC);
    printf("Errores : %d\n", errorcount);

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
    mpz_clear(u);
    mpz_clear(s);

    gmp_randclear(state);

    return 0;
}

int floorlog(int num)
{
    int i = 0;
    while(num >= 1)
    {
        num = num / 3;
        i++;
    }
    return i;
}

// Test de Pocklington con criterio de tamaño cúbico
//

int pocklingtonTest(mpz_t n, mpz_t p, mpz_t r, mpz_t base, mpz_t criterion, mpz_t mcd, mpz_t k, mpz_t s, mpz_t u)
{
    int i, bases[2] = {2, 3};

    //gcd con 105
    mpz_gcd_ui(mcd, n, 105);
    if(mpz_cmp_ui(mcd, 1) != 0)
    {
        //printf("gcd with 105 failed\n");
        return 0;
    }

    for (i = 0; i < 2; i++)
    {
        mpz_set_ui(base, bases[i]);
        mpz_powm(criterion, base, r, n);
        mpz_sub_ui(criterion, criterion, 1);
        mpz_gcd(mcd, criterion, n);

        if(mpz_cmp_ui(mcd, 1) != 0)
        {
            //printf("Pocklington test failed\n");
            return 0;
        }

        mpz_add_ui(criterion, criterion, 1);
        mpz_powm(criterion, criterion, p, n);
        if(mpz_cmp_ui(criterion, 1) == 0)
        {
            //printf("Pocklington test passed\n");

            //Revisión de tercera condición del criterio:
            // r = u*p + s, con u impar o que se cumpla que u es par y s^2-4u no es cuadrado perfecto
            mpz_fdiv_qr(u, s, k, p);

            if(mpz_get_ui(s) & 1)
            {
                //printf("s is odd\n");
                return 1;
            }
            else
            {
                //printf("s is even\n");
                mpz_mul(criterion, s, s);
                mpz_submul_ui(criterion, u, 4);
                if(mpz_perfect_square_p(criterion) == 0)
                    return 1;
                else
                    return 0;
            }
        }
    }
    return 0;
}

void bitcount(mpz_t n)
{
    mpz_t aux;
    mpz_init(aux);
    mpz_set(aux, n);
    int i = 0;
    while(mpz_cmp_ui(aux, 0) > 0)
    {
        mpz_fdiv_q_ui(aux, aux, 2);
        i++;
    }
    printf("bits = %d\n", i);

    mpz_clear(aux);
}

void randomNBitOddNumber(mpz_t num, int nbits, gmp_randstate_t state)
{
    mpz_rrandomb(num, state, nbits);
    if((mpz_get_ui(num) & 1) == 0)
        mpz_add_ui(num, num, 1);
}


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