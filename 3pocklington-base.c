#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <time.h>

// Función floor de logaritmo base 3
int floorlog(int num);

//Test de Pocklington
int pocklingtonTest(mpz_t n, mpz_t p, mpz_t r, mpz_t base, mpz_t criterion, mpz_t mcd, mpz_t k, mpz_t s, mpz_t u, int counter[2]);

// Retorna largo en bits del número
void bitcount(mpz_t n);

// Retorna número impar aleatorio de n-bits
void randomNBitOddNumber(mpz_t num, int nbits, gmp_randstate_t state);


//Test de Miller-Rabin
int millerrabintest(mpz_t n, int exp, mpz_t t, mpz_t nMinus1, mpz_t base);

//Test de Pocklington cuadrático
int pocklingtonTest2(mpz_t n, mpz_t p, mpz_t r, mpz_t base, mpz_t criterion, mpz_t mcd);


int main()
{
    //Declaración de variables
    // n es el candidato
    // n = pr + 1
    // r = 2k

    mpz_t k, n, p, r, randNumb, millerrabin, mrbase, nmenos1, base, criterion, mcd, u, s, psqr;

    //Declaración de estado para rng
    gmp_randstate_t state;

    //Variables de tiempo
    clock_t startTotal, endTotal, startTest, endTest;

    //Variables enteras
    //phi_n son las bases para demostrar que el primer candidato es primo
    int i, nbits, aux, exp, m, proof = 0, mrfactor, errorcount = 0;
    int phi[4] = {2, 3, 5, 7}, stages[10];
    int j, numtests, counter[2] = {0, 0}, tries = 0, initial = 0, expver, sizefail = 0;
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
    mpz_init(psqr);


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

        //Ciclo para generar primos más grandes acotado por el log base 3 de nbits
        aux = floorlog(nbits);

        //Generación de un número primo aleatorio de 27 bits 
        //Empezar con 2^27
        
        randomNBitOddNumber(p, 27, state);

        //El primer primo debe ser primo demostrado
        //Se puede hacer con Miller-Rabin con 4 bases phi_4

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
                    randomNBitOddNumber(p, 27, state);
                    break;
                }
            }
            if(i == 4)
                proof = 1;
        }


        //exponente base 3^3 y multiplicado por 2 para ajustar el tamaño del segundo primo generado
        //expver = 27*2;
        exp = mpz_sizeinbase(p, 2)*2;

        for(i = 3; i < aux; i++)
        {
            mpz_mul(psqr, p, p);
            do
            {
                do
                {
                    do
                    {
                        mpz_rrandomb(k, state, exp-1);
                    }while(mpz_cmp(k, psqr) > 0);

                    mpz_gcd(r, k, p);
                }while(mpz_cmp_ui(r, 1) != 0);

                mpz_mul_ui(r, k, 2);
                mpz_mul(n, r, p);
                mpz_add_ui(n, n, 1);
            }while(!pocklingtonTest(n, p, r, base, criterion, mcd, k, s, u, counter)); 
            
            //Revisa el criterio de tamaño k < p^2+1
            //mpz_mul(r, p, p);
            //mpz_mul_ui(r, r, 2);
            //mpz_add_ui(r,r,1);
            //if(mpz_cmp(k, r) > 0)
            //    printf("ILLEGAL------------------------------------------------------------------\n");
            
            //printf("Bits de k");
            //bitcount(k);
            //printf("Bits de p");
            //bitcount(p);

            mpz_set(p, n);
            //printf("Bits del primo n ");
            //bitcount(n);

            expver = expver*3;
            exp = mpz_sizeinbase(p, 2)*2;
        }

        // Ajuste de bits para el último primo
        //printf("exp = %d\n", exp);
        //printf("expver = %d\n", expver);
        m = nbits - exp/2;

        if(m >= exp >> 1)
        {

        mpz_mul(psqr, p, p);

        do
        {
            do
            {
                do
                {
                    mpz_rrandomb(k, state, m-1);
                }while(mpz_cmp(k, psqr) > 0);
                mpz_gcd(r, k, p);
            }while(mpz_cmp_ui(r, 1) != 0);

            mpz_mul_ui(r, k, 2);
            mpz_mul(n, r, p);
            mpz_add_ui(n, n, 1);

        }while(!pocklingtonTest(n, p, r, base, criterion, mcd, k, s, u, counter));

        if(mpz_cmp(k, psqr) > 0)
            sizefail += 1;
        }
        else
        {
        do
        {
            mpz_rrandomb(k, state, m-1);
            while(mpz_cmp(k, p) > 0)
                mpz_rrandomb(k, state, m-1);
            mpz_mul_ui(r, k, 2);
            mpz_mul(n, r, p);
            mpz_add_ui(n, n, 1);
        }while(!pocklingtonTest2(n, p, r, base, criterion, mcd));

        if(mpz_cmp(k, p) > 0)
            sizefail += 1;

        }

        //printf("Bits de k ");
        //bitcount(k);
        //printf("Bits de p ");
        //bitcount(p);

        mpz_set(p, n);

        endTest = clock();

        avgtime += ((double)endTest - startTest) / CLOCKS_PER_SEC;

        if(mpz_probab_prime_p(p, 15) == 0)
            errorcount += 1;

        //printf("Bits del primo n: ");
        //bitcount(p);
        //printf("Tiempo de búsqueda para número primo de %d bits: %f\n segundos", nbits, ((double)end - start) / CLOCKS_PER_SEC);
    }

    endTotal = clock();


    printf("Tiempo de búsqueda promedio para primo de %d bits: %f segundos\n", nbits, avgtime / numtests);
    printf("Tiempo de ejecución total para %d pruebas: %f segundos\n", numtests, ((double)endTotal - startTotal) / CLOCKS_PER_SEC);
    printf("Errores : %d\n", errorcount);
    printf("Fallas de tamaño: %d\n", sizefail);
    //printf("s impares: %d\n", counter[0]);
    //printf("s pares: %d\n", counter[1]);
    //printf("Intentos de Pocklington: %d\n", tries);

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
    mpz_clear(psqr);

    gmp_randclear(state);

    return 0;
}

int floorlog(int num)
{
    int i = 0;
    while(num >= 3)
    {
        num = num / 3;
        i++;
    }
    return i;
}

// Test de Pocklington con criterio de tamaño cúbico
//

int pocklingtonTest(mpz_t n, mpz_t p, mpz_t r, mpz_t base, mpz_t criterion, mpz_t mcd, mpz_t k, mpz_t s, mpz_t u, int counter[2])
{
    int i, bases[2] = {2, 3};

    //gcd con 105
    //mpz_gcd_ui(mcd, n, 105);
    //if(mpz_cmp_ui(mcd, 1) != 0)
    //{
        //printf("gcd with 105 failed\n");
    //    return 0;
    //}

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
            //gmp_printf("k = %Zd * p + %Zd\n", u, s);

            if(mpz_get_ui(u) & 1)
            {
                //printf("s is odd\n");
                //counter[0] += 1;
                return 1;
            }
            else
            {
                //printf("s is even\n");
                //counter[1] += 1;
                mpz_mul(criterion, s, s);
                //mpz_mul(mcd, s, s);
                mpz_submul_ui(criterion, u, 4);
                //if(mpz_cmp(criterion, mcd) == 0)
                //    printf("u es igual a cero\n");
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

int pocklingtonTest2(mpz_t n, mpz_t p, mpz_t r, mpz_t base, mpz_t criterion, mpz_t mcd)
{
    int i, bases[2] = {2, 3};

    //Se intenta demostrar que se cumple el criterio de Pocklington con las bases 2 y 3
    for (i = 0; i < 2; i++)
    {

        //Primera parte del criterio:
        //base^((n-1)/p)
        mpz_set_ui(base, bases[i]);
        mpz_powm(criterion, base, r, n);
        mpz_powm(criterion, criterion, p, n);
        if(mpz_cmp_ui(criterion, 1) != 0)
        {
            //printf("Pocklington test passed\n");
            return 0;
        }

        mpz_powm(criterion, base, r, n);
        mpz_sub_ui(criterion, criterion, 1);
        mpz_gcd(mcd, criterion, n);

        if(mpz_cmp_ui(mcd, 1) == 0)
        {
            //printf("Primera parte de Pocklington falló \n");
            return 1;
        }

        //Segunda parte del criterio:
        //criterio^p = criterio^(n-1) -> Test de Fermat
        //Si el candidato lo pasa es primo,  no requiere verificación de tamaño por construcción del candidato
        //Si no lo pasa, se intenta con otra base o se rechaza el número
        
   }
   return 0;
}
