#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <time.h>

//Función floor de logaritmo base 2
int floorlog(int num);

//Test de Pocklington
int pocklingtonTest(mpz_t n, mpz_t p, mpz_t r, mpz_t base, mpz_t criterion, mpz_t mcd);

// Retorna largo en bits del número
void bitcount(mpz_t n);

// Retorna un número impar aleatorio de n-bits
void randomNBitOddNumber(mpz_t num, int nbits, gmp_randstate_t state);

//Test de Miller-Rabin
int millerrabintest(mpz_t n, int exp, mpz_t t, mpz_t nMinus1, mpz_t base);

//Generador de números primos de nbits en base a Pocklington cuadrático
void primeGenPocklington2(int nbits, mpz_t n, gmp_randstate_t state);


int main()
{

    //Declaración de variables
    mpz_t *n;
    mpz_t mcd, s,t, prime, r, base, criterion;

    //Declaración de estado para rng
    gmp_randstate_t state;

    //Variables de tiempo
    clock_t startTotal, endTotal, startTest, endTest;

    //Variables enteras
    //phi_n son las bases para demostrar que el primer candidato es primo con test de Miller-Rabin

    int i, j, k, nbits, numtests, errorcount = 0, primecount, trycount = 0, flag, totalprimes = 2;
    double avgtime = 0;

    //Arreglo de mpz_t de tamaño 2
    n = (mpz_t *)malloc(2 * sizeof(mpz_t));

    for(i = 0; i < 2; i++)
        mpz_init(n[i]);

    mpz_init(mcd);
    mpz_init(s);
    mpz_init(t);
    mpz_init(prime);
    mpz_init(r);
    mpz_init(base);
    mpz_init(criterion);


    //Inicialización de estado para rng
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, time(NULL));

    //Input de número de bits y cantidad de primos para generar
    printf("Enter the number of bits: ");
    scanf("%d", &nbits);

    printf("Enter the number of tests:\n");
    scanf("%d", &numtests);

    //mpz_init_set_str(filter, "232862364358497360900063316880507363070", 10);

    startTotal= clock();
    
    for(i = 0; i < numtests; i++)
    {
        startTest = clock();

        primeGenPocklington2(nbits/2, n[0], state);
        primeGenPocklington2(nbits/2, n[1], state);
        primecount = 2;

        //for(j = 0; j < primecount-1; j++)
        flag = 0;
        j = 0;
        while(flag != 1)
        {
            for(j = primecount-2; j >= 0; j--)
            {
                mpz_gcdext(mcd, s, t, n[j], n[primecount - 1]);
                //gmp_printf("s = %Zd, t = %Zd, mcd = %Zd\n", s, t, mcd);

                if(mpz_sgn(s) == -1)
                {
                    mpz_mul(prime, s, n[j]);
                    mpz_mul_si(prime, prime, -2);
                    mpz_add_ui(prime, prime, 1);
                    mpz_mul_si(r, s, -2);

                    if(pocklingtonTest(prime, n[j], r, base, criterion, mcd) == 1)
                    {
                        //gmp_printf("Primo--+strong generado: %Zd\n", prime);

                        //mpz_add_ui(mcd, prime, 1);
                        //mpz_mul(r, n[primecount - 1], t);
                        //mpz_mul_ui(r, r, 2);

                        //if(mpz_cmp(mcd, r) == 0)
                        //    printf("N+1 = 2tq -> Primo--+strong\n");
                        //bitcount(prime);
                        flag = 1;
                        break;
                    }
                }
                else
                {
                    mpz_mul(prime, t, n[primecount - 1]);
                    mpz_mul_si(prime, prime, -2);
                    mpz_add_ui(prime, prime, 1);
                    mpz_mul_si(r, t, -2);

                    if(pocklingtonTest(prime, n[primecount - 1], r, base, criterion, mcd) == 1)
                    {
                        //gmp_printf("Primo--+strong generado: %Zd\n", prime);

                        //mpz_add_ui(mcd, prime, 1);
                        //mpz_mul(r, n[j], s);
                        //mpz_mul_ui(r, r, 2);

                        //if(mpz_cmp(mcd, r) == 0)
                        //    printf("N+1 = 2sp -> Primo--+strong\n");

                        //bitcount(prime);
                        flag = 1;
                        break;
                    }
                }
                //else
                //gmp_printf("Primo generado = %Zd\n", n);
    
                //bitcount(p);

                //printf("Tiempo de búsqueda para número primo de %d bits: %f\n segundos", nbits, ((double)end - start) / CLOCKS_PER_SEC);
                
            }
            trycount += 1;
            //printf("N = 2*sp + 1 No es primo \n");
            primecount += 1;
            //if(primecount > 100)
            //{
            //    printf("No se encontró primo--+strong\n");
            //    flag = 1;
            //    break;
            //}

            if(primecount >= totalprimes)
            {
                n = (mpz_t *)realloc(n, primecount * sizeof(mpz_t));
                mpz_init(n[primecount - 1]);
                totalprimes += 1;
            }

            primeGenPocklington2(nbits/2, n[primecount - 1], state);
            //printf("Primo generado\n");

        }
        endTest = clock();
        //printf("primo p--+ generado\n");
        avgtime += ((double)endTest - startTest) / CLOCKS_PER_SEC;
        if(mpz_probab_prime_p(prime, 15) == 0)
            errorcount += 1;
    }

    endTotal = clock();

    printf("Tiempo de búsqueda promedio para primo de %d bits: %f segundos\n", nbits, avgtime / numtests);
    printf("Tiempo de ejecución total para %d pruebas: %f segundos\n", numtests, ((double)endTotal - startTotal) / CLOCKS_PER_SEC);
    printf("Errores: %d\n", errorcount);

    printf("Primos generados para encontrar un primo--+strong: %d\n", trycount);

    //Liberación de memoria
    for(i = 0; i < totalprimes; i++)
        mpz_clear(n[i]);
    
    mpz_clear(mcd);
    mpz_clear(s);
    mpz_clear(t);
    mpz_clear(prime);
    mpz_clear(r);
    mpz_clear(base);
    mpz_clear(criterion);


    free(n);

    gmp_randclear(state);


    return 0;
}

int floorlog(int num)
{
    int i = 0;
    while(num > 1)
    {
        num = num / 2;
        i++;
    }
    return i;
}

// Test de Pocklington con criterio de tamaño cuadrático

int pocklingtonTest(mpz_t n, mpz_t p, mpz_t r, mpz_t base, mpz_t criterion, mpz_t mcd)
{
    int i, bases[2] = {2, 3};

    //gcd con 105 para descartar candidatos rápido
    //mpz_gcd_ui(mcd, n, 105);
    //mpz_gcd(mcd, n, filter);
    //if(mpz_cmp_ui(mcd, 1) != 0)
    //{
        //printf("gcd with 105 failed\n");
    //    return 0;
    //}

    //Se intenta demostrar que se cumple el criterio de Pocklington con las bases 2 y 3
    for (i = 0; i < 2; i++)
    {

        //Primera parte del criterio:
        //base^((n-1)/p)
        mpz_set_ui(base, bases[i]);
        mpz_powm(criterion, base, r, n);
        mpz_sub_ui(criterion, criterion, 1);
        mpz_gcd(mcd, criterion, n);

        if(mpz_cmp_ui(mcd, 1) != 0)
        {
            //printf("Primera parte de Pocklington falló \n");
            return 0;
        }

        //Segunda parte del criterio:
        //criterio^p = criterio^(n-1) -> Test de Fermat
        //Si el candidato lo pasa es primo,  no requiere verificación de tamaño por construcción del candidato
        //Si no lo pasa, se intenta con otra base o se rechaza el número
        mpz_add_ui(criterion, criterion, 1);
        mpz_powm(criterion, criterion, p, n);
        if(mpz_cmp_ui(criterion, 1) == 0)
        {
            //printf("Pocklington test passed\n");
            return 1;
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

void primeGenPocklington2(int nbits, mpz_t n, gmp_randstate_t state)
{
    //Declaración de variables
    // n es el candidato
    // n = pr + 1
    // r = 2k
    // El resto de variables se utiliza dentro de los tests de primalidad

    mpz_t k, p, r, randNumb, millerrabin, mrbase, nmenos1, base, criterion, mcd;


    int i, aux, exp, m, proof = 0;
    int phi[4] = {2, 3, 5, 7};
    unsigned long mrfactor;

    //Inicialización de variables de gmp
    mpz_init(k);
    mpz_init(p);
    mpz_init(r);
    mpz_init(millerrabin);
    mpz_init(mrbase);
    mpz_init(nmenos1);
    mpz_init(randNumb);
    mpz_init(base);
    mpz_init(criterion);
    mpz_init(mcd);



        //Generación de un número primo aleatorio de 32 bits 
        //Empezar con candidato de 2^31

        randomNBitOddNumber(p, 32, state);

        //El primer primo se demuestra con Miller-Rabin
        //El ciclo termina cuando se encuentra un primo demostrado

        while(!proof)
        {
            mpz_gcd_ui(randNumb, p, 15015);
            if(mpz_cmp_ui(randNumb, 1) == 0)
            {
                mpz_sub_ui(nmenos1, p, 1);
                mrfactor = mpz_scan1(nmenos1, 0);
                mpz_tdiv_q_2exp(millerrabin, p, mrfactor);
 
                for(i = 0; i < 4; i++)
                {
                    mpz_set_ui(mrbase, phi[i]);
                    if(millerrabintest(p, mrfactor, millerrabin, nmenos1, mrbase) != 1)
                    {
                        randomNBitOddNumber(p, 32, state);
                        break;
                    }
                }
                if(i == 4)
                    proof = 1;
            }
            else
                randomNBitOddNumber(p, 32, state);
        }

        //gmp_printf("Primer primo p demostrado = %Zd\n", p);

        //Ciclo para generar primos más grandes desde 2^32 hasta 2^nbits

        aux = floorlog(nbits);

        for(i = 5; i < aux; i++)
        {
            exp = 1 << i;
            //printf("exponente = 2^%d = %d\n", i, exp);

            //Generación de candidato y preparación de parámetros del test de Pocklington

            do
            {
                do
                {
                
                    //randomNBitOddNumber(k, exp, state);
                    mpz_rrandomb(k, state, exp-1);
                    mpz_mul_ui(r, k, 2);
                    mpz_mul(n, r, p);
                    mpz_add_ui(n, n, 1);
                }while(mpz_probab_prime_p(n, 1) == 0);
            }while(!pocklingtonTest(n, p, r, base, criterion, mcd));
            //{
            //    randomNBitOddNumber(k, exp, state);
            //    mpz_mul_ui(r, k, 2);
            //    mpz_mul(n, r, p);
            //    mpz_add_ui(n, n, 1);
            //}
            //gmp_printf("prime n = %Zd = %Zd * %Zd + 1\n", n, p, r);
            mpz_set(p, n);
        }

        //bitcount(n);

        //Ajuste de bits

        m = 1 << aux;
        //printf("logfloor de nbits = %d\n", aux);
        //printf("m = %d\n", m);
        m = nbits - m;
        //printf("m = %d\n", m);

        do
        {
            do
            {
            
                //randomNBitOddNumber(k, m, state);
                mpz_rrandomb(k, state, m-1);
                mpz_mul_ui(r, k, 2);
                mpz_mul(n, r, p);
                mpz_add_ui(n, n, 1);
            } while (mpz_probab_prime_p(n, 1) == 0);
        }while(!pocklingtonTest(n, p, r, base, criterion, mcd));

        //mpz_set(p, n);
        //gmp_printf("Número primo n = %Zd \n", p);

    //Liberación de memoria
    mpz_clear(k);
    mpz_clear(p);
    mpz_clear(r);
    mpz_clear(randNumb);
    mpz_clear(millerrabin);
    mpz_clear(mrbase);
    mpz_clear(nmenos1);
    mpz_clear(base);
    mpz_clear(criterion);
    mpz_clear(mcd);
}