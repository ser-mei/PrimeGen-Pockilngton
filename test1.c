#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <time.h>

int millerrabintest(mpz_t n, int exp, mpz_t t, mpz_t nMinus1, mpz_t base);

int main()
{
    mpz_t randnum, gcd, prod, nmenos1, millerrabin, mrbase;
    int i, j, nbits, numtests, compcount1 = 0, compcount2 = 0, compcount3 = 0, mrfactor;
    int bases[26] = {2,3,5,7,11,13,17,19,23,29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101};
    clock_t start, end;
    double avgtime1 = 0, avgtime2 = 0, avgtime3 = 0;

    gmp_randstate_t state;

    mpz_init(randnum);
    mpz_init(gcd);
    mpz_init(prod);
    mpz_init(nmenos1);
    mpz_init(millerrabin);
    mpz_init(mrbase);

    gmp_randinit_mt(state);

    gmp_randseed_ui(state, time(NULL));

    mpz_set_ui(prod, 1);

    for(i = 0; i < 26; i++)
        mpz_mul_ui(prod, prod, bases[i]);
    
    gmp_printf("Product of bases: %Zd\n", prod);

    printf("Enter the number of bits: ");
    scanf("%d", &nbits);

    printf("Enter the number of tests:\n");
    scanf("%d", &numtests);

    for(i = 0; i < numtests; i++)
    {
        mpz_rrandomb(randnum, state, nbits);

        start = clock();

        for(j = 0; j < 26; j++)
        {
            //mpz_mod_ui(gcd, randnum, bases[j]);
            if(mpz_divisible_ui_p(randnum, bases[j]) != 0)
            {
                compcount1 += 1;
                break;
            }
        }

        end = clock();
        avgtime1 += (double)(end - start) / CLOCKS_PER_SEC;

        start = clock();
        mpz_gcd(gcd, randnum, prod);
        if(mpz_cmp_ui(gcd, 1) != 0)
            compcount2 += 1;

        end = clock();
        avgtime2 += (double)(end - start) / CLOCKS_PER_SEC;

        start = clock();

        mpz_sub_ui(nmenos1, randnum, 1);
        mrfactor = mpz_scan1(nmenos1, 0);
        mpz_tdiv_q_2exp(millerrabin, randnum, mrfactor);
 
        mpz_set_ui(mrbase, 2);
        if(millerrabintest(randnum, mrfactor, millerrabin, nmenos1, mrbase) != 1)
            compcount3 += 1;

        end = clock();
        avgtime3 += (double)(end - start) / CLOCKS_PER_SEC;

    }

    printf("Tiempo promedio para trial division: %f\n", avgtime1 / numtests);
    printf("Número de compuestos encontrados: %d\n", compcount1);
    printf("Tiempo promedio para gcd: %f\n", avgtime2 / numtests);
    printf("Número de candidatos descartados: %d\n", compcount2);
    printf("Tiempo promedio para Miller-Rabin: %f\n", avgtime3 / numtests);
    printf("Número de compuestos encontrados: %d\n", compcount3);

    mpz_clear(randnum);
    mpz_clear(gcd);
    mpz_clear(prod);

    gmp_randclear(state);

    return 0;
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