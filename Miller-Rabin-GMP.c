#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <time.h>

int main()
{
    mpz_t num;
    int nbits, ntests;
    double avgtime;
    clock_t start, end;
    gmp_randstate_t state;

    mpz_init(num);

    printf("Enter the number of bits: ");
    scanf("%d", &nbits);

    printf("Enter the number of tests: ");
    scanf("%d", &ntests);

    gmp_randinit_mt(state);
    gmp_randseed_ui(state, time(NULL));

    for(int i = 0; i < ntests; i++)
    {
        start = clock();

        mpz_rrandomb(num, state, nbits);
        if((mpz_get_ui(num) & 1) == 0)
            mpz_add_ui(num, num, 1);

        while(mpz_probab_prime_p(num, 25) == 0)
        {
            mpz_rrandomb(num, state, nbits);
            if((mpz_get_ui(num) & 1) == 0)
                mpz_add_ui(num, num, 1);
        }
        end = clock();

        avgtime += (double)(end - start) / CLOCKS_PER_SEC;
    }

    printf("Average time for a %d bits prime: %f segundos \n", nbits, avgtime/ntests);

    mpz_clear(num);

    gmp_randclear(state);

    return 0;
}
