#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <time.h>

int main()
{
    mpz_t num;
    int nbits, ntests, minsize, maxsize, mrtests;
    double avgtime = 0;
    clock_t start, end;
    FILE *f;

    gmp_randstate_t state;

    mpz_init(num);

    printf("Enter the initial number of bits: ");
    scanf("%d", &minsize);

    printf("Enter the final number of bits: ");
    scanf("%d", &maxsize);

    printf("Enter the number of primes to generate: ");
    scanf("%d", &ntests);

    printf("Enter the number of Miller-Rabin tests: ");
    scanf("%d", &mrtests);

    gmp_randinit_mt(state);
    //gmp_randseed_ui(state, time(NULL));
    gmp_randseed_ui(state, 1234567890);

    f = fopen("GMP-Generator-Times.txt", "w");

    if(f == NULL)
    {
        printf("Error opening file!\n");
        return 1;
    }
    else
    {
        fprintf(f, "Generador de GMP, optmimizado con Baillie.PSW - Pruebas con %d números generados, utiilzando %d cantidad de tests de  Miller-Rabin\n", ntests, mrtests);
        fprintf(f, "nbits\t ktests\t avgtime\n");
    }

    for(nbits = minsize; nbits <= maxsize; nbits += 500)
    {
        avgtime = 0;

        for(int i = 0; i < ntests; i++)
        {
            start = clock();

            mpz_rrandomb(num, state, nbits);
            if((mpz_get_ui(num) & 1) == 0)
                mpz_add_ui(num, num, 1);

            while(mpz_probab_prime_p(num, mrtests) == 0)
            {
                mpz_rrandomb(num, state, nbits);
                if((mpz_get_ui(num) & 1) == 0)
                    mpz_add_ui(num, num, 1);
            }
            end = clock();

            avgtime += (double)(end - start) / CLOCKS_PER_SEC;
        }

        printf("Average time for a %d bits prime: %f segundos \n", nbits, avgtime/ntests);

        fprintf(f, "%d\t %d\t %f\n", nbits, ntests, avgtime/ntests);

    }

    mpz_clear(num);

    gmp_randclear(state);

    fclose(f);

    return 0;
}
