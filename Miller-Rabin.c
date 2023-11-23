#include <stdio.h>
#include <gmp.h>
#include <time.h>


int millerrabintest(mpz_t n, int exp, mpz_t t, mpz_t nMinus1, mpz_t base, mpz_t criterion);

void nBitNumber(mpz_t rand_num, int nbits, gmp_randstate_t state);



int main()
{

    //Variables enteras
    //flag se utiliza en el test de Miller-Rabin y flag2 para saber si se encontró un primo probable fuerte
    int nbits, trials, flag = 1, flag2 = 0, power, numtest;
    double avgtime = 0;
    int count = 0;

    clock_t start, end;

    mpz_t rand_num, rand_base, n_minus_one, factor, criterion, mcd;

    gmp_randstate_t state;

    // Initialize random number and state
    mpz_init(rand_num);
    mpz_init(rand_base);
    mpz_init(n_minus_one);
    mpz_init(factor);
    mpz_init(criterion);
    mpz_init(mcd);

    gmp_randinit_mt(state);

    // Seed the random state
    gmp_randseed_ui(state, time(NULL));


    // Get the number of bits as input
    printf("Enter the number of bits: ");
    scanf("%d", &nbits);

    printf("Enter the number of tests: ");
    scanf("%d", &numtest);


    trials = 25;

    for(int j = 0; j < numtest; j++)
    {

    start = clock();

    flag2 = 0;

    while(!flag2)
    {

        flag = 1;

        // Generate random number with n-bits
        //mpz_urandomb(rand_num, state, nbits);
        nBitNumber(rand_num, nbits, state);

        //if the number is even, add 1 to make it odd
        //if ((mpz_get_ui(rand_num) & 1) == 0) 
        //    mpz_add_ui(rand_num, rand_num, 1);

        mpz_gcd_ui(mcd, rand_num, 105);   

    if(mpz_cmp_ui(mcd, 1) != 0)
        flag = 0;
    else
    {
    //printf("Enter the number of trials: ");
    //scanf("%d", &trials);

    //Declaración de n-1
    mpz_sub_ui(n_minus_one, rand_num, 1);

    //La posición del primer 1 menos significante es la cantidad de 2 que se pueden factorizar 
    power = mpz_scan1(n_minus_one, 0);
    mpz_tdiv_q_2exp(factor, n_minus_one, power);

    for(int i = 0; i < trials; i++)
    {
        mpz_urandomb(rand_base, state, nbits);
        mpz_add_ui(rand_base, rand_base, 2);

        if(millerrabintest(rand_num, power, factor, n_minus_one, rand_base, criterion) == 0)
        {
            flag = 0;
            break;
        }
    }

    //if(mpz_probab_prime_p(rand_num, trials) > 0)
    //    flag = 1;
    //else
    //    flag = 0;

    }

        if(flag == 1)
        {
            //gmp_printf("Probable prime number: %Zd\n", rand_num);
            flag2 = 1;
        }

    }

    end = clock();
    avgtime += (double)(end - start) / CLOCKS_PER_SEC;

    if(mpz_probab_prime_p(rand_num, 25) == 0)
        count += 1;
//    else
//        printf("GMP says it's composite\n");
    }

    //printf("Time: %f segundos \n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("Average time for a %d bits prime: %f segundos \n", nbits, avgtime/numtest);
    printf("Number of errors: %d\n", count);



    // Clear memory
    mpz_clear(rand_num);
    mpz_clear(rand_base);
    mpz_clear(n_minus_one);
    mpz_clear(factor);
    mpz_clear(criterion);
    mpz_clear(mcd);

    gmp_randclear(state);

    return 0;
}


int millerrabintest(mpz_t n, int exp, mpz_t t, mpz_t nMinus1, mpz_t base, mpz_t criterion)
{
    int i;

    mpz_powm(criterion, base, t, n);

    if(mpz_cmp_ui(criterion, 1) == 0 || mpz_cmp(criterion, nMinus1) == 0)
    {
        return 1;
    }

    for(i = 0; i < exp; i++)
    {
        mpz_powm_ui(criterion, criterion, 2, n);
        if(mpz_cmp(criterion, nMinus1) == 0)
        {
            return 1;
        }
    }
    return 0;
}

void nBitNumber(mpz_t rand_num, int nbits, gmp_randstate_t state)
{
    mpz_rrandomb(rand_num, state, nbits);
    if((mpz_get_ui(rand_num) & 1) == 0)
        mpz_add_ui(rand_num, rand_num, 1);
}
