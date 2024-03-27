#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <time.h>

//Función floor de logaritmo base 2
int floorlog(int num, int *initial);

//Test de Pocklington
int pocklingtonTest(mpz_t n, mpz_t p, mpz_t r, mpz_t base, mpz_t criterion, mpz_t mcd);

// Retorna un número impar aleatorio de n-bits
void randomNBitOddNumber(mpz_t num, int nbits, gmp_randstate_t state);

//Test de Miller-Rabin
int millerrabintest(mpz_t n, int exp, mpz_t t, mpz_t nMinus1, mpz_t base);

//Generador de números primos de nbits en base a Pocklington cuadrático
void primeGenPocklington2(int nbits, mpz_t n, gmp_randstate_t state, mpz_t k, mpz_t p, mpz_t r, mpz_t randNumb, mpz_t millerrabin, mpz_t mrbase, mpz_t nmenos1, mpz_t base, mpz_t criterion, mpz_t mcd, mpz_t filter);

int compositefilter(mpz_srcptr n);


int main()
{

    //Declaración de variables
    mpz_t *n;
    mpz_t mcd, s,t, prime, r, base, criterion;
    mpz_t k, p, randNumb, millerrabin, mrbase, nmenos1, filter;

    int ln2;

    //Declaración de estado para rng
    gmp_randstate_t state;

    //Variables de tiempo
    clock_t startTotal, endTotal, startTest, endTest;

    //Variables enteras
    //phi_n son las bases para demostrar que el primer candidato es primo con test de Miller-Rabin

    int i, j, nbits, numtests, errorcount = 0, primecount, trycount = 0, flag, totalprimes = 2, euclidCount = 0, primeDiff = 0, primeDiffCount = 0;
    double avgtime = 0, avgtries, avgprimes;

    //FILE *f;

    mpz_init(mcd);
    mpz_init(s);
    mpz_init(t);
    mpz_init(prime);
    mpz_init(r);
    mpz_init(base);
    mpz_init(criterion);
    mpz_init(k);
    mpz_init(p);
    mpz_init(randNumb);
    mpz_init(millerrabin);
    mpz_init(mrbase);
    mpz_init(nmenos1);
    mpz_init(filter);

    //Inicialización de estado para rng
    gmp_randinit_mt(state);
    //gmp_randseed_ui(state, time(NULL));
    gmp_randseed_ui(state, 1234567890);

    //Input de número de bits y cantidad de primos para generar
    printf("Enter the number of bits: ");
    scanf("%d", &nbits);

    printf("Enter the number of tests:\n");
    scanf("%d", &numtests);

    mpz_set_str(filter, "116431182179248680450031658440253681535", 10);

    //Arreglo de mpz_t de tamaño nbits 
    n = (mpz_t *)malloc(nbits * sizeof(mpz_t));

    for(i = 0; i < nbits; i++)
        mpz_init(n[i]);


    startTotal= clock();
    
    for(i = 0; i < numtests; i++)
    {
        startTest = clock();

        primeGenPocklington2(nbits/2, n[0], state, k, p, r, randNumb, millerrabin, mrbase, nmenos1, base, criterion, mcd, filter);
        primeGenPocklington2(nbits/2, n[1], state, k, p, r, randNumb, millerrabin, mrbase, nmenos1, base, criterion, mcd, filter);
        primecount = 2;
        trycount = 0;

        //for(j = 0; j < primecount-1; j++)
        flag = 0;
        j = 0;

        while(flag != 1)
        {
            for(j = primecount-2; j >= 0; j--)
            {
                //primeDiff = mpz_cmp(n[j], n[primecount - 1]);
                mpz_gcdext(mcd, s, t, n[j], n[primecount - 1]);
                //gmp_printf("s = %Zd, t = %Zd, mcd = %Zd\n", s, t, mcd);

                if(mpz_sgn(s) == -1)
                {
                    mpz_mul(prime, s, n[j]);
                    mpz_mul_si(prime, prime, -2);
                    mpz_add_ui(prime, prime, 1);
                    mpz_mul_si(r, s, -2);

                    mpz_gcd(mcd, prime, filter);

                    if(mpz_cmp_ui(mcd, 1) == 0)
                    {

                        if(compositefilter(prime) == 1)
                        {
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
                                //mpz_abs(s,s);
                                //if(mpz_cmp(n[j], s) > 0)
                                //    euclidCount += 1;
                                //if(primeDiff > 0)
                                //    primeDiffCount += 1;
                                break;
                            }
                        }
                   }
                }
                else
                {
                    mpz_mul(prime, t, n[primecount - 1]);
                    mpz_mul_si(prime, prime, -2);
                    mpz_add_ui(prime, prime, 1);
                    mpz_mul_si(r, t, -2);

                    mpz_gcd(mcd, prime, filter);

                    if(mpz_cmp_ui(mcd, 1) == 0)
                    {

                        if(compositefilter(prime) != 0)
                        {
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
                            //mpz_abs(t,t);
                            //if(mpz_cmp(n[primecount-1], t) > 0)
                            //    euclidCount += 1;
                            //if(primeDiff < 0)
                            //    primeDiffCount += 1;
                                break;
                            }
                        }
                    }
                //else
                //gmp_printf("Primo generado = %Zd\n", n);
    
                //bitcount(p);

                //printf("Tiempo de búsqueda para número primo de %d bits: %f\n segundos", nbits, ((double)end - start) / CLOCKS_PER_SEC);
                
                }

                trycount += 1;
            }
            //printf("N = 2*sp + 1 No es primo \n");
            primecount += 1;
            //if(primecount > 100)
            //{
            //    printf("No se encontró primo--+strong\n");
            //    flag = 1;
            //    break;
            //}

            if(primecount > totalprimes)
            {
                //n = (mpz_t *)realloc(n, primecount * sizeof(mpz_t));
                //mpz_init(n[primecount - 1]);
                totalprimes += 1;
            }
            if(totalprimes > nbits)
            {
                n = (mpz_t *)realloc(n, primecount * sizeof(mpz_t));
                mpz_init(n[primecount - 1]);
            }

            primeGenPocklington2(nbits/2, n[primecount - 1], state, k, p, r, randNumb, millerrabin, mrbase, nmenos1, base, criterion, mcd, filter);
            //printf("Primo generado\n");

        }

        endTest = clock();
        //printf("primo p--+ generado\n");
        
        avgtime += ((double)endTest - startTest) / CLOCKS_PER_SEC;
        avgtries += (double)trycount / numtests;
        avgprimes += (double)primecount / numtests;

        if(mpz_probab_prime_p(prime, 15) == 0)
            errorcount += 1;
    }

    endTotal = clock();

    printf("Tiempo de búsqueda promedio para primo de %d bits: %f segundos\n", nbits, avgtime / numtests);
    printf("Tiempo de ejecución total para %d pruebas: %f segundos\n", numtests, ((double)endTotal - startTotal) / CLOCKS_PER_SEC);
    printf("Errores: %d\n", errorcount);
    printf("Promedio de candidatos: %f\n", avgtries);
    printf("Promedio de primos generados: %f\n", avgprimes);
    printf("Primos máximos generados: %d\n", totalprimes);
    printf("--------------------------------------------------------------------------------------\n");

    //printf("Primos generados para encontrar un primo--+strong: %d y %d\n", trycount, totalprimes);

    //printf("Euclid: %d y Correspondencia de diferencia de tamaño: %d\n", euclidCount, primeDiffCount);

    //Liberación de memoria
    for(i = 0; i < nbits; i++)
        mpz_clear(n[i]);

    if(totalprimes > nbits)
    {
        for(i = nbits; i < totalprimes; i++)
            mpz_clear(n[i]);
    }
    
    mpz_clear(mcd);
    mpz_clear(s);
    mpz_clear(t);
    mpz_clear(prime);
    mpz_clear(r);
    mpz_clear(base);
    mpz_clear(criterion);
    mpz_clear(k);
    mpz_clear(p);
    mpz_clear(randNumb);
    mpz_clear(millerrabin);
    mpz_clear(mrbase);
    mpz_clear(nmenos1);


    free(n);

    gmp_randclear(state);


    return 0;
}

int floorlog(int num, int *initial)
{
    int i = 0;
    while(num > 31)
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

void primeGenPocklington2(int nbits, mpz_t n, gmp_randstate_t state, mpz_t k, mpz_t p, mpz_t r, mpz_t randNumb, mpz_t millerrabin, mpz_t mrbase, mpz_t nmenos1, mpz_t base, mpz_t criterion, mpz_t mcd, mpz_t filter)
{
    //Declaración de variables
    // n es el candidato
    // n = pr + 1
    // r = 2k
    // El resto de variables se utiliza dentro de los tests de primalidad

    int i, aux, exp, m, proof = 0, bsize, initial[1], j, numtests, initialStart;
    int phi[4] = {2, 3, 5, 7};
    unsigned long mrfactor;

            aux = floorlog(nbits, initial);
            initialStart = initial[0] + 2;

            proof = 0;

            //Generación de un número primo aleatorio de 32 bits 
            //Empezar con candidato de 2^31

            randomNBitOddNumber(p, initialStart, state);

            //El primer primo se demuestra con Miller-Rabin
            //El ciclo termina cuando se encuentra un primo demostrado

            while(!proof)
            {
                mpz_gcd(mcd, p, filter);
                if(mpz_cmp_ui(mcd, 1) == 0)
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
                else
                    randomNBitOddNumber(p, initialStart, state);
            }

            bsize = mpz_sizeinbase(p, 2);

            for(i = 1; i < aux; i++)
            {
                do
                {
                    do
                    {
                        do
                        {
                            do
                            {
                                mpz_rrandomb(k, state, bsize-1);
                            }while(mpz_cmp(k, p) > 0);

                            mpz_mul_ui(r, k, 2);
                            mpz_mul(n, r, p);
                            mpz_add_ui(n, n, 1);
                            mpz_gcd(mcd, n, filter);
                        }while(mpz_cmp_ui(mcd, 1) != 0);

                    }while(compositefilter(n) == 0);

                }while(!pocklingtonTest(n, p, r, base, criterion, mcd));

                mpz_set(p, n);
                bsize = mpz_sizeinbase(p, 2);
            }

            m = nbits - bsize;

            do
            {
                do
                {
                    do
                    {

                        do
                        {
                            mpz_rrandomb(k, state, m-1);
                        }while(mpz_cmp(k, p) > 0);

                        mpz_mul_ui(r, k, 2);
                        mpz_mul(n, r, p);
                        mpz_add_ui(n, n, 1);
                        mpz_gcd(mcd, n, filter);
                    }while(mpz_cmp_ui(mcd, 1) != 0);

                }while(compositefilter(n) == 0);

        }while(!pocklingtonTest(n, p, r, base, criterion, mcd));
}

int compositefilter(mpz_srcptr n)
{
    unsigned long int ln2;
    
    unsigned long L = 1;
    unsigned long M,next_max;
    long i,j,prev_tested,stop;

     unsigned long primes[982] = {103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,1153,1163,1171,1181,1187,1193,1201,1213,1217,1223,1229,1231,1237,1249,1259,1277,1279,1283,1289,1291,1297,1301,1303,1307,1319,1321,1327,1361,1367,1373,1381,1399,1409,1423,1427,1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,1489,1493,1499,1511,1523,1531,1543,1549,1553,1559,1567,1571,1579,1583,1597,1601,1607,1609,1613,1619,1621,1627,1637,1657,1663,1667,1669,1693,1697,1699,1709,1721,1723,1733,1741,1747,1753,1759,1777,1783,1787,1789,1801,1811,1823,1831,1847,1861,1867,1871,1873,1877,1879,1889,1901,1907,1913,1931,1933,1949,1951,1973,1979,1987,1993,1997,1999,2003,2011,2017,2027,2029,2039,2053,2063,2069,2081,2083,2087,2089,2099,2111,2113,2129,2131,2137,2141,2143,2153,2161,2179,2203,2207,2213,2221,2237,2239,2243,2251,2267,2269,2273,2281,2287,2293,2297,2309,2311,2333,2339,2341,2347,2351,2357,2371,2377,2381,2383,2389,2393,2399,2411,2417,2423,2437,2441,2447,2459,2467,2473,2477,2503,2521,2531,2539,2543,2549,2551,2557,2579,2591,2593,2609,2617,2621,2633,2647,2657,2659,2663,2671,2677,2683,2687,2689,2693,2699,2707,2711,2713,2719,2729,2731,2741,2749,2753,2767,2777,2789,2791,2797,2801,2803,2819,2833,2837,2843,2851,2857,2861,2879,2887,2897,2903,2909,2917,2927,2939,2953,2957,2963,2969,2971,2999,3001,3011,3019,3023,3037,3041,3049,3061,3067,3079,3083,3089,3109,3119,3121,3137,3163,3167,3169,3181,3187,3191,3203,3209,3217,3221,3229,3251,3253,3257,3259,3271,3299,3301,3307,3313,3319,3323,3329,3331,3343,3347,3359,3361,3371,3373,3389,3391,3407,3413,3433,3449,3457,3461,3463,3467,3469,3491,3499,3511,3517,3527,3529,3533,3539,3541,3547,3557,3559,3571,3581,3583,3593,3607,3613,3617,3623,3631,3637,3643,3659,3671,3673,3677,3691,3697,3701,3709,3719,3727,3733,3739,3761,3767,3769,3779,3793,3797,3803,3821,3823,3833,3847,3851,3853,3863,3877,3881,3889,3907,3911,3917,3919,3923,3929,3931,3943,3947,3967,3989,4001,4003,4007,4013,4019,4021,4027,4049,4051,4057,4073,4079,4091,4093,4099,4111,4127,4129,4133,4139,4153,4157,4159,4177,4201,4211,4217,4219,4229,4231,4241,4243,4253,4259,4261,4271,4273,4283,4289,4297,4327,4337,4339,4349,4357,4363,4373,4391,4397,4409,4421,4423,4441,4447,4451,4457,4463,4481,4483,4493,4507,4513,4517,4519,4523,4547,4549,4561,4567,4583,4591,4597,4603,4621,4637,4639,4643,4649,4651,4657,4663,4673,4679,4691,4703,4721,4723,4729,4733,4751,4759,4783,4787,4789,4793,4799,4801,4813,4817,4831,4861,4871,4877,4889,4903,4909,4919,4931,4933,4937,4943,4951,4957,4967,4969,4973,4987,4993,4999,5003,5009,5011,5021,5023,5039,5051,5059,5077,5081,5087,5099,5101,5107,5113,5119,5147,5153,5167,5171,5179,5189,5197,5209,5227,5231,5233,5237,5261,5273,5279,5281,5297,5303,5309,5323,5333,5347,5351,5381,5387,5393,5399,5407,5413,5417,5419,5431,5437,5441,5443,5449,5471,5477,5479,5483,5501,5503,5507,5519,5521,5527,5531,5557,5563,5569,5573,5581,5591,5623,5639,5641,5647,5651,5653,5657,5659,5669,5683,5689,5693,5701,5711,5717,5737,5741,5743,5749,5779,5783,5791,5801,5807,5813,5821,5827,5839,5843,5849,5851,5857,5861,5867,5869,5879,5881,5897,5903,5923,5927,5939,5953,5981,5987,6007,6011,6029,6037,6043,6047,6053,6067,6073,6079,6089,6091,6101,6113,6121,6131,6133,6143,6151,6163,6173,6197,6199,6203,6211,6217,6221,6229,6247,6257,6263,6269,6271,6277,6287,6299,6301,6311,6317,6323,6329,6337,6343,6353,6359,6361,6367,6373,6379,6389,6397,6421,6427,6449,6451,6469,6473,6481,6491,6521,6529,6547,6551,6553,6563,6569,6571,6577,6581,6599,6607,6619,6637,6653,6659,6661,6673,6679,6689,6691,6701,6703,6709,6719,6733,6737,6761,6763,6779,6781,6791,6793,6803,6823,6827,6829,6833,6841,6857,6863,6869,6871,6883,6899,6907,6911,6917,6947,6949,6959,6961,6967,6971,6977,6983,6991,6997,7001,7013,7019,7027,7039,7043,7057,7069,7079,7103,7109,7121,7127,7129,7151,7159,7177,7187,7193,7207,7211,7213,7219,7229,7237,7243,7247,7253,7283,7297,7307,7309,7321,7331,7333,7349,7351,7369,7393,7411,7417,7433,7451,7457,7459,7477,7481,7487,7489,7499,7507,7517,7523,7529,7537,7541,7547,7549,7559,7561,7573,7577,7583,7589,7591,7603,7607,7621,7639,7643,7649,7669,7673,7681,7687,7691,7699,7703,7717,7723,7727,7741,7753,7757,7759,7789,7793,7817,7823,7829,7841,7853,7867,7873,7877,7879,7883,7901,7907,7919,7927,7933,7937,7949,7951,7963,7993, 8009};

    ln2 = mpz_sizeinbase(n, 2);

    prev_tested = 0;
    stop = 982;
    for ( i = 0 ; i < 982; i++ )
    {
        if ( primes[i] < ln2 )
        {
            if ( primes[i] < (18446744073709551615/L) )
            {
                L *= primes[i];
            }
            else
            {
                M = mpz_fdiv_ui(n,L);
                for ( j = prev_tested ; j < i ; j++ )
                {
                    if ( ( M % (primes[j]) ) == 0 )
                    {
                        return 0;
                    }
                }
                prev_tested = i;
                L = primes[i];
            }
        }
        else
        {
            stop = i;
            break;
        }
    }
    M = mpz_fdiv_ui(n,L);
    for ( j = prev_tested ; j < stop ; j++ )
    {
        if ( ( M % (primes[j]) ) == 0 )
        {
            return 0;
        }
    }
    return 1;
}
