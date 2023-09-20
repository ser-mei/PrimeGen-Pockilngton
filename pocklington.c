#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <time.h>

//Floor function
int floorlog(int num);

//Pocklington test
int pocklingtonTest(mpz_t n, mpz_t p, mpz_t r);

//Bit count
void bitcount(mpz_t n);

//Random n bit number
void randomNBitNumber(mpz_t num, int nbits);

int main()
{
    //Declaración de variables
    mpz_t k, n, p, r, randNumb;

    //Declaración de estado para rng
    gmp_randstate_t state;

    //Variables de tiempo
    clock_t start, end;

    //Variables enteras
    int i, nbits, aux, exp, m;

    //Inicialización de variables de gmp
    mpz_init(k);
    mpz_init(n);
    mpz_init(p);
    mpz_init(r);
    mpz_init(randNumb);

    //Inicialización de estado para rng
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, time(NULL));

    //Input de número de bits
    printf("Enter the number of bits: ");
    scanf("%d", &nbits);

    //Inicio de medición del tiempo
    start = clock();

    //Generación de un número primo aleatorio de 32 bits 
    //randomNBitNumber(p, 32);
    mpz_set_ui(p, 4294967295);
    mpz_urandomb(randNumb, state, 31);
    mpz_sub(p, p, randNumb);
    //El primer primo debe ser primo demostrado
    //Se puede partir con un gcd con 105
    //luego se puede hacer con Miller-Rabin con 4 bases phi_4

    while(mpz_probab_prime_p(p, 15) == 0)
    {
        mpz_urandomb(randNumb, state, 31);
        while(mpz_even_p(randNumb) == 0)
            mpz_urandomb(randNumb, state, 31);
            //randomNBitNumber(p, 32);
        mpz_set_ui(p, 4294967295);
        mpz_sub(p, p, randNumb);
    }

    //Agregar medición de tiempo
    
    gmp_printf("Prime p = %Zd\n", p);

    //Ciclo para generar primos más grandes
    aux = floorlog(nbits);
    for(i = 5; i < aux; i++)
    {
        exp = 1 << i;
        //printf("exponente = 2^%d = %d\n", i, exp);
        //mpz_ui_pow_ui(exp, 2, i);
        //mpz_sub_ui(exp, exp, 1);
        mpz_urandomb(k, state, exp);
        //mpz_urandomb(k, state, i);

        mpz_mul_ui(r, k, 2);
        mpz_mul(n, r, p);
        mpz_add_ui(n, n, 1);
        gmp_printf("Candidate n = %Zd = %Zd * %Zd + 1\n", n, p, r);

        pocklingtonTest(n, p, r);

        while(mpz_probab_prime_p(n, 15) == 0)
        {
            mpz_urandomb(k, state, exp);
            mpz_mul_ui(r, k, 2);
            mpz_mul(n, r, p);
            mpz_add_ui(n, n, 1);
            pocklingtonTest(n, p, r);
        }
        mpz_set(p, n);
    }
    //if(mpz_probab_prime_p(n, 15) > 0)
    //    gmp_printf(" Probably Prime n = %Zd\n", n);
    //else
    //    gmp_printf("Composite n = %Zd\n", n);

    bitcount(n);

    m = 1 << aux;
    printf("logfloor de nbits = %d\n", aux);
    printf("m = %d\n", m);
    m = nbits - m;
    printf("m = %d\n", m);


    //mpz_ui_pow_ui(m, 2, aux);
    //mpz_sub(m, r, m);
    //mpz_sub_ui(m, m ,1);

    mpz_urandomb(k, state, m);
    mpz_mul_ui(r, k, 2);
    mpz_mul(n, r, p);
    mpz_add_ui(n, n, 1);

    pocklingtonTest(n, p, r);
    while(mpz_probab_prime_p(n, 15) == 0)
    {
        mpz_urandomb(k, state, m);
        mpz_mul_ui(r, k, 2);
        mpz_mul(n, r, p);
        mpz_add_ui(n, n, 1);
        pocklingtonTest(n, p, r);
    }
    mpz_set(p, n);

    gmp_printf("Número primo n = %Zd \n", p);

    end = clock();

    if(mpz_probab_prime_p(p, 15) > 0)
        gmp_printf(" Probably Prime n = %Zd\n", p);
    else
        gmp_printf("Composite n = %Zd\n", p);


    bitcount(p);

    printf("Tiempo de búsqueda para número primo de %d bits: %f\n", nbits, ((double)end - start) / CLOCKS_PER_SEC);


    //Liberación de memoria
    mpz_clear(k);
    mpz_clear(n);
    mpz_clear(p);
    mpz_clear(r);
    mpz_clear(randNumb);

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

int pocklingtonTest(mpz_t n, mpz_t p, mpz_t r)
{
    mpz_t base, fermat, nmenos1, mcd;
    mpz_init(base);
    mpz_init(fermat);
    mpz_init(nmenos1);
    mpz_init(mcd);

    mpz_set(mcd, n);
    mpz_set_ui(base, 2);
    mpz_sub_ui(nmenos1, n, 1);
    mpz_powm(fermat, base, nmenos1, n);

    //Hacer el gcd con 105

    if(mpz_cmp_ui(fermat, 1) == 0)
    {
        printf("Fermat test passed\n");
        mpz_powm(mcd, base, r, n);
        mpz_sub_ui(mcd, mcd, 1);
        mpz_gcd(mcd, mcd, n);
        if(mpz_cmp_ui(mcd, 1) == 0)
        {
            printf("Pocklington criterion passed\n");
            return 1;
        }
        else
        {
            printf("Pocklington test failed\n");
            return 0;
        }
    }
    else
        return 0;

    mpz_clear(base);
    mpz_clear(fermat);
    mpz_clear(nmenos1);
    mpz_clear(mcd);
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

void randomNBitNumber(mpz_t num, int nbits)
{
    mpz_t randNumb;
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, time(NULL));

    mpz_init(randNumb);

    //mpz_set_ui(num, 4294967295);
    mpz_ui_pow_ui(num, 2, nbits);
    mpz_sub_ui(num, num, 1);
    //gmp_printf("num = %Zd es 4294967295?\n", num);
    mpz_urandomb(randNumb, state, nbits-1);
    mpz_sub(num, num, randNumb);
    //gmp_printf("num = %Zd\n", num);

    mpz_clear(randNumb);

}