#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gmp.h>

// ----- Cria crivo simples até sqrt(fim) -----
void sieve_base(unsigned long limit, bool *base)
{
    for (unsigned long i = 0; i <= limit; i++)
        base[i] = true;

    base[0] = base[1] = false;

    for (unsigned long i = 2; i * i <= limit; i++)
        if (base[i])
            for (unsigned long j = i * i; j <= limit; j += i)
                base[j] = false;
}

// ----- Conta primos no intervalo gigantesco [inicio, fim] -----
unsigned long count_primes_big_range(const mpz_t inicio, const mpz_t fim)
{
    mpz_t one, temp, start, idx;
    mpz_inits(one, temp, start, idx, NULL);
    mpz_set_ui(one, 1);

    // Se fim < 2 → zero primos
    if (mpz_cmp_ui(fim, 2) < 0)
        return 0;

    // Ajusta início
    mpz_t A;
    mpz_init_set(A, inicio);
    if (mpz_cmp_ui(A, 2) < 0)
        mpz_set_ui(A, 2);

    // sqrt(fim)
    unsigned long limit = mpz_get_ui(fim);
    limit = (unsigned long)mpz_sqrt_ui(fim);

    // Crivo base até sqrt(fim)
    bool *base = malloc((limit + 1) * sizeof(bool));
    sieve_base(limit, base);

    // Quantidade de números no intervalo
    mpz_t range;
    mpz_init(range);
    mpz_sub(range, fim, A);
    mpz_add_ui(range, range, 1);

    size_t size = mpz_get_ui(range);

    bool *prime = malloc(size * sizeof(bool));
    for (size_t i = 0; i < size; i++)
        prime[i] = true;

    // --- Marca múltiplos para cada primo base ---
    for (unsigned long p = 2; p <= limit; p++)
    {

        if (!base[p])
            continue;

        // start = ceil(A / p) * p
        mpz_set_ui(temp, p);

        mpz_fdiv_q(start, A, temp);  // start = A/p
        mpz_mul(start, start, temp); // start = floor(A/p)*p

        if (mpz_cmp(start, A) < 0)
            mpz_add(start, start, temp);

        // Evita marcar o próprio primo no intervalo
        if (mpz_cmp_ui(start, p) == 0)
            mpz_add(start, start, temp);

        // Marca múltiplos: prime[(j - A)]
        mpz_t j;
        mpz_init_set(j, start);

        while (mpz_cmp(j, fim) <= 0)
        {
            mpz_sub(idx, j, A);
            size_t index = mpz_get_ui(idx);
            prime[index] = false;
            mpz_add(j, j, temp);
        }

        mpz_clear(j);
    }

    // Conta primos
    unsigned long count = 0;
    for (size_t i = 0; i < size; i++)
        if (prime[i])
            count++;

    // Libera memória
    free(base);
    free(prime);

    mpz_clears(one, temp, start, idx, A, range, NULL);

    return count;
}

int main()
{
    mpz_t inicio, fim;
    mpz_init(inicio);
    mpz_init(fim);

    // Intervalo gigantesco
    mpz_set_str(inicio, "1000000000000000000000000", 10); // 10^24
    mpz_set_str(fim, "1000000000000000000100000", 10);    // 10^24 + 10^5

    unsigned long result = count_primes_big_range(inicio, fim);

    gmp_printf("Quantidade de primos no intervalo [%Zd, %Zd] = %lu\n",
               inicio, fim, result);

    mpz_clear(inicio);
    mpz_clear(fim);

    return 0;
}
