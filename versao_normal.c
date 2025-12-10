#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include "mpi.h"

// ------------------ CRIVO BASE ------------------------
void sieve_base(int limit, bool base[])
{
    for (int i = 0; i <= limit; i++)
        base[i] = true;
    base[0] = base[1] = false;

    for (int i = 2; i * i <= limit; i++)
        if (base[i])
            for (int j = i * i; j <= limit; j += i)
                base[j] = false;
}

// ------------------ CRIVO SEGMENTADO -------------------
unsigned int count_primes_range(unsigned long long inicio,
                                unsigned long long fim)
{
    if (fim < 2)
        return 0;
    if (inicio < 2)
        inicio = 2;

    int limit = floor(sqrt(fim));
    bool *base = malloc((limit + 1) * sizeof(bool));
    sieve_base(limit, base);

    unsigned long long size = fim - inicio + 1;
    bool *prime = malloc(size * sizeof(bool));
    for (unsigned long long i = 0; i < size; i++)
        prime[i] = true;

    for (int p = 2; p <= limit; p++)
    {
        if (base[p])
        {
            long long start = (inicio / p) * p;
            if (start < inicio)
                start += p;
            if (start == p)
                start += p;

            for (long long j = start; j <= fim; j += p)
                prime[j - inicio] = false;
        }
    }

    unsigned int count = 0;
    for (unsigned long long i = 0; i < size; i++)
        if (prime[i])
            count++;

    free(base);
    free(prime);
    return count;
}

int main(int argc, char *argv[])
{
    unsigned long long inicio_global = 1;
    unsigned long long fim_global = 100000000ULL; // 1e8

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Barrier(MPI_COMM_WORLD);
    double t_paralelo_ini = MPI_Wtime();

    // Subintervalos ...
    unsigned long long inicio_local, fim_local;
    unsigned long long bloco = (fim_global - inicio_global + 1) / 4;

    if (rank == 0)
    {
        inicio_local = inicio_global;
        fim_local = inicio_global + bloco - 1;
    }
    else if (rank == 1)
    {
        inicio_local = inicio_global + bloco;
        fim_local = inicio_global + 2 * bloco - 1;
    }
    else if (rank == 2)
    {
        inicio_local = inicio_global + 2 * bloco;
        fim_local = inicio_global + 3 * bloco - 1;
    }
    else
    {
        inicio_local = inicio_global + 3 * bloco;
        fim_local = fim_global;
    }

    printf("Rank %d executando intervalo [%llu, %llu]\n",
           rank, inicio_local, fim_local);

    unsigned int local_count = count_primes_range(inicio_local, fim_local);
    
    // Resultado final
    unsigned int total_count = 0;
    MPI_Reduce(&local_count, &total_count, 1,
               MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    double t_paralelo_fim = MPI_Wtime();

    if (rank == 0)
    {
        double tempo_paralelo = t_paralelo_fim - t_paralelo_ini;

        printf("\n=== RESULTADO PARALELO ===\n");
        printf("Total de primos: %u\n", total_count);
        printf("Tempo paralelo: %.6f segundos\n\n", tempo_paralelo);

        printf("Executando versão sequencial...\n");

        double t_seq_ini = MPI_Wtime();
        unsigned int seq = count_primes_range(inicio_global, fim_global);
        double t_seq_fim = MPI_Wtime();

        double tempo_sequencial = t_seq_fim - t_seq_ini;

        printf("\n=== RESULTADO SEQUENCIAL ===\n");
        printf("Total de primos: %u\n", seq);
        printf("Tempo sequencial: %.6f segundos\n\n", tempo_sequencial);

        double speedup = tempo_sequencial / tempo_paralelo;
        double eficiencia = speedup / size;

        printf("=== MÉTRICAS ===\n");
        printf("Speedup: %.4f\n", speedup);
        printf("Eficiência: %.4f (%.2f%%)\n",
               eficiencia, eficiencia * 100.0);
    }

    MPI_Finalize();
    return 0;
}
