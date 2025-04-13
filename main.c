#include <sys/mman.h>
#include <fcntl.h>
#include <stdio.h>
#include <pthread.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#define GENE_COUNT 10
#define STATE_COUNT (1 << GENE_COUNT)
#define THREAD_COUNT 16

// #define THREAD_COUNT 1
// #define DEBUG_INDIVIDUAL_STATES
// #define DEBUG_INDIVIDUAL_STATES_ALL

struct WorkerArg
{
  int id;
  int start;
  int end;
  int (*data)[GENE_COUNT];
  int cyclic_counts;
  // Helps with big array alignment
  __attribute__((aligned(64))) int attractors[STATE_COUNT];
  // Add padding to avoid cache sharing between numa nodes
  __attribute__((aligned(64))) int _padding[8];
};

// Global lookup table
int gene_lookup_table[STATE_COUNT][GENE_COUNT];

void helper_print_line(int *line)
{
  printf("  Line: ");
  for (int i = 0; i < GENE_COUNT; i++)
  {
    printf("%d ", line[i]);
  }
  printf("\n");
}
void helper_print_state(int state)
{
  for (int i = 0; i < GENE_COUNT; i++)
  {
    printf("%d", (state >> i) & 1);
  }
}

void helper_print_results(struct WorkerArg *arg)
{
  int total = 0;

  for (int i = 0; i < STATE_COUNT; i++)
  {
    total += arg->attractors[i];
  }

  for (int i = 0; i < STATE_COUNT; i++)
  {
    if (arg->attractors[i] > 0)
    {
      printf("  Attractor %3d: %9d (%.2f%%)\n", i, arg->attractors[i], ((float)arg->attractors[i] * 100 / (float)total));
      printf("    State: ");
      helper_print_state(i);
      printf("\n");
    }
  }
  printf("Cyclic: %d (%.2f%%)\n", arg->cyclic_counts, ((float)arg->cyclic_counts * 100 / (float)total));
  printf("Total: %d\n", total + arg->cyclic_counts);
}

void helper_print_states(int *states, int len)
{
  for (int i = 0; i < len; i++)
  {
    printf("%6d: ", i);
    helper_print_state(states[i]);
    printf("\n");
  }
}

void aggregate_args(struct WorkerArg *args, int num_threads, struct WorkerArg *result)
{
  memset(result->attractors, 0, sizeof(result->attractors));
  result->cyclic_counts = 0;

  for (int i = 0; i < num_threads; i++)
  {
    for (int j = 0; j < STATE_COUNT; j++)
    {
      result->attractors[j] += args[i].attractors[j];
    }
    result->cyclic_counts += args[i].cyclic_counts;
  }
}

static inline int gene_update_value(int id, int state)
{
  switch (id)
  {
  case 0:
    return (state >> 8 & 1) | (state >> 1 & 1) | (state >> 3 & 1);
  case 1:
    return (state >> 0 & 1) & !((state >> 2 & 1) & (state >> 8 & 1));
  case 2:
    return ((state >> 0 & 1) | ((state >> 8 & 1) & (state >> 9 & 1))) & !((state >> 1 & 1) & (state >> 4 & 1));
  case 3:
    return ((state >> 4 & 1) | (state >> 3 & 1)) & !(state >> 9 & 1);
  case 4:
    return (((state >> 4 & 1) & (state >> 5 & 1)) | (state >> 3 & 1) | (state >> 6 & 1) | ((state >> 0 & 1) & (state >> 1 & 1) & (state >> 3 & 1))) & !((state >> 8 & 1) | (state >> 9 & 1));
  case 5:
    return (state >> 4 & 1) & !(state >> 8 & 1);
  case 6:
    return (state >> 4 & 1);
  case 7:
    return (state >> 3 & 1) | (state >> 4 & 1);
  case 8:
    return ((state >> 0 & 1) | (state >> 9 & 1) | (state >> 8 & 1) | (state >> 2 & 1)) & !((state >> 5 & 1) | (state >> 7 & 1));
  case 9:
    return ((state >> 8 & 1) | (state >> 9 & 1)) & !((state >> 3 & 1) | ((state >> 4 & 1) & (state >> 0 & 1) & (state >> 1 & 1)) | (state >> 7 & 1));
  }
  return 0;
}

void precompute_gene_lookup_table()
{
  for (int state = 0; state < STATE_COUNT; state++)
  {
    for (int gene = 0; gene < GENE_COUNT; gene++)
    {
      int value = gene_update_value(gene, state);
      gene_lookup_table[state][gene] = value << gene;
    }
  }
}

void *worker(struct WorkerArg *arg)
{
  printf("Thread %3d: %7d ~%7d\n", arg->id, arg->start, arg->end);
  int previous_states[STATE_COUNT];
  int previous_states_len;

  for (int i = arg->start; i < arg->end; i++)
  {
    int *line = arg->data[i];

    for (int initial_state = 0; initial_state < STATE_COUNT; initial_state++)
    {
      previous_states_len = 0;
      int state = initial_state;

      for (;;)
      {
        previous_states[previous_states_len++] = state;

        for (int time = 0; time < GENE_COUNT; time++)
        {
          int next_state = 0;
          // I have absolutely no idea how, but this makes code faster (30% faster)
          int NOP = 1;

          for (int gene = 0; gene < GENE_COUNT; gene++)
          {
            if (line[gene] == time)
            {
              next_state |= gene_lookup_table[state][gene];

              // next_state &= ~(1 << gene);
              // next_state |= gene_lookup_table[state][gene];

              NOP = 0;
            }
            else
            {
              next_state |= (state & (1 << gene));
            }
          }

          state = next_state;

          if (NOP)
          {
            break;
          }
        }

        for (int check_cyclic = 0; check_cyclic < previous_states_len; check_cyclic++)
        {
          if (previous_states[check_cyclic] == state)
          {
            if (check_cyclic + 1 == previous_states_len)
              arg->attractors[state]++;
            else
              arg->cyclic_counts++;

#ifdef DEBUG_INDIVIDUAL_STATES
            printf("%5d | ", initial_state);
            helper_print_state(initial_state);
            printf(" -> ");
            helper_print_state(state);
            printf(" [%5d] *%d", state, check_cyclic);

            if (check_cyclic + 1 < previous_states_len)
            {
              printf(" (cyclic %5d)", previous_states_len - check_cyclic);
            }

#ifdef DEBUG_INDIVIDUAL_STATES_ALL
            printf("\n");
            helper_print_states(previous_states, previous_states_len);
#else
            printf("\n");
#endif

#endif

            goto break_loop;
          }
        }
      }

    break_loop:;
    }
  }

  printf("Thread %3d done\n", arg->id);

  return NULL;
}

int main(int argc, char **argv)
{
  int fd = open("groups.bin", O_RDONLY);
  if (fd < 0)
  {
    perror("open");
    exit(1);
  }

  off_t file_size = lseek(fd, 0, SEEK_END);

  if (file_size < 0)
  {
    perror("lseek");
    close(fd);
    exit(1);
  }

  int TOTAL_LINES = file_size / (sizeof(int) * GENE_COUNT);

  printf("Total lines: %d\n", TOTAL_LINES);

  int(*data)[GENE_COUNT] = mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
  if (data == MAP_FAILED)
  {
    perror("mmap");
    close(fd);
    exit(1);
  }

  precompute_gene_lookup_table();

  struct WorkerArg args[THREAD_COUNT];
  pthread_t threads[THREAD_COUNT];
  int chunk_size = TOTAL_LINES / THREAD_COUNT;

  for (int i = 0; i < THREAD_COUNT; i++)
  {
    args[i].id = i;
    args[i].start = i * chunk_size;
    args[i].end = (i == THREAD_COUNT - 1) ? TOTAL_LINES : (i + 1) * chunk_size;
    args[i].data = data;
    memset(args[i].attractors, 0, sizeof(args[i].attractors));
    args[i].cyclic_counts = 0;
    pthread_create(&threads[i], NULL, (void *(*)(void *))worker, &args[i]);
  }

  for (int i = 0; i < THREAD_COUNT; i++)
  {
    pthread_join(threads[i], NULL);
  }

  struct WorkerArg final_result;
  aggregate_args(args, THREAD_COUNT, &final_result);

  helper_print_results(&final_result);

  printf("done\n");

  munmap(data, file_size);
  close(fd);
  return 0;
}
