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

struct WorkerArg
{
  int id;
  int start;
  int end;
  int fd;
  int attractors[STATE_COUNT];
  int cyclic_counts;
};

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
  printf("    State: ");
  for (int i = 0; i < GENE_COUNT; i++)
  {
    printf("%d", (state >> i) & 1);
  }
  printf("\n");
}

void helper_print_results(struct WorkerArg *arg)
{
  int total = 0;
  for (int i = 0; i < STATE_COUNT; i++)
  {
    if (arg->attractors[i] > 0)
    {
      total += arg->attractors[i];
      printf("  Attractor %3d: %9d\n", i, arg->attractors[i]);
      helper_print_state(i);
    }
  }
  printf("Cyclic: %d\n", arg->cyclic_counts);
  printf("Total: %d\n", total + arg->cyclic_counts);
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

static inline int gene_update_value(int id, short state)
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

void *worker(struct WorkerArg *arg)
{
  printf("Thread %3d: %7d ~%7d\n", arg->id, arg->start, arg->end);
  int previous_states[STATE_COUNT];
  int previous_states_len;

  for (int i = arg->start; i < arg->end; i++)
  {
    int line[GENE_COUNT];
    if (!~pread(arg->fd, line, sizeof(line), i * sizeof(line)))
    {
      perror("pread");
      exit(1);
    }

    for (int initial_state = 0; initial_state < STATE_COUNT; initial_state++)
    {
      previous_states_len = 0;
      int state = initial_state;

      for (;;)
      {
        previous_states[previous_states_len++] = state;

        for (int time = 0; time < GENE_COUNT; time++)
        {
          int staging_state = state;

          for (int gene = 0; gene < GENE_COUNT; gene++)
          {
            if (line[gene] == time)
            {
              int update_value = gene_update_value(gene, staging_state) != 0;
              staging_state = (staging_state & ~(1 << gene)) | (update_value << gene);
            }
          }
          state = staging_state;
        }

        // new test
        // if (previous_states_len > 1 && previous_states[previous_states_len - 1] == state)
        // {
        //   arg->attractors[state]++;
        //   break;
        // }

        for (int check_cyclic = 0; check_cyclic < previous_states_len; check_cyclic++)
        {
          if (previous_states[check_cyclic] == state)
          {
            if (check_cyclic + 1 == previous_states_len)
              arg->attractors[state]++;
            else
              arg->cyclic_counts++;
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
  int fd = open("groups.bin", O_RDONLY, 0);
  off_t file_size = lseek(fd, 0, SEEK_END);
  lseek(fd, 0, SEEK_SET);
  int TOTAL_LINES = file_size / (4 * GENE_COUNT);
  printf("Total lines: %d\n", TOTAL_LINES);

  struct WorkerArg args[THREAD_COUNT];
  pthread_t threads[THREAD_COUNT];
  int chunk_size = TOTAL_LINES / THREAD_COUNT;

  for (int i = 0; i < THREAD_COUNT; i++)
  {
    args[i].id = i;
    args[i].start = i * chunk_size;
    args[i].end = (i == THREAD_COUNT - 1) ? TOTAL_LINES : (i + 1) * chunk_size;
    args[i].fd = fd;
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

  close(fd);
  return 0;
}
