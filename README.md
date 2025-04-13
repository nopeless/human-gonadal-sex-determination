# So, what happened

Well, I was talking to one of the researchers and they said it took 2 hours to run a certain section of the program and I said I can do it under a minute so heres that

Also somehow I couldn't get the exact original results but the ratio is close enough so I guess its fine

# More explanation

the state is stored in a binary format in an integer

due to how the original paper works, GATA4 is the **least** significant bit, but binaries are printed **backwards**

```c
// internal
int state = 0b10000000011;
//            ^--------++ CTNNB1
//                     ^+ NR5A1
//                      ^ GATA4

// printed
// 1100000001

// GENE 0 is GATA4, also appears first on the binary table
```

# How to run

Prepare precomputed binary file

```sh
$ python preprocess_csv.py
```

Run the program (edit some macros as needed)

```less
```
