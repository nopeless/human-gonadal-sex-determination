# Paper Reproduction and Technical Error Analysis: "Dynamical Robustness of a Boolean Model for Human Gonadal Sex Determination"

## Incident

After a casual conversation with Erika Vivanco about her work related to the SRY gene she had introduced me to her paper *Dynamical Robustness of a Boolean Model for Human Gonadal Sex Determination (Vivanco 2024)*. She had explained that the R code used for the production of results in the paper took her 2 hours to run through the required dataset. This conversation sparked my curiosity and I have made a bet with her that I could write a more efficient implementation of the code in C and cut the required time to less than a minute. Without much experience in Computational Systems Biology, I have decided to take this as an opportunity to learn more about the field and explore gene expressions and the interplay with computer modeling systems, also hoping some validation from a researcher in the field as well.

## Technical

The nature of Boolean Networks have been summarized extrodinarily well in section 1.2 (page 2~3) of the paper. Both synchronous and block sequential update schedules have been implemented by implementing a block sequential update schedule that includes a synchronous update schedule as the first case (trivially also present in `matlab.csv:1`). The problem space thus becomes finding a way to evaluate `2**10` possible state configurations over $895029$ possible update schedules (subset), totalling to $916509696$ configurations. The exact way this scales over `N` genes is additionally Fubini Numberes (ordered Bell numbers) times 2 to the power of N. However, we work with a much smaller set (0.875%) of the total configurations as many configurations eventually result in duplicate states.

$F_n=\sum_{k=1}^{n} k! \, S(n,k)$

$S(n,k)=\frac{1}{k!}\sum_{j=0}^{k}(-1)^{k-j}
\binom{k}{j}j^n$

The asymptotic formula is $F_n \approx \frac{n!}{2(\ln 2)^{n+1}}$ with upper bound of $O(n^n)$

Thus, manually simulating every possible configuration becomes infeasible quickly; The paper acknowledged this and performed an exhaustive search by first reducing the original 19-node boolean network into a 10-node boolean network. The paper additionally notes that the reduced boolean network produces the same fixed points faithfully.

Assuming duplicate group state configuration is $1-1/n^2$ and original R code had been used, an 11-node boolean network would have taken $2 \times \frac{1622632573}{102247563} \times \frac{10^2}{11^2} \approx 26$ (*constants are* $F_{11}/F_{10}$) hours to evaluate, with 12 nodes taking about 2 weeks.

## Implementation, optimizations and modelling

The optimization strategy taking in this project was to parallelize the evaluation of the boolean network and utilize cached lookup to avoid computation of genes entirely. The relatively small state size ($2^{10}$) allows a precomputed lookup table to be used to evaluate the next state of a gene given its current state.

For precomputing gene interactions, `Network10.txt` provided by the researcher was parsed via a python program to generate a C integer binary expression to be pasted into the C source code. Additionally the group state configurations that were precomputed and provided in `matlab.csv` was also parsed by a python program to generate a C friendly binary file `groups.bin` to be read by worker threads. Reducing the file size and removing parsing overhead allowed rapid iteration of the source code by making the initializtion step single `mmap` call.

Additionally, the state is stored in a binary format in an integer just as the original paper, but in this implementation GATA4 is the **least** significant bit and binaries are printed **backwards**

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

Certain micro optimizations were utilized in hotloops to reduce the number of branches and improve instruction pipelining. Every potential heap allocation was moved to the stack. `gcc` with assembly dump was used to identify potential sites for cache misses and branch misprediction.

As a result, in a machine with a 4GHz CPU with 16 threads, the evaluation of the boolean network completes in 5 seconds, a $\approx 1400$ times speed up.

```less
gcc -Wall -fdiagnostics-color=always main.c -O3 && time sudo nice -n -20 ./a.out
Total lines: 895029
Thread   0:       0 ~  55939
Thread   1:   55939 ~ 111878
Thread   2:  111878 ~ 167817
Thread   3:  167817 ~ 223756
Thread   4:  223756 ~ 279695
Thread   6:  335634 ~ 391573
Thread   5:  279695 ~ 335634
Thread   7:  391573 ~ 447512
Thread   8:  447512 ~ 503451
Thread   9:  503451 ~ 559390
Thread  10:  559390 ~ 615329
Thread  11:  615329 ~ 671268
Thread  12:  671268 ~ 727207
Thread  13:  727207 ~ 783146
Thread  14:  783146 ~ 839085
Thread  15:  839085 ~ 895029
Thread   0 done
Thread   1 done
Thread   2 done
Thread   4 done
Thread   5 done
Thread   3 done
Thread   9 done
Thread   6 done
Thread  10 done
Thread   8 done
Thread   7 done
Thread  15 done
Thread  13 done
Thread  14 done
Thread  11 done
Thread  12 done
  Attractor   0:  53329046 (5.82%)
    State: 0000000000
  Attractor 251: 514512013 (56.18%)
    State: 1101111100
  Attractor 773: 348027806 (38.00%)
    State: 1010000011
Cyclic: 640831 (0.07%)
Total: 916509696
done

real    0m5.306s
user    0m0.001s
sys     0m0.006s
```

## Verification

The result of the evaluation was initially done by comparing a screenshot of an excel sheet of the original results for specific numbers. This heavily aided in the early developing and debugging stages of this project. However, a single number proved to be insufficient to track down specific bugs in the implementation and thus I have requested the researcher for additional artifacts that may assist in the verification of the implementation. One part of her implementation particularly was a red herring as she simply detected whether there exists a known attractor after 100 updates. The implementation in this project had a strict check for every new state configuration. The result was never affected by the implementation because the max length of a state chain was 8 and max length of a cycle was 2 in every configuration.

## Discovery of the Error

While debugging the mismatch between percentages of final attractors and cyclic attractors, I have asked her for a complete evaluation trace for scheme 2 (update cycle of `1111111112`). While tracing the evaluation step, I noticed that my implementation reaches the same final state as the original implementation, but by traversing a different path. After troubleshooting for hours, I have started to read the original R code to understand the logic behind the implementation. During this process, I have discovered that the original implementation had a bug in the evaluation of the boolean network. Specifically, `condition5.2` had an incorrect application of de-morgan's law

from `Network10.txt:6`: `SOX9,((SOX9 & FGF9) | SRY | PGD2 | (GATA4 & NR5A1 & SRY)) & !(WNT4 | CTNNB1)`

```R
    condition5.2 = ((final_bin[,5,,j] == 1 & final_bin[,6,,j] == 1) |
                      (final_bin[,1,,j] == 1 & final_bin[,2,,j] == 1 & final_bin[,4,,j] == 1) |
                      final_bin[,4,,j] == 1 |
                      final_bin[,7,,j] == 1) &
      (final_bin[,9,,j] == 0 | final_bin[,10,,j] == 0)
```

Resulting in an incorrect gene expression for `SOX9` based on `WNT4` and `CTNNB1`. This was demonstrated with a simple written truth table evaluator from [one of my other projects](https://github.com/nopeless/python-generic-shunting-yard-parser-evaluator) during the discussion

```txt
 CTNNB1  WNT4  !(WNT4 | CTNNB1)
   True  True             False
   True False             False
  False  True             False
  False False              True
```

Additionally, there is a similar error from `condition10.2`, resulting in an incorrect gene expression for `CTNNB1`

from `Network10.txt:11`: `CTNNB1,(WNT4 | CTNNB1) & !(SRY | (SOX9 & GATA4 & NR5A1) | DKK1)`

The above can be written as

$\text{CTNNB1} = (\text{WNT4} \lor \text{CTNNB1}) \land \neg(\text{SRY} \lor (\text{SOX9} \land \text{GATA4} \land \text{NR5A1}) \lor \text{DKK1})$

The R code

```R
    condition10.2 = (final_bin[,9,,j] == 1 | final_bin[,10,,j] == 1) &
      (final_bin[,4,,j] == 0 |
         (final_bin[,1,,j] == 0 & final_bin[,2,,j] == 0 & final_bin[,5,,j] == 0) |
         final_bin[,8,,j] == 0)
```

The above can be written as

$\text{CTNNB1} = (\text{WNT4} \lor \text{CTNNB1}) \land \left[\neg \text{SRY} \lor(\neg \text{GATA4} \land \neg \text{NR5A1} \land \neg \text{SOX9}) \lor \neg \text{DKK1} \right]$


This resulted in incorrect expression for `CTNNB1` gene with regards to `SRY`, `SOX9`, `GATA4`, `NR5A1`, and `DKK1` and necessary inhibitory genes did not get accounted for.

**Despite the bug**, the overall results didn't change. The yield was 5.82%, 38.00%, 56.18% for the 3 attractors stated in the paper where the incorrect yield was 1.85%, 43.85%, 54.30% (the attractor basin had no change at all). Critically, the results with the fixed bug do not change the findings of the paper, since the conclusions of the paper are aimed towards robustness and not exact values

## Limitation of the implementation

There are so many avenues for improving this design, mostly due to the fact that while the number of final states grows less than linearly, the number of initial group states grows superexponentially. A much better approach to create the fastest and scalable boolean network evaluator would be to work backwards from the final states and find all possible initial states that lead to the final states. This would result in many initial states being presented as a combination of segments for each final state, easily compressible and representible for future computer processing.

## Conclusion

The original paper started with a straightforward explanation of Boolean Networks. Despite having no prior experience in this field, I searched online to understand the underlying material and eventually produce a faithful implementation. In the end, not only was I able to achieve a x1,400 speed up, I also found 2 implementation errors in the original R code. This implementation of a boolean model could potentially aid research in non-CS fields where scientists that are not accustomed to heavy optimization, as many other implementations of boolean networks are stochastic and asynchronous.

I am so delighted from going beyond simply reproducing the same scientific results. Understanding the model, debugging its behavior, and examining discrepancies with the head researcher was a unique experience I'll never forget. I am grateful for the opportunity to explore this field and have discussions with Erika, who made this project possible; It sparked my interest in computational systems biology and I hope to have similar research opportunities in the future.

## References

*Vivanco E, Goles E, Montalva-Medel M, Poupin MJ. Dynamical robustness of a Boolean model for the human gonadal sex determination. Comput Biol Chem. 2024;113:108225. doi:10.1016/j.compbiolchem.2024.108225*

---

*Note: Originally this repository had jokes and a rather sloppy explanation of my work. I have decided to clean up the README to present it as an evidence of my enthusiasm and potential for research. I am a joyful and casual person at heart and like to have fun, and this small research project is one of my thousands of endeavors to learn more about things around me. You can browse the previous version of this read me by looking at the commits. If you have stumbled upon this repository as my friend, please don't be alienated by the new tone :\)*
