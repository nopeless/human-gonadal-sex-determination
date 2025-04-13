import re
import itertools

lines = open("Network10.txt").readlines()[1:]

# print(lines)

genes = [re.match(r"^\w+(?=,)", line).group(0) for line in lines]

gene_table = dict(zip(genes, itertools.count()))


# Use this to index genes
def EXPR_REPLACER(n: int):
    # return f"(state >> {n} & 1)"
    return f"final_bin[,{n+1},,j] == 1"


# Use this to generate functions
def LINE_REPLACER(case, expr: str):
    # return f"case {case}:\nreturn {expr};"
    return f"condition{case+1}.2 = {expr}"


import sys

if len(sys.argv) < 2:
    for case, line in enumerate(lines):
        definition = re.search(r"(?<=,).+", line).group(0)

        # Clean up work for better reading
        definition = re.sub(r"\s+", " ", definition)

        # Use number defs
        definition = re.sub(
            r"\w+", lambda x: EXPR_REPLACER(gene_table[x.group(0)]), definition
        )
        print(LINE_REPLACER(case, definition))

else:
    # NOTE: this is the most shittiest code ive written in my entire life
    # This is purely for debugging purposes and I absolute do not care
    def eval_state(user_input):
        direction = user_input[0]
        if direction not in ["L", "R"]:
            raise ValueError("Direction must be 'L' or 'R'.")

        state = [{"0": False, "1": True}[c] for c in user_input[1:]]

        if len(state) != len(genes):
            raise ValueError(
                f"State length {len(state)} does not match number of genes {len(genes)}."
            )

        if direction == "R":
            state = state[::-1]

        def print_state(state):
            nonlocal direction
            if direction == "R":
                state = state[::-1]

            print(direction + "".join(["1" if x else "0" for x in state]))

        def print_gene_configuration(state):
            for i, gene in enumerate(state):
                print(f"{genes[i]:>6}: {str(gene)}")

        print("=== INITIAL ===")
        print_gene_configuration(state)
        print_state(state)

        seen = {
            tuple(state): 0,
        }

        last_result = state

        for iteration in range(1, 1000):
            print(f"=== Iteration {iteration} ===")

            # Evaluate all genes
            result = []
            for gene, gene_index in gene_table.items():

                eval_string = lines[gene_index].split(",")[1].strip()
                # Substitute all symbols with true false
                eval_string = re.sub(
                    r"\w+",
                    lambda x: f"{last_result[gene_table[x.group(0)]]}",
                    eval_string,
                )
                # substitute all & with and
                eval_string = re.sub(r"&", " and ", eval_string)
                # substitute all | with or
                eval_string = re.sub(r"\|", " or ", eval_string)
                # substitute all ! with not
                eval_string = re.sub(r"!", " not ", eval_string)

                spaces = r"\s+"

                activated = eval(eval_string)

                print(
                    f" => {gene:>6}: {str(activated):<5} <- {re.sub(spaces, ' ', eval_string)}"
                )

                result.append(activated)

            result = tuple(result)
            print_state(result)

            if (lidx := seen.get(result, None)) is not None:
                if lidx + 1 < iteration:
                    # There was a cycle
                    print(f"!!! Cycle detected starting from iteration {lidx + 1} !!!")

                break

            seen[result] = iteration
            last_result = result

        print("=== RESULT ===")
        print(f"Iterations: {iteration - 1}")
        print(f"Cycle length: {iteration - lidx - 1}")

        # Get last iteration
        final = result

        if direction == "R":
            final = final[::-1]

        print_gene_configuration(final)
        print_state(final)

    eval_state(sys.argv[1])
