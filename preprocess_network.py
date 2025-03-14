import re
import itertools


# Use this to index genes
def EXPR_REPLACER(n: int):
    return f"(state >> {n} & 1)"


# Use this to generate functions
def LINE_REPLACER(case, expr: str):
    return f"case {case}:\nreturn {expr};"


lines = open("Network10.txt").readlines()[1:]

# print(lines)

genes = [re.match(r"^\w+(?=,)", line).group(0) for line in lines]

gene_table = dict(zip(genes, itertools.count()))

for case, line in enumerate(lines):
    definition = re.search(r"(?<=,).+", line).group(0)

    # Clean up work for better reading
    definition = re.sub(r"\s+", " ", definition)

    # Use number defs
    definition = re.sub(
        r"\w+", lambda x: EXPR_REPLACER(gene_table[x.group(0)]), definition
    )

    print(LINE_REPLACER(case, definition))
