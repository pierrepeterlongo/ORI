#!/usr/bin/env python3

import json
import sys



# def reformat_head_name(name):
#     """
#     Convert:
#     head_11368_GCA_000091005_1_fasta
#     into:
#     head_11368_GCA_000091005.1.fasta
#     """
#     if name.endswith("_fasta"):
#         name = name[:-6]  # remove "_fasta"
#     name = name.replace("_1", ".1").replace("_2", ".2")
#     return name + ".fasta"


# def convert_score_to_fraction(score, total):
#     """
#     Convert proportion (float) into integer fraction string.
#     Example:
#     0.627906976744186 -> 54/86
#     """
#     numerator = round(score * total)
#     return f"{numerator}/{total}", numerator


def transform_json(data):
    for dataset_name, queries in data.items():
        for query_id, hits in queries.items():

            # Header line
            print(f"*query{query_id} {len(hits)}")

            # Sort hits by score (descending)
            sorted_hits = sorted(
                hits.items(),
                key=lambda x: x[1],
                reverse=True
            )

            for head_name, score in sorted_hits:
                # fraction_str, numerator = convert_score_to_fraction(score, TOTAL_LENGTH)
                # reformatted_name = reformat_head_name(head_name)

                print(f"{head_name} undex/undex {score:.6f}")


def main():
    if len(sys.argv) != 2:
        print("Usage: python transform.py input.json")
        sys.exit(1)

    input_file = sys.argv[1]

    with open(input_file, "r") as f:
        data = json.load(f)

    transform_json(data)


if __name__ == "__main__":
    main()
