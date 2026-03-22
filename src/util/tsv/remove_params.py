#!/usr/bin/env python3
import sys

def remove_blocks(input_file, output_file):
    with open(input_file) as f:
        lines = f.readlines()

    result = []
    skip = False
    delimiter_count = 0
    for line in lines:
        if "**********" in line and delimiter_count < 2:
            delimiter_count += 1
            skip = not skip
            continue
        if not skip:
            result.append(line)

    with open(output_file, "w") as f:
        f.writelines(remove_empty_lines(result))

def remove_empty_lines(lines):
    return [line for line in lines if line.strip() != ""]

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.tsv output.tsv")
        sys.exit(1)
    remove_blocks(sys.argv[1], sys.argv[2])
