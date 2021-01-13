#!/usr/bin/env python3.8
import sys

sp, input_file, sep = sys.argv

out_content = []
with open(input_file, 'r') as fin:
    for line in fin:
        out_content.append(sep.join([e.strip() for e in line.strip().split(sep)]))

with open(input_file, 'w') as fout:
    fout.write("\n".join(out_content))
