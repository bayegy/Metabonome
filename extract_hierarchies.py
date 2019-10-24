#!/usr/bin/env python3

import sys
import json

input_json = sys.argv[1]
out_file = sys.argv[2]

"""
def children(tree):
    children_list = []

    def get_all_children(tree):
        for child in tree['children']:
            if 'children' not in child:
                children_list.append(child['name'].strip().split(' ')[0])
            else:
                get_all_children(child)
    get_all_children(tree)
    return children_list
"""


def get_all_children(tree):
    for child in tree['children']:
        if 'children' not in child:
            yield child['name'].strip().split(' ')[0]
        else:
            yield from get_all_children(child)


with open(input_json) as infile, open(out_file, 'w') as outfile:
    root = json.load(infile)["children"]
    for tree in root:
        for child in get_all_children(tree):
            outfile.write("{}\t{}\n".format(child, tree['name']))
