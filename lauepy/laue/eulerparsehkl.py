#!/usr/bin/env python
# coding: utf-8

import re


def extract_hkl(filename):
    file_name = open(filename).read()

    find_hkls = re.findall(r"\[\s+\d+\]\s+\(.+?\)\s+\((.+?)\)", file_name)
    hkls = []
    for hkl in find_hkls:
        hkl = re.findall(r"\s?(\-?\d+)", hkl)
        hkl_list = []
        for element in hkl:
            element = int(element)
            hkl_list.append(element)
        hkls.append(hkl_list)

    return hkls


if __name__ == "__main__":
    print(extract_hkl("Index.txt"))
