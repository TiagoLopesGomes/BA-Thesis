#! /usr/bin/python3

import sys, getopt
import itertools
import numpy as np
import operator as op
import collections
import math
from functools import reduce


def ncr(n, r): # counts Newton
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer / denom

def permutation_groups(n): # it returns all possible permutation groups for given length
    x=itertools.combinations_with_replacement('ACDEFGHIJKLMNPQRSTWY', r=n)
    S = {}
    for i in x:
        S[''.join(i)] = 0
    return S

def readData(file): # reading file with n-meres and their counts
    f = open(file).readlines()
    S={}
    for i in f:
        i=i.split("\t")
        S[i[0]] = int(i[1])
    return S

def countGroupLen(data, n): # it counts total occurances of all seqs in permutation group
    G = permutation_groups(n)
    C = collections.Counter()
    for i in G:
        C[i]=0
    for k,v in data.items():
        k = ''.join(sorted(k))
        C[k] += v
    return C


def calculateM(s): # for given string calculates all possible permutations (m)
    a = list(s)
    n = len(a)
    cnts = [1 for i in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if a[i] == a[j]:
                cnts[i] += 1
    m=1
    j=0
    while n > 1:
        m *= ncr(n, cnts[j])
        n-=cnts[j]
        j+=1
    return m

def calculateZScore(Nobs, n, p): # calculate one z-score
    return round( (Nobs - n * p) * 1.0 /math.sqrt(n * p * (1 - p)), 2)

def ZScores(data, n): # calculate z-scores for all given sequences
    # P = permutation_groups(n)
    C = countGroupLen(data, n)
    Z = {}
    cnt = 0
    for perm, Nobs in data.items():
        cnt += 1
        if cnt % 100000 == 0: print(round(cnt/len(data.items())*100, 1), '%')
        m=calculateM(perm)
        p = 1/m
        group_len = C[''.join(sorted(perm))]
        if n == 0 or p == 0 or p == 1:
            continue
        z = calculateZScore(Nobs, group_len, p)
        Z[perm] = z
    return Z




def main(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "hi:n:o:")
    except getopt.GetoptError:
        print('ERROR! Pass correct arguments\n', argv[0],'-i <inputfile> -n <n_meres> [-o <outplutfile>]')
        sys.exit()
    if '-h' in argv:
        print('''This program is for anylysis of short n-meres-peptides in proteins.\n
                    It returns z-scores of n-meres-peptides based on the size of its permutation group\n
                    Input file should contain nmeres with its counts in dataset.\n
                    Output file is n-mer with its z-score''')
        sys.exit()

    if '-n' not in argv or '-i' not in argv:
        print('ERROR! Pass correct arguments\n', argv[0], '-i <inputfile> -n <n_meres> [-o <outplutfile>]')
        sys.exit()

    outputfile = argv[0] + "_ZScores"

    for opt, arg in opts:
        if opt == '-i':
            inputfile = arg
        if opt == '-n':
            n = int(arg)
        if opt == '-o':
            outplutfile = arg


    data = readData(inputfile)
    Z = ZScores(data, n)
    o = open(outplutfile, 'w')
    for i in sorted(Z.items(), key=lambda kv: kv[1], reverse= True):
        o.write(i[0])
        o.write('\t')
        o.write(str(i[1]))
        o.write('\n')
    # G = countGroupLen(data, 5)
    # for k,v in G.items():
    #     print(k, '\t', v)


if __name__ == '__main__':
    main(sys.argv)