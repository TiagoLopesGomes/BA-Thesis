#! /usr/bin/python3


import sys
import re
import time
import numpy as np
import pandas as pd
import getopt


def runningmean(list,
                n):  # function implementing running mean with increasing and decresing slicing window of length n at the beginning and the end of a list
    if n % 2 == 0:
        raise Exception("n musi byc nieparzyste")
    R = np.empty((len(list)))
    sum = 0
    for i in range(int(n / 2) + 1):
        sum += list[i]
    for i in range(len(list) - int(n / 2) - 1):
        R[i] = sum / min(n, i + int(n / 2) + 1)
        sum += list[i + int(n / 2) + 1]
        if i - int(n / 2) >= 0:
            sum -= list[i - int(n / 2)]
    for i in range(len(list) - int(n / 2) - 1, len(list)):
        R[i] = sum / min(n, len(list) - i + int(n / 2))
        sum -= list[i - int(n / 2)]

    return R


def slicer(M):
    # slice the sequence of proteins into idp-like and nonidp fragments
    # M is a matrix from iupred; M=[aa,iu,anch] where aa is aminoacid, iu - idp score; anch - achore score for predicting binding site of idp
    idp = ""
    non = ""
    currentidp = ""
    currentnon = ""

    toidp = False
    tonon = False

    runmean = runningmean(M[1], 15)

    for i in range(len(runmean)):
        if runmean[i] > 0.5:
            currentidp += M[0][i]
            if tonon and len(currentnon) > 10:
                non += currentnon + "\n"
                currentnon = ""
            toidp = True
            tonon = False
        elif runmean[i] < 0.5:
            currentnon += M[0][i]
            if toidp and len(currentidp) > 10:
                idp += currentidp + "\n"
                currentidp = ""
            tonon = True
            toidp = False
    if len(currentnon) > 10:
        non += currentnon + "\n"
    if len(currentidp) > 10:
        idp += currentidp + "\n"
    return [idp, non]


def divider(data, out):
    # this divide sequences to idp and nonidp given data from iupred
    like = open("IDP_" + out + ".fasta", "a+")
    nonlike = open("nonIDP_" + out + ".fasta", "a+")

    for protId, M in data.items():
        if protId == "#":
            continue
        IDP, nonIDP = slicer(M)
        like.write(protId + "\n")
        nonlike.write(protId + "\n")
        like.write(IDP + "\n")
        nonlike.write(nonIDP + "\n")

    like.close()
    nonlike.close()


def readData(file):
    # s = "[\n#]{18}#.*\n#.*\n#.*\n#.*\n|#.*\n#.*\n#.*\n|[\n#]{18}"
    s = "###\n"
    f = open(file).read()
    results = re.split(s, f)
    R = {}
    cnt = 0
    del results[-1]
    for rec in results:
        cnt += 1
        # teraz zajmuje sie jednym białkiem
        rec = rec.split("\n")
        id = rec[0]
        del rec[0]
        del rec[-1]
        aa = np.empty([len(rec)], dtype=str)
        iu = np.empty([len(rec)], dtype=float)
        anch = np.empty([len(rec)], dtype=float)
        for i in range(len(rec)):
            aaline = rec[i].split("\t")
            aa[i] = aaline[1]
            iu[i] = aaline[2]
            anch[i] = aaline[3]
        M = [aa, iu, anch]
        R[id] = M
    return R # returns the dictionary where key is the name of a protein and value a list of data from iupred


def main(argv):

    try:
        opts, args = getopt.getopt(argv[1:], "hi:o:")
    except getopt.GetoptError:
        print('ERROR! Pass correct arguemnts\n', argv[0], '-i <inputfile> -o <outplutfile>')
        sys.exit()

    outputfile = "seqs"
    data = None

    for opt, arg in opts:
        if opt == '-h':
            print("Script takes the data from iupred2a and divide it to two files. In one there are fragments of proteins which were classified as IDP-like,"
                  "and in the other there are fragments classified as nonIDP. One protein sequence can be splited many time, in particular into very small fragments."
                  "So the script takes advantage of running mean to make the predictive function smooth.")
            print(argv[0], '-i <inputfile> -o <outplutfile>')
        if opt == '-i':
            inputfile = arg
            data = readData(inputfile)

        if opt == '-o':
            outputfile = arg

    divider(data, outputfile)


if __name__ == "__main__":
    start = time.time()
    main(sys.argv)
    end = time.time()
    print(end - start)
