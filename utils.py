from Bio import SeqIO
from collections import defaultdict
import pandas as pd
import re


def count_len(parStr):
    return parStr.count('M') + parStr.count('S')


def get_lower5_str(pres, pree, rnaLen, ctFile):
    ctPairs = [ map(int, l.strip().split()[4:6][::-1]) 
                for l in open(ctFile, 'rt').readlines()[1:1 + rnaLen] ]    
    lstems5 = [ (pos, pair) for pos, pair in ctPairs if pos < pres and pair > rnaLen / 2 ][::-1]
    # lstems5: pairs from the very bottom of pre-miRNA 5' end to pri-miRNA 5'end
    if len(lstems5) < 2:
        return
    
    lstr5 = ''
    for (up5, up3), (down5, down3) in zip(lstems5[:-1], lstems5[1:]): 
        if up5 == down5 + 1: # in stem
            lstr5 = 'M' + lstr5
        elif up3 + 1 == down3: # if 5' bulge
            bulgeSize = up5 - down5 - 1
            if any([ (ctPairs[p - 1][1] != 0) for p in range(down5 + 1, up5)]): 
                # if the bulge folds by itself
                lstr5 = 'X' * bulgeSize + 'M' + lstr5
            else:
                lstr5 = 'A' * bulgeSize + 'M' + lstr5
        else: # if internal loop
            bulge5 = up5 - down5 - 1
            bulge3 = down3 - up3 - 1
            dist = min(bulge5, bulge3)
            if any([ (ctPairs[p - 1][1] != 0) for p in range(down5 + 1, up5)]): 
                # if 5' side folds by itself
                lstr5 = 'X' * bulge5 + 'M' + lstr5
            elif any([ (ctPairs[p - 1][1] != 0) for p in range(up3 + 1, down3) ]): 
                # if 3' side folds by itself
                lstr5 = 'X' * bulge5 + 'M' + lstr5
            else:
                lstr5 = 'S' * dist + 'A' * (bulge5 - dist) + 'M' + lstr5

    lastMS = max(min(lstr5.find('M'), lstr5.find('S')), 0)
    lstr5 = lstr5[:lastMS] + 'M' + lstr5[lastMS:] # because the last pair does not enter the above loop 
    final5 = lstems5[-1][0]
    lstr5 = 'F' * (final5 - 1) + lstr5
    return lstr5


def get_lower3_str(pres, pree, rnaLen, ctFile):
    ctPairs = [ map(int, l.strip().split()[4:6][::-1]) 
                for l in open(ctFile, 'rt').readlines()[1:1 + rnaLen] ]    
    lstems5 = [ (pos, pair) for pos, pair in ctPairs if pos < pres and pair > rnaLen / 2 ][::-1]
    # lstems5: pairs from the very bottom of pre-miRNA 5' end to pri-miRNA 5'end
    if len(lstems5) < 2:
        return

    lstems3 = [ (pos, pair) for pos, pair in ctPairs if pos >= lstems5[0][1] and 0 < pair < rnaLen / 2 ]
    # lstems3: pairs from the very bottom of pre-miRNA 3' end to pri-miRNA 3'end
    if len(lstems3) < 2:
        return
    
    lstr3 = ''
    for (up3, up5), (down3, down5) in zip(lstems3[:-1], lstems3[1:]):        
        if up3 + 1 == down3: # in stem
            lstr3 = lstr3 + 'M'
        elif up5 == down5 + 1: # if 3' bulge
            bulgeSize = down3 - up3 - 1
            if any([ (ctPairs[p - 1][1] != 0) for p in range(up3 + 1, down3) ]): 
                # if the bulge folds by itself
                lstr3 = lstr3 + 'M' + 'X' * bulgeSize
            else:
                lstr3 = lstr3 + 'M' + 'A' * bulgeSize
        else: # if internal loop
            bulge5 = up5 - down5 - 1
            bulge3 = down3 - up3 - 1
            dist = min(bulge5, bulge3)
            if any([ (ctPairs[p - 1][1] != 0) for p in range(up3 + 1, down3) ]): 
                # if 5' side folds by itself
                lstr3 = lstr3 + 'M' + 'X' * bulge3
            elif any([ (ctPairs[p - 1][1] != 0) for p in range(down5 + 1, up5)]): 
                # if 3' side folds by itself
                lstr3 = lstr3 + 'M' + 'X' * bulge3
            else:
                lstr3 = lstr3 + 'M' + 'A' * (bulge3 - dist) + 'S' * dist

    lastMS = max(lstr3.rfind('M'), lstr3.rfind('S'))
    lstr3 = lstr3[:lastMS+1] + 'M' + lstr3[lastMS+1:] # because the last pair does not enter above loop 
    final3 = lstems3[-1][0]
    lstr3 = lstr3 + 'F' * (rnaLen - final3)
    return lstr3


def get_upper_str(pres, pree, rnaLen, ctFile):
    ctPairs = [ map(int, l.strip().split()[4:6][::-1]) 
                for l in open(ctFile, 'rt').readlines()[1:1 + rnaLen] ]    
    lstems5 = [ (pos, pair) for pos, pair in ctPairs if pos < pres and pair > rnaLen / 2 ][::-1]
    # lstems5: pairs from the very bottom of pre-miRNA 5' end to pri-miRNA 5'end
    if len(lstems5) < 2:
        return
    
    ustems5 = [ (pos, pair) for pos, pair in ctPairs 
                if lstems5[0][0] <= pos < pair and rnaLen - 30 < pos + pair < rnaLen + 30 ]
    # ustems5: pairs from the very bottom of pre-miRNA 5' end to the start of loop
    if len(ustems5) < 2:
        return

    ustr5, ustr3 = '', ''
    for (down5, down3), (up5, up3) in zip(ustems5[:-1], ustems5[1:]):
        if up5 == down5 + 1 and up3 + 1 == down3: # in stem
            ustr5 = ustr5 + 'M'
            ustr3 = 'M' + ustr3

        elif up5 == down5 + 1: # if 3' bulge
            bulgeSize = down3 - up3 - 1
            ustr5 = ustr5 + 'M'
            if any([ (ctPairs[p - 1][1] != 0) for p in range(up3 + 1, down3) ]): 
                # if the bulge folds by itself
                ustr3 = 'X' * bulgeSize + 'M' + ustr3
            else:
                ustr3 = 'A' * bulgeSize + 'M' + ustr3

        elif up3 + 1 == down3: # if 5' bulge
            bulgeSize = up5 - down5 - 1
            ustr3 = 'M' + ustr3
            if any([ (ctPairs[p - 1][1] != 0) for p in range(down5 + 1, up5) ]): 
                # if the bulge folds by itself
                ustr5 = ustr5 + 'M' + 'X' * bulgeSize
            else:
                ustr5 = ustr5 + 'M' + 'A' * bulgeSize

        else: # if internal loop
            bulge5 = up5 - down5 - 1
            bulge3 = down3 - up3 - 1
            if bulge5 < 0 or bulge3 < 0: # if 
                up5 = down5
                up3 = down3
                break
                
            dist = min(bulge5, bulge3)
            if any([ (ctPairs[p - 1][1] != 0) for p in range(down5 + 1, up5)]): 
                # if 5' side folds by itself
                ustr5 = ustr5 + 'M' + 'X' * bulge5
                ustr3 = 'X' * bulge3 + 'M' + ustr3
            elif any([ (ctPairs[p - 1][1] != 0) for p in range(up3 + 1, down3) ]): 
                # if 3' side folds by itself
                ustr5 = ustr5 + 'M' + 'X' * bulge5
                ustr3 = 'X' * bulge3 + 'M' + ustr3
            else:
                ustr5 = ustr5 + 'M' + 'A' * (bulge5 - dist) + 'S' * dist    
                ustr3 = 'S' * dist + 'A' * (bulge3 - dist) + 'M' + ustr3
                
        #if (up5, up3) == ustems5[-1]:
    loop = 'L' * (up3 - up5 - 1)
    return ustr5[1:] + 'M' + loop + 'M' + ustr3[:-1]


def parse_ct_to_str(pres, pree, rnaLen, ctFile):
    lower5 = get_lower5_str(pres, pree, rnaLen, ctFile)
    lower3 = get_lower3_str(pres, pree, rnaLen, ctFile)
    upper = get_upper_str(pres, pree, rnaLen, ctFile)
    if not all([lower5, lower3, upper]):
        return 'None'
    whole = lower5 + upper + lower3
    if len(whole) != rnaLen:
        return 'None'
    return whole


def measure_overhang(parStr, start, end):
    if parStr=='None' or re.search('[XF]', parStr[start-2:start+2]+parStr[end-2:end+2]):
        return -99
    lstr5 = parStr[:start - 1]
    lstr3 = parStr[end:]
    len5 = count_len(lstr5)
    len3 = count_len(lstr3)
    overhang = len5 - len3
    return overhang


