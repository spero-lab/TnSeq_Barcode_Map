#!/usr/bin/env python 3

# DEFINING FUNCTIONS 
def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter) -33

def hamdist_qs(qs_line: str,ham_dist_threshold: int, qs_threshold: int):
    '''
    Assesses whether or not current qs_line meets qs-threshold depending on hamming_distance (# of pos.s permitted to fall below qs-threshold)
    Input: qs_line (str)[normally from R2 or R3 for demultiplex.], ham_dist threshold, qs_threshold 
    Output: a boolean (T/F) saying whether or not given qs_line passed the quality-score threshold 
    #IMPORTANT: REQUIRES bioinfo.convert_phred ( a fxn that converts a qs-line from ASCII-->actual phred scores)
    '''
    #Converting QS line (ASCII-->Phred (assumes +33 encoding))
    conv_qs=[]
    for letter in qs_line:
        phred_score=convert_phred(letter)
        conv_qs.append(phred_score)
    #Calc. Ham_Dist (a counter for # of bases w/ a qs that fell BELOW threshold)
        #Init. Ham_Dist 
    ham_dist=int(0)
        #cac. ham_dist based on qs threshold 
    for score in conv_qs:
        if score <=qs_threshold:
            ham_dist+=1
        #Filtering based on threshold
    if ham_dist >= ham_dist_threshold:
        return False
    elif ham_dist < ham_dist_threshold:
        return True

def hamdist_barcode(seq_line: str,ham_dist_threshold: int,tranposon_barcode: str):
    '''
    Checks if transposon barcode is in the sequence, with leniency for mismatches
    Input: seq_line (str), ham_dist_threshold (max # of mismatched bases permitted to still be considered a "matched read"), transposon_barcode (LAST ~12 bases of transposon border seq.)
    Output: 
        - if tranposon seq. detected, will return the 0-based POS where the detected tranposon-border seq. ENDS in the seq. 
        - if tranposon seq. NOT detected, will return 'None' 
    '''
    # Looping over substrings the length of the transposon_barcode (sliding a k-mer window across the seq.) --> looking for a match 
        # Calc. k of kmer-window = length of barcode 
    barcode_len=len(tranposon_barcode)
        # Loop through each kmer-window in seq.
    for i in range(len(seq_line)-barcode_len+1): # +1 to make it 1-based
        barcode_window = seq_line[i:i+barcode_len]
        # Calc. Ham_Dist (# of mismatched bases) 
            # Init. Ham_Dist 
        ham_dist=int(0)
            # FOR every aligned pos. of a given barcode_window --> calc. ham_dist (# of mismatched bases)
                # zip(a,b) --> returns a tuple for EVERY aligned position b/w strings a & b
        for a, b in zip(barcode_window, tranposon_barcode):
            #IF the tuple's bases mismatch --> Increment ham_dist +1 for given barcode_window 
            if a != b:
                ham_dist+=1
        # If hamming distance passes threshold (doesn't exceed permitted # of mismatches), return & (2) the position (0-based) where transposon barcode ends 
        if ham_dist <= ham_dist_threshold:
            return i + barcode_len #(0-based cutoff pos.)
    # If none of the kmer windows pass ham_dist_threshold, return None (False) 
    return None 