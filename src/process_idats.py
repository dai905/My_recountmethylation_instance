#!/usr/bin/env python3

""" process_idats.py

    Authors: Sean Maden, Abhi Nellore
    
    Functions to preprocess idats before being read into minfi.
    Functions:
        * expand_idats: Expand compressed idat files.
"""

import os, sys, re, gzip, shutil
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
import settings; settings.init()

def expand_idats(idatspath, compext = ".*idat.gz$", expext = ".*idat$"):
    """ Detect and expand available idat files.
        
        Arguments:
            * idatspath : Path to instance directory containing downloaded 
                            IDATs (valid file path).
            * compext : Regular expression pattern for extension of compressed
                            IDAT files (string, regex pattern).
            * expext : Regular expression pattern for extension of expanded 
                            IDAT files (string, regex pattern).

        Returns:
            * ridatd dictionary containing expanded IDAT info.

    """
    idats_fnlist = os.listdir(idatspath)
    rexpanded1 = re.compile(expext); rcompressed1 = re.compile(compext)
    idats_fnlist_filt1 = list(filter(rexpanded1.match, idats_fnlist))
    idats_fnlist_filt2 = list(filter(rcompressed1.match, idats_fnlist))
    idats_fnlist_filt = [fn for fn in idats_fnlist_filt2 
                            if not fn in idats_fnlist_filt1]
    ridatd = {} # return dictionary
    print("Expanding "+str(len(idats_fnlist_filt))+"compressed IDATs...")
    for compidat in idats_fnlist_filt:
        idat_fn = os.path.splitext(compidat)[0]
        statuslist = []; ridatd[compidat] = []
        with gzip.open(os.path.join(idatspath, compidat), 'rb') as f_in:
            with open(os.path.join(idatspath, idat_fn), 'wb') as f_out:
                try:
                    shutil.copyfileobj(f_in, f_out); ridatd[compidat].append(1)
                except shutil.Error as se:
                    ridatd[compidat].append(se)
    return ridatd

if __name__ == "__main__":
    """ Process downloaded IDATs

    Process downloaded IDATs for a recountmethylation instance. Prepares IDATs
    for compilation.

    """
    expand_idats(idatspath = settings.idatspath)
