#!/usr/bin/env python3
""" Mod notes: 20230424 b
    Authors: Liuhan
    Notes: Create a new .py named rsheet_Win due to hardlink limitations on os

"""

""" Mod notes: 20230424 a
    Authors: Liuhan
    Notes: Make hardlinks in the EPIC_idat data folder, a symbolic link to folder on lab server

    Details:
        1) change to absolute directory
        2) Save *.hlink* data in EPIC_idat data folder which is a link to data folder on lab server
        3) delete redundant scripts and variables
        4) fix .idat files with bad names - meaning without 'GSM******.' at the very beginning. 
    e.g. GSM6723534_203286230174_R06C01_Red.idat
        These files are not associated with hard links yet since they can not pass
    {
     gsm_idats = [i for i in instpath_idatspathlist
        if i.split(".")[0] == gsmid and
        os.path.exists(os.path.join(targetpath, i))]
     }

    Procedures:
        1) First of all, make sure our cwd is always recountmethylation_instance
        2) Specify their directory using absolute dir
        3) Store them in EPIC_datafolder
    Functions:
        1) rmdb_fpaths_EPIC: Make hardlinks to idat files.
        2) new_idat_hlinks_EPIC:  Make hardlinks to idat files

"""

""" rsheet.py

    Author: Sean Maden

    Description:
    Prepare an R sheet for analysis, using RMDB.

    Functions:
    * rmdb_files_for_processing: compile valid file info from RMDB
    * compile_rsheet: compile R sheet from files and RMDB pointers.

"""

import pymongo, sys, os, datetime, inspect, re, json
os.chdir('C:\\Users\\liuhand\\Documents\\GitHub\\My_recountmethylation_instance\\src')
sys.path.insert(0, os.path.join("src"))
import settings; settings.init()
from utilities import gettime_ntp, getlatest_filepath, get_queryfilt_dict

def new_idat_hlinks_EPIC(gsmid, ts, igrn_fn, ired_fn):
    """ new_idat_hlinks_EPIC

        Make new hlink files and return new hlink paths.

        Arguments
        * gsmid: Valid sample GSM ID (str).
        * ts: Timestamp (int or str).
        * igrn_fn: Filename of green intensity data file (IDAT, str).
        * ired_fn: Filename of red intensity data file (IDAT, str).

        Returns
        * rlist (list): List of path to new grn [0] and red [1] hlinked idats
    """
    targetpath = 'Z:\\Liuhan Dai\\DNA Methylation Detection\\20230408 figure 6d\\recountmethylation_EPIC_idat'
    print("generating new hlink filenames...")
    gnewhlinkfn = '.'.join([gsmid, ts, 'hlink',
        '.'.join(igrn_fn.split('.')[2:])
        ]
    )
    rnewhlinkfn = '.'.join([gsmid, ts, 'hlink',
        '.'.join(ired_fn.split('.')[2:])
        ]
    )
    # create new hardlink
    print("make new hlinks with filenames...")
    rlist = []
    # Alternatively we can use win32file module to create hardlink
    grn_hlink = os.link(os.path.join(targetpath,igrn_fn),
            os.path.join(targetpath, gnewhlinkfn)
        )
    red_hlink = os.link(os.path.join(targetpath,ired_fn),
            os.path.join(targetpath, rnewhlinkfn)
        )
    rlist.append(os.path.join(targetpath, gnewhlinkfn))
    rlist.append(os.path.join(targetpath, rnewhlinkfn))
    return rlist

def new_idat_hlinks(gsmid, ts, igrn_fn, ired_fn):
    """ new_idat_hlinks

        Make new hlink files and return new hlink paths.

        Arguments
        * gsmid: Valid sample GSM ID (str).
        * ts: Timestamp (int or str).
        * igrn_fn: Filename of green intensity data file (IDAT, str).
        * ired_fn: Filename of red intensity data file (IDAT, str).

        Returns
        * rlist (list): List of path to new grn [0] and red [1] hlinked idats
    """
    print("generating new hlink filenames...")
    gnewhlinkfn = '.'.join([gsmid, ts, 'hlink',
        '.'.join(igrn_fn.split('.')[2:])
        ]
    )
    rnewhlinkfn = '.'.join([gsmid, ts, 'hlink',
        '.'.join(ired_fn.split('.')[2:])
        ]
    )
    # create new hardlink
    print("make new hlinks with filenames...")
    rlist = []
    grn_hlink = os.link(os.path.join(settings.idatspath,igrn_fn),
            os.path.join(settings.idatspath, gnewhlinkfn)
        )
    red_hlink = os.link(os.path.join(settings.idatspath,ired_fn),
            os.path.join(settings.idatspath, rnewhlinkfn)
        )
    rlist.append(os.path.join(settings.idatspath, gnewhlinkfn))
    rlist.append(os.path.join(settings.idatspath, rnewhlinkfn))
    return rlist

def rmdb_fpaths_EPIC():
    """ rmdb_fpaths_EPIC

        Get filepaths for existant sample idats and msrap outfiles.

        Returns:
        * hlinklist, list of new hlink files created at settings.idatspath.

    """
    # Fix cwd as follows:
    if os.getcwd() != 'C:\\Users\\liuhand\\Documents\\GitHub\\My_recountmethylation_instance\\src':
        os.chdir('C:\\Users\\liuhand\\Documents\\GitHub\\My_recountmethylation_instance\\src')

    timestamp = gettime_ntp()

    targetpath = 'Z:\\Liuhan Dai\\DNA Methylation Detection\\20230408 figure 6d\\recountmethylation_EPIC_idat'


    # list all previously expanded idat files directy from idats dir
    instpath_allidats = os.listdir(os.path.join(targetpath))

    # We discover a bunch of .idat files that are not properly named somehow. Let us fix them.
    bad_allidats = [i for i in instpath_allidats if len(i.split('.')) <= 2] # find those bad ones

    if bad_allidats:
    # loop through all files that have bad names
        for i in bad_allidats:
            # construct new names
            new_name = i.split('_')[0] + '.' + i
            # rename the file
            os.rename(os.path.join(targetpath, i), os.path.join(targetpath, new_name))
            
        print('find .*idat.* files with bad names and fix them......')
        print('Please reload files')
        return

    # expanded idat hlinks
    instpath_hlink = list(filter(re.compile('.*hlink.*').match,
        instpath_allidats))
    inst_hlinkID = [i.split('.')[-2] for i in instpath_hlink]
    # expanded idats that do not contain .hlink and .gz substring
    instpath_expidat = [i for i in instpath_allidats
        if '.hlink' not in i and '.gz' not in i]

    # expanded idats without hlinks where it doesn not have GSMID and Color Code in instpath_expidat
    instpath_nohlink = [i for i in instpath_expidat
        if i.split('.')[-2] not in inst_hlinkID]
    print("Detected " + str(len(instpath_expidat)) +
        " expanded IDATs, and " + str(len(instpath_nohlink)) +
        " expanded IDATs without hlinks " + str(len(instpath_hlink)) +
        " expanded IDATs with hlinks.")
    print("Getting GSM IDs for IDATs without hlinks...")

    gsmlist = list(set([i.split(".")[0] for i in instpath_nohlink]))

    instpath_idatspathlist = [i for i in instpath_nohlink
        if i.split(".")[0] in gsmlist]; hlinklist=[]
    
    n = 0; nn = 0;
    for gsmid in gsmlist:
        print("Processing GSM ID " + gsmid + "...");
        ired_fn = ""; igrn_fn = ""; basename_grn = ""; basename_red = ""
        gsm_idats = [i for i in instpath_idatspathlist
            if i.split(".")[0] == gsmid and
            os.path.exists(os.path.join(targetpath, i))]
        try:
            if len(gsm_idats) >= 2:
                igrn_fn = list(filter(re.compile(".*Grn\.idat$").match, gsm_idats))[0]
                ired_fn = list(filter(re.compile(".*Red\.idat$").match, gsm_idats))[0]
                basename_grn = "_".join(igrn_fn.split(".")[-2].split("_")[0:-1])
                basename_red = "_".join(ired_fn.split(".")[-2].split("_")[0:-1])
                if (basename_grn==basename_red and not basename_grn == ""
                    and not basename_red == ""):
                    print("Making new IDAT hlinks for GSM ID " + gsmid)
                    rlist = new_idat_hlinks_EPIC(gsmid = gsmid, ts = timestamp,
                        igrn_fn = igrn_fn, ired_fn = ired_fn)
                    hlinklist.append(rlist)
                else:
                    nn = nn + 1
            else:                
                n = n + 1
        except:
            print("Couldn't find Red and Grn IDATs for GSM ID " + gsmid)
        print("Finished with GSM ID " + gsmid)
    print("Made " + str(len(hlinklist)) + " new IDAT hlinks. Returning...")
    print(str(n) + ' gsm_idats only have one color channel.')
    print(str(nn) + ' gsm_idats do not have concensus name format.')
    return hlinklist


if __name__ == "__main__":
    """ rsheet.py

        Coordinate MongoDB documents on valid IDATs and SOFT files.

        First generates any missing hlink files.

        Second, passes to compile_rsheet a dictionary of GSM IDs corresponding
        to newly created hlink files for validation.

    """
    try:
        hlinklist = rmdb_fpaths_EPIC()
        #hlinklist = rmdb_fpaths()
        #gsmlist = [os.path.basename(i[0]).split(".")[0] for i in hlinklist]
        #if gsmfpd:
        #    print("Successfully ran rmdb_files_for_processing(). Continuing...")
        #    rsheet = compile_rsheet(gsmdict)
        #    if rsheet:
        #        print("Successfully ran compile_rsheet() on new gsmdd dict.")
        #print("Completed all tasks!")
    except:
        print("Error running rsheet.py! Processes not completed.")
