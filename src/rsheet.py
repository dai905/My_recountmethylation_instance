#!/usr/bin/env python3

""" Mod notes: 20230424
    Authors: Liuhan
    Notes: Make hardlinks in the EPIC_idat data folder, a symbolic link to folder on lab server
    
    Details:
        1) change to absolute directory
        2) Save *.hlink* data in EPIC_idat data folder which is a link to data folder on lab server
        3) delete redundant scripts and variables
    
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
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
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
    grn_hlink = os.link(os.path.join('EPIC_idat',igrn_fn),
            os.path.join('EPIC_idat', gnewhlinkfn)
        )
    red_hlink = os.link(os.path.join('EPIC_idat',ired_fn),
            os.path.join('EPIC_idat', rnewhlinkfn)
        )
    rlist.append(os.path.join('EPIC_idat', gnewhlinkfn))
    rlist.append(os.path.join('EPIC_idat', rnewhlinkfn))
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
    if os.getcwd() != '/home/liuha/recountmethylation_instance':
        os.chdir('/home/liuha/recountmethylation_instance')
        
    timestamp = gettime_ntp()

    # list all previously expanded idat files directy from idats dir
    instpath_allidats = os.listdir(os.path.join('EPIC_idat'))
    
   
    # expanded idats
    instpath_expidat = list(filter(re.compile('.*\.idat$').match, 
        instpath_allidats))
    # idat hlinks
    instpath_hlink = list(filter(re.compile('.*hlink.*').match, 
        instpath_expidat))
    # expanded idats without hlinks
    instpath_nohlink = [i for i in instpath_expidat 
        if not i in instpath_hlink]
    print("Detected " + str(len(instpath_expidat)) + 
        " expanded IDATs, and " + str(len(instpath_nohlink)) + 
        " expanded IDATs without hlinks.")
    print("Getting GSM IDs for IDATs without hlinks...")
    gsmlist = list(set([i.split(".")[0] for i in instpath_nohlink]))
    
    instpath_idatspathlist = [i for i in instpath_nohlink 
        if i.split(".")[0] in gsmlist]; hlinklist=[]
    for gsmid in gsmlist:
        print("Processing GSM ID " + gsmid + "..."); 
        ired_fn = ""; igrn_fn = ""; basename_grn = ""; basename_red = ""
        gsm_idats = [i for i in instpath_idatspathlist
            if i.split(".")[0] == gsmid and 
            os.path.exists(os.path.join('EPIC_idat', i))]
        try:
            igrn_fn = list(filter(re.compile(".*Grn\.idat$").match, gsm_idats))[0]
            ired_fn = list(filter(re.compile(".*Red\.idat$").match, gsm_idats))[0]
            basename_grn = "_".join(igrn_fn.split(".")[2].split("_")[0:-1])
            basename_red = "_".join(ired_fn.split(".")[2].split("_")[0:-1])
            if (basename_grn==basename_red and not basename_grn == "" 
                and not basename_red == ""):
                print("Making new IDAT hlinks for GSM ID " + gsmid)
                rlist = new_idat_hlinks_EPIC(gsmid = gsmid, ts = timestamp, 
                    igrn_fn = igrn_fn, ired_fn = ired_fn)
                hlinklist.append(rlist)
        except:
            print("Couldn't find Red and Grn IDATs for GSM ID " + gsmid)
        print("Finished with GSM ID " + gsmid)
    print("Made " + str(len(hlinklist)) + " new IDAT hlinks. Returning...")
    return hlinklist

def rmdb_fpaths():
    """ rmdb_fpaths

        Get filepaths for existant sample idats and msrap outfiles.

        Returns:
        * hlinklist, list of new hlink files created at settings.idatspath.

    """
    timestamp = gettime_ntp()
    # connect to RMDB mongodb
    #client = pymongo.MongoClient(mdb_host, mdb_port)
    #mdbcon = client.recount_methylation; mdb_idatscon = mdbcon.gsm.idats
    #mdb_idatrecords = list(mdb_idatscon.find())
    # list all previously expanded idat files directy from idats dir
    instpath_allidats = os.listdir(settings.idatspath)
    # compressed idats
    instpath_compidats = list(filter(re.compile('.*\.idat.gz$').match, 
        instpath_allidats))
    # expanded idats
    instpath_expidat = list(filter(re.compile('.*\.idat$').match, 
        instpath_allidats))
    # idat hlinks
    instpath_hlink = list(filter(re.compile('.*hlink.*').match, 
        instpath_expidat))
    # expanded idats without hlinks
    instpath_nohlink = [i for i in instpath_expidat 
        if not i in instpath_hlink]
    print("Detected " +str(len(instpath_compidats))+ 
        " compressed IDATs, " + str(len(instpath_expidat)) + 
        " expanded IDATs, and " + str(len(instpath_nohlink)) + 
        " expanded IDATs without hlinks.")
    print("Getting GSM IDs for IDATs without hlinks...")
    instpath_nohlink_gsm = list(set([i.split(".")[0] for i in instpath_nohlink]))
    print("Getting IDAT filepaths from MongoDB records, "+
        "for GSM IDs lacking hlinks...")
    #mdb_nohlink_gsmlist = list(set([i["gsmid"] for i in mdb_idatrecords 
    #    if i["gsmid"] in instpath_nohlink_gsm]))
    gsmlist = list(set([i.split(".")[0] for i in instpath_nohlink]))
    instpath_idatspathlist = [i for i in instpath_nohlink 
        if i.split(".")[0] in gsmlist]; hlinklist=[]
    for gsmid in gsmlist:
        print("Processing GSM ID " + gsmid + "..."); 
        ired_fn = ""; igrn_fn = ""; basename_grn = ""; basename_red = ""
        gsm_idats = [i for i in instpath_idatspathlist
            if i.split(".")[0] == gsmid and 
            os.path.exists(os.path.join(settings.idatspath, i))]
        try:
            igrn_fn = list(filter(re.compile(".*Grn\.idat$").match, gsm_idats))[0]
            ired_fn = list(filter(re.compile(".*Red\.idat$").match, gsm_idats))[0]
            basename_grn = "_".join(igrn_fn.split(".")[2].split("_")[0:-1])
            basename_red = "_".join(ired_fn.split(".")[2].split("_")[0:-1])
            if (basename_grn==basename_red and not basename_grn == "" 
                and not basename_red == ""):
                print("Making new IDAT hlinks for GSM ID " + gsmid)
                rlist = new_idat_hlinks(gsmid = gsmid, ts = timestamp, 
                    igrn_fn = igrn_fn, ired_fn = ired_fn)
                hlinklist.append(rlist)
        except:
            print("Couldn't find Red and Grn IDATs for GSM ID " + gsmid)
        print("Finished with GSM ID " + gsmid)
    print("Made " + str(len(hlinklist)) + " new IDAT hlinks. Returning...")
    return hlinklist

def compile_rsheet(gsmfpathdict):
    """ compile_rsheet

        Takes dictionary of GSM IDs. Compiles valid GSM IDs, filenames, and 
        path into an rsheet object.
        
        Arguments:
            * gsmfpathdict: gsm paths dict obj output from rmdb_fpaths() 
                (dict).
        
        Returns:
            * lsheet, produces rsheet file as a side effect.

    """
    timestamp = gettime_ntp()
    print("Getting equery filter...")
    eqd = get_queryfilt_dict()
    gsmvalidlist = list(set([gsmid for gselist in list(eqd.values()) 
        for gsmid in gselist
    ]))
    sheetspath = settings.sheetspath; sheetfn_ext = settings.sheetfnstem
    os.makedirs(sheetspath, exist_ok = True)
    sheets_fpath = os.path.join(sheetspath, ".".join([timestamp, sheetfn_ext]))
    # table written as list of row strings
    print("Forming table list for rsheet..."); lsheet = []
    lsheet.append(" ".join(["gsmid",
        "gseid",
        "idats_fn",
        "msrapmd_fn",
        "msrapmd_flatjson",
        "SENTRIX_ID",
        "ARRAY_ID",
        "Basename"]))
    lsheet[0] = lsheet[0]+"\n"; print("Forming filtered GSM dictionary...")
    gsmvalid_fpathlist = {key:value for (key,value) in gsmfpathdict.items() 
        if key in gsmvalidlist}
    if gsmvalid_fpathlist:
        print("Starting iterations on gsm filepaths list of len = "
            +str(len(list(gsmvalid_fpathlist.keys()))))
        for gsmindex, gsmid in enumerate(gsmvalid_fpathlist, 1):
            print("Beginning GSM num "+str(gsmindex)+", id: "+str(gsmid))
            gsmvalid_fp = [fp for fp in gsmvalid_fpathlist[gsmid] 
                if not fp==None
                and not fp==False]
            if gsmvalid_fp:
                print("Getting GSE ID...")
                gseid = ';'.join(list(set([gsek for gsek in list(eqd.keys()) 
                            if gsmid in eqd[gsek]
                            ]
                        )      
                    )
                )
                print("GSE id found: "+str(gseid))
                gsm_fpaths = gsmvalid_fp
                gsmi_redidatpath = [fp for fp in gsm_fpaths if "_Red.idat" in fp]
                gsmi_grnidatpath = [fp for fp in gsm_fpaths if "_Grn.idat" in fp]
                gsmi_msrappath = [fp for fp in gsm_fpaths if "soft.msrapout" in fp]
                if gsmi_redidatpath and gsmi_grnidatpath:
                    print("Found paired channel idats for GSM...")
                    gsmi_redidatpath = gsmi_redidatpath[0]
                    gsmi_grnidatpath = gsmi_grnidatpath[0]
                    # idat filenames
                    grn_idatfn = os.path.basename(gsmi_grnidatpath)
                    red_idatfn = os.path.basename(gsmi_redidatpath)
                    # sample basename (common stem of channel array filenames)
                    print("Getting sample basename..."); gsmimdd = []
                    gsmi_basename = "_".join(red_idatfn.split("_")[0:3]) # basename
                    if gsmi_msrappath:
                        print("Detected metadata file for GSM.")
                        gsmi_msrappath = gsmi_msrappath[0]
                        gsmi_msrappath_var = os.path.basename(gsmi_msrappath)   
                        print("Forming flattened sample metadata...")
                        try:
                            # load msrap mapped terms json file
                            with open(gsmi_msrappath, 'r') as msrapmd:
                                gsmimdd = json.load(msrapmd)
                        except json.decoder.JSONDecodeError:
                            print("Error, cannot load non-json file: "+gsmi_msrappath)
                            gsmi_msrappath_var = "NA"
                    else:
                        gsmi_msrappath_var = "NA"
                    if gsmimdd and not gsmi_msrappath_var == "NA":
                        gsmi_md = gsmimdd[0]; gmd = []
                        # coerce json metadata to flat string
                        for key in list(gsmi_md.keys()):
                            print(str(key))
                            kval = ''
                            if key == 'sample type':
                                print("key = 'sample type'")
                                print("kval = "+str(gsmi_md[key]))
                                gmd.append(str("sampletype="+str(gsmi_md[key])))
                            if key == 'real-value properties':
                                print("key = "+str(key))
                                kval = gsmi_md[key]
                                print("kval : "+str(kval))
                                for index, val in enumerate(kval):
                                    subkval = kval[index]
                                    print("subkval : "+str(subkval))
                                    gmd.append(str(subkval['property_id']))
                            if key == 'mapped ontology terms':
                                print("key = "+str(key))
                                kval = gsmi_md[key]
                                gmd.append(";".join([term for term in kval]))
                            if key == 'sample-type confidence':
                                print("key = "+str(key))
                                gmd.append('sampletypeconf='+str(gsmi_md[key]))
                        gsmi_mdvar = ";".join(gmd) # long metadata string
                    else:
                        gsmi_mdvar = "NA"     
                    # form table row entry for gsmid as long string
                    print("Adding row to table list...")
                    lgsmi = " ".join([gsmid, # gsm id
                        gseid, # gse id
                        ";".join([red_idatfn,grn_idatfn]), # idat filenames
                        gsmi_msrappath_var, # metadata filename
                        gsmi_mdvar, # flattened json file
                        grn_idatfn.split("_")[-2], # sentrix id
                        grn_idatfn.split("_")[-3], # array id
                        gsmi_basename # minfi path Basename, for arrays
                    ])
                    lgsmi = lgsmi+"\n"
                    lsheet.append(lgsmi)
                else:
                    print("Error: GSM is missing one or more valid filepaths. Continuing...")
    else:
        print("No valid GSM IDs detected. Check idats and pipeline folders.")
        return None
    print("Finished processing the GSM files dictionary, writing new rsheet...")
    with open(sheets_fpath,'+w') as fsheet:
        for gsmitem in lsheet:
            fsheet.write(gsmitem)
    return lsheet


def rmdb_fpaths_old(rmhlinks=False):
    """ rmdb_fpaths
        Get filepaths for existant sample idats and msrap outfiles.
        Arguments:
        * rmhlinks : Whether to remove old hardlinks and form new ones, 
                regardless of whether current hlinks exist (boolean).
        Returns:
        * gsmdocdict (dict.) : Dictionary of validated filepaths.
    """
    timestamp = gettime_ntp()
    # connect to RMDB mongodb
    client = pymongo.MongoClient(settings.rmdbhost, settings.rmdbport)
    dbcon = client.recount_methylation; idatscon = dbcon.gsm.idats
    softcon = dbcon.gse.soft; idatslist = list(idatscon.find())
        # grab unique gsm ids
    idatslist = [record for record in idatslist if 'gsmid' in record.keys()]
    gsmindex = list(set([record['gsmid'] for record in idatslist]))
    print("from idats db, found n = "+str(len(gsmindex))+" gsm ids")
        # fname catch patterns for re
    grnidatcatch = settings.grnidat_expcatch
    redidatcatch = settings.redidat_expcatch
    msrapoutcatch = settings.msrapoutfnpattern
        # filter all records for gsm on most recent update datetime
    gsm_fpaths_dd = {}
    # list all previously expanded idat files directy from idats dir
    allidatslist = os.listdir(settings.idatspath)
    allidatslist = list(filter(re.compile('.*\.idat$').match, allidatslist))
    print("found n = "+str((len(allidatslist)))+" expanded idat filenames...")
    # grab and filter idats and msrap outfiles lists
    if rmhlinks:
        print("Beginning sample iterations with hlink removal.")
    else:
        print("Beginning sample iterations without hlink removal.")
    for gi, gsmid in enumerate(gsmindex, 1):
        print("Getting fpaths for gsm: "+str(gsmid)+", num: "+str(gi), end="\r")
        gsm_fpaths_dd[gsmid] = []
        # all idat records for the GSM id
        recordsgsm = [record for record in idatslist if record['gsmid']==gsmid]
        # filter records by channel type,
        # note most records are for compressed files
        idatsrec_gsmgrn = [record for record in recordsgsm 
            if isinstance(record['date'],datetime.datetime)
            and re.search('.*Grn\.idat.*',os.path.basename(record['filepath']))
        ]
        idatsrec_gsmred = [record for record in recordsgsm 
            if isinstance(record['date'],datetime.datetime)
            and re.search('.*Red\.idat.*',os.path.basename(record['filepath']))
        ]
        if idatsrec_gsmgrn and idatsrec_gsmred:
            # get latest records for each channel
            irec_filtgrn = sorted(idatsrec_gsmgrn, key=lambda k: k['date'])[-1]
            irec_filtred = sorted(idatsrec_gsmred, key=lambda k: k['date'])[-1]
            # valid record file basenames
            igrnrec_bn = os.path.basename(irec_filtgrn['filepath'])
            iredrec_bn = os.path.basename(irec_filtred['filepath'])
            # check for expanded versions of compressed files
            igrn_fn = [fn for fn in allidatslist 
                if igrnrec_bn[:-3] in fn
            ]
            ired_fn = [fn for fn in allidatslist 
                if iredrec_bn[:-3] in fn
            ]
            if igrn_fn and ired_fn:
                igrn_fn = igrn_fn[0]
                ired_fn = ired_fn[0]
                hllist = []
                if rmhlinks:
                    # remove old hard links to sample idats
                    grnhl_torm = [fn for fn in allidatslist
                        if "hlink" in fn 
                        and '.'.join(igrn_fn.split('.')[2:]) in fn
                    ]
                    redhl_torm = [fn for fn in allidatslist
                        if "hlink" in fn 
                        and '.'.join(ired_fn.split('.')[2:]) in fn
                    ]
                    if grnhl_torm:
                        for hlfn in grnhl_torm:
                            os.remove(os.path.join(settings.idatspath, 
                                    hlfn)
                                )
                    if redhl_torm:
                        for hlfn in redhl_torm:
                            os.remove(os.path.join(settings.idatspath, 
                                    hlfn)
                                )
                    # new hlinks
                    hllist = new_idat_hlinks(gsmid, ts=timestamp, 
                            igrn_fn=igrn_fn, ired_fn=ired_fn
                        )
                else:
                    # check if hlinks exist, create new ones otherwise
                    grnhllist = [fn for fn in allidatslist
                        if "hlink" in fn 
                        and '.'.join(igrn_fn.split('.')[2:]) in fn
                    ]
                    redhllist = [fn for fn in allidatslist
                        if "hlink" in fn 
                        and '.'.join(ired_fn.split('.')[2:]) in fn
                    ]
                    # get matching grn and red hlink fn's if they exist
                    status_hlink = None
                    grnfnpass = None
                    redfnpass = None
                    if grnhllist and redhllist:
                        grnhllistfilt = list(set(grnhllist))
                        redhllistfilt = []
                        for ghl in grnhllistfilt:
                            for rhl in redhllist:
                                # check that base array ids identical
                                if ghl[:-9]==rhl[:-9]:
                                    redhllistfilt.append(rhl)
                                else:
                                    redhllistfilt.append("")
                        rhlfiltsub = [rhl[:-9] for rhl in redhllistfilt]
                        grnhllistfilt = [ghl for ghl in grnhllistfilt 
                            if ghl[:-9] in rhlfiltsub]
                        redhllistfilt = [rhl for rhl in redhllistfilt
                            if not rhl==""]
                        if grnhllistfilt and redhllistfilt:
                            grnfnpass = grnhllistfilt[0]
                            redfnpass = redhllistfilt[0]
                            # pass hlinks to return dictionary
                            hllist.append(os.path.join(settings.idatspath, grnfnpass))
                            hllist.append(os.path.join(settings.idatspath, redfnpass))
                        else:
                            # make new hlinks
                            hllist = new_idat_hlinks(gsmid, ts=timestamp, 
                                igrn_fn=igrn_fn, ired_fn=ired_fn)
                    else:
                        # make new hlinks
                        hllist = new_idat_hlinks(gsmid, ts=timestamp, 
                            igrn_fn=igrn_fn, ired_fn=ired_fn)
                # finally, pass listed hlinks to return dictionary
                gsm_fpaths_dd[gsmid].append(hllist[0])  
                gsm_fpaths_dd[gsmid].append(hllist[1])    
            else:
                gsm_fpaths_dd[gsmid].append(None)
                gsm_fpaths_dd[gsmid].append(None)
        else:
            gsm_fpaths_dd[gsmid].append(False)
        # check for valid MetaSRA-pipeline filepaths
        try:
            msraplatest = getlatest_filepath(filepath=settings.gsmmsrapoutpath,
                filestr=gsmid, embeddedpattern=True, tslocindex=0, 
                returntype='returnlist'
            )
            if msraplatest and len(msraplatest)==1:
                gsm_fpaths_dd[gsmid].append(msraplatest[0])
        except:
            gsm_fpaths_dd[gsmid].append(False)
        print("Finished with sample num "+str(gi), end="\r")
    print("Finished sample iterations. Returning...")
    # return gsmid dictionary with lists of filtered results or valid fpaths
    return gsm_fpaths_dd

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
