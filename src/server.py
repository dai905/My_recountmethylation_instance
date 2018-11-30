#!/usr/bin/env python3

""" server.py
Recount Methylation Server
Authors: Sean Maden, Abhinav Nellore
Website: 

DESCRIPTION
    Server script manages an instance of recount-methylation database.

OVERVIEW
    A recount-methylation instance consists of files (namely edirect query 
    results, experiment metadata in soft format, and methylation array intensity 
    data or 'idat' files) obtained from edirect queries and ftp-called downloads 
    from the Gene Expression Omnibus (GEO). RMDB is a recount-methylation Mongo 
    database that aggregates file metadata as documents, including experiement 
    (GSE) and sample (GSM) ids, ftp addresses and file paths to downloaded 
    files, and a datetime-formatted date corresponding to last file update. 
    Files are versioned using NTP timestamps in filenames. 

    For best results, we recommend users attempt an initial setup of their 
    recount-methylation instance using default generated directory trees and 
    filenames, and do not directly change locations or names of files initially 
    downloaded.

SERVER PROCESSES
    The server.py script manages process queues, error handling, and 
    coordination of recount-methylation. It currently uses Celery distributed 
    task queue to queue jobs synchronously. Jobs are brokered using RabbitMQ, 
    and queue details are backed up locally in a SQLite db.

DEPENDENCIES AND SETUP
    1. Recount-methylation primarily uses Python 3 for download handling and 
        file management. R is used for SOFT-to-JSON conversion, and for 
        preprocessing arrays. MetaSRA-pipeline is used for standardizing 
        experiment metadata, and we recommend installing a fork of the original
        repo (link to GitHub repo). 

    2. To run recount-methylation server.py, ensure the following resources are 
        installed and running:
        * Celery (http://www.celeryproject.org/) 
            -- you may run this in a new terminal or shell window to monitor 
                jobs in realtime (in mac terminal use:
                'celery worker -A gse_celerytask -l INFO' 
                from dir: './recount-methylation-server/src/')
        * RabbitMQ (https://www.rabbitmq.com/)
        * MongoDB (https://www.mongodb.com/)

"""

import subprocess
import glob
import sys
import os
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
import edirect_query
from edirect_query import gsm_query, gse_query, querydict, gsequery_filter


def getlatest_filepath(filepath,filestr):
    """ Get path the latest version of a file, based on timestamp
        Arguments
            * filepath (str) : path to dir to search
            * filestr (str) : pattern of file to search
        Returns
            * latest_file_path (str) : path to latest version of file, OR
            * 0 : search turned up no files at location
    """
    filelist = glob.glob('.'.join([os.path.join(filepath, filestr), '*']))
    if filelist:
        if len(filelist) > 1:
            # sort on timestamp
            filelist.sort(key=lambda x: int(x.split('.')[1]))
            latest_file_path = filelist[-1]
        else:
            latest_file_path = filelist[0]
        return latest_file_path
    else:
        return 0       

def firsttime_run(filedir = 'recount-methylation-files'):
    """ On first setup, run new equeries and query filter
        Arguments
            * filedir (str): dir name for db files 
        Returns
            * gseidlist (list): list of valid GSE IDs
    """
    equery_dest = os.path.join(filedir,'equery')
    temp = os.path.join(filedir,'temp')
    gse_query(dest = equery_dest, temp = temp)
    gsm_query(dest = equery_dest, temp = temp)
    gseqfile = getlatest_filepath(equery_dest,'gse_edirectquery') 
    gsmqfile = getlatest_filepath(equery_dest,'gsm_edirectquery')
    gsequery_filter(gsequery = gseqfile,gsmquery = gsmqfile,
        target = equery_dest,splitdelim = ' '
        )
    gsefiltpath = getlatest_filepath(equery_dest,'gsequery_filt')
    if gsefiltpath and not gsefiltpath == 0:
        gsefiltd = querydict(query = gsefiltpath,splitdelim=' ')
        gseidlist = list(gsefilt.keys())
        print("GSE id list of len "+str(len(gseidlist))+" found. Returning...")
        return gseidlist
    else:
        print("Error retrieving gse query filtered file. Returning...")
        return 0

def scheduled_run(eqfilt_path=False, filedir = 'recount-methylation-files'):
    """ Tasks performed on regular schedule, after first setup
        Arguments
            * filedir : file to search for filtered equery files
            * eqfilt_path : path to GSE equery filtered file
        Returns
            * 0 (int) : If error encountered, or
            * gse_list (list) : list of valid GSE IDs
    """
    eqpath = os.path.join(filedir,'equery')
    gsefilt_latest = getlatest_filepath(eqpath,'gsequery_filt')
    if not gsefilt_latest or gsefilt_latest == 0:
        print("No gse query filt file found, checking for gse and gsm "
            +"queries...")
        # run fresh gse/gsm queries if not found
        gsequery_latest = getlatest_filepath(eqpath,'gse_edirectquery')
        if not gsequery_latest or gsequery_latest == 0:
            gse_query(dest = eqpath,
                temp=os.path.join(filedir,'temp')
                )
        gsmquery_latest = getlatest_filepath(eqpath,'gsm_edirectquery')
        if not gsmquery_latest or gsmquery_latest == 0:
            gsm_query(dest = eqpath,
                temp = os.path.join(filedir,'temp')
                )
        # run new gse filt
        print("Running filter on gse query...")
        gsequery_filter(gsequery=gsequery_latest,
                gsmquery = gsmquery_latest,
                target = eqpath, splitdelim = ' '
                )
        gsefilt_latest = getlatest_filepath(eqpath,'gsequery_filt')
    else:
        print("Gse query filt file found. Continuing...")
    if gsefilt_latest:
        print("Getting dictionary from gse filt file...")
        gsefiltd = querydict(query = gsefilt_latest, splitdelim = ' ')
        gseidlist = list(gsefiltd.keys())
        print("GSE id list of len "+str(len(gseidlist)) + " found. Returning..")
        return gseidlist
    else: 
        print("Error, attempt to produce gseidlist failed. Check for a valid "
            +"gse query filt file and gse/gsm query files.")
        return 0

def get_queryfilt_dict(filesdir='recount-methylation-files',eqtarget='equery'):
    """
    """
    eqpath = os.path.join('recount-methylation-files','equery')
    gsefilt_latest = getlatest_filepath(eqpath,'gsequery_filt')
    if gsefilt_latest and not gsefilt_latest == 0:
        gsefiltd = querydict(query = gsefilt_latest, splitdelim = ' ')
        return gsefiltd
    else:
        print("Error, no gse filtered file found at location.")

def on_restart():
    """ Handling server.py restart, inc. interruptions when processing queue
        Arguments
        Returns
    """

def main(files_dir='recount-methylation-files'):
    """ Script to run on call from cl.
        Arguments
        Returns
    """

if __name__ == "__main__":
    from gse_celerytask import gse_task
    gselist = [] # queue input, gse-based
    qstatlist = [] # job status object, also stored at sqlite db
    files_dir = 'recount-methylation-files'
    if os.path.exists(files_dir):
        print(files_dir+" found. Running scheduled_run...")
        gselist = scheduled_run()
    else:    
        print(files_dir+" not found. Creating dir and running firsttime_run...")
        os.makedirs(files_dir)
        gselist = firsttime_run()
    if gselist:
        gqd = get_queryfilt_dict()
        for gse in gselist:
            qstatlist.append(gse_task(gse_id = gse,
                gsefiltdict = gqd)
            )

    #queuerun = run_gsequeue(gse_list = gselist)
        
        
