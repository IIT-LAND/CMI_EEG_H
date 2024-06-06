# _01_cmi_rawplots.py
"""
_01_cmi_rawplots.py

Goes through cleaned raw data and makes some plots so we can visualize the raw data

Example usage:

# Run on laptop and generate scripts only
python _01_cmi_rawplots.py laptop script

# Run on laptop and run batches
python _01_cmi_rawplots.py laptop batch

# Run on server and generate scripts only
python _01_cmi_rawplots.py server script

# Run on server and run batches
python _01_cmi_rawplots.py server batch


--- Written by mvlombardo
"""

# import modules ---------------------------------------------------------------
import os
import sys
import fnmatch


# Parse input arguments --------------------------------------------------------
# First argument tells us whether we are running it on the lab's server or our own laptop
if sys.argv[1]=="server":
    run_on_server = True
elif sys.argv[1]=="laptop":
    run_on_server = False

# Second argument tells if you are running to generate scripts only or to run batches as well
if sys.argv[2]=="script":
    script_only = True
elif sys.argv[2]=="batch":
    script_only = False

# specify number of computational threads for use n MATLAB
maxNumCompThreads = 1


# Directories ------------------------------------------------------------------
if run_on_server:
    rootdir = "/media/DATA/RAW/cmihbn"
else:
    rootdir = "/Users/mlombardo/Dropbox/cmi_test"

datadir = "%s/data/raw" % rootdir
codedir = "%s/code" % rootdir
roottmpdir = "%s/tmp" % rootdir
tmpdir = "%s/tmp" % datadir


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Main code
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# generate my list of subjects from subjects already present in the raw data directory
sublist = os.listdir(datadir)
# sublist.remove('.DS_Store')
sublist.remove('tmp')

# loop over subjects
for sub in sublist:

    # subject-specific paths in raw and preproc directories
    subpath = "%s/%s" % (datadir,sub)

    # generate a list of files inside subpath
    flist = os.listdir(subpath)

    # I want my flist to be just the tasks of interest. Other files like *.sh
    # files might exist in subpath, so here I need to find them and remove them
    # from flist before moving on. The remaining flist list is a list of strings
    # for each task directory
    pattern2find = "*.sh"
    matches = fnmatch.filter(flist, pattern2find)
    for match in matches:
        flist.remove(match)

    # loop over tasks in flist
    for task in flist:

        # start my bash script
        script = []
        script.append("#!/bin/bash")
        script.append("")
        script.append("# %s Plotting" % sub)
        script.append("")
        script.append("# %s" % task)

        # file to preprocess
        datafile = "%s/%s/%s_%s.set" % (subpath,task,sub,task)

        # MATLAB stuff to make the plots
        outfile =  "%s/%s/%s_%s_raw_tsplot.jpg" % (subpath,task,sub,task)
        if os.path.isfile(outfile):
            print("%s exists already. Moving on..." % outfile)
        else:
            matlab_command2run = "cd('%s'); maxNumCompThreads(%d); tic; try; cmi_tsplot('%s'); catch ME; disp(ME.message); toc; exit; end; toc; exit;" % (codedir,maxNumCompThreads,datafile)
            script.append('matlab -nodesktop -r "%s"' % matlab_command2run)

        outfile =  "%s/%s/%s_%s_raw_carpetplot.jpg" % (subpath,task,sub,task)
        if os.path.isfile(outfile):
            print("%s exists already. Moving on..." % outfile)
        else:
            matlab_command2run = "cd('%s'); maxNumCompThreads(%d); tic; try; cmi_carpetplot('%s'); catch ME; disp(ME.message); toc; exit; end; toc; exit;" % (codedir,maxNumCompThreads,datafile)
            script.append('matlab -nodesktop -r "%s"' % matlab_command2run)

        outfile =  "%s/%s/%s_%s_raw_tsnrcarpetplot.jpg" % (subpath,task,sub,task)
        if os.path.isfile(outfile):
            print("%s exists already. Moving on..." % outfile)
        else:
            matlab_command2run = "cd('%s'); maxNumCompThreads(%d); tic; try; cmi_tsnrplot('%s'); catch ME; disp(ME.message); toc; exit; end; toc; exit;" % (codedir,maxNumCompThreads,datafile)
            script.append('matlab -nodesktop -r "%s"' % matlab_command2run)

        # Writing bash script to file ----------------------------------------------
        print("Saving rawplots bash script for %s %s ..." %(sub,task))
        fname2save = "%s/%s/_01_rawplots_%s_%s.sh" % (subpath,task,sub,task)
        fname = open(fname2save,"w")
        fname.write("\n".join(script)+"\n")
        fname.close()

        os.system("mkdir -p %s/_01_rawplots_batches" % codedir)
        fname2save = "%s/_01_rawplots_batches/_01_rawplots_%s_%s.sh" % (codedir,sub,task)
        fname = open(fname2save,"w")
        fname.write("\n".join(script)+"\n")
        fname.close()

        # if you're not just generating scripts and want to run the bash scripts
        # then do this below...
        if not script_only:

            if run_on_server:

                print("Inserting batch into the queue ...")

                # options for sbatch
                outlog = "%s/%s/%s_%s_rawplots_out.log" % (subpath,task,sub,task)
                errlog = "%s/%s/%s_%s_rawplots_error.log" % (subpath,task,sub,task)
                sbatch_opts = "-J %s --tasks-per-node=1 --cpus-per-task=1 --time=02:00:00 --no-requeue --mem=8g --output=%s --error=%s" % (sub,outlog,errlog)

                # run batch
                os.system("sbatch %s %s" % (sbatch_opts,fname2save))

            else:

                # run bash script
                os.system("bash %s" % fname2save)
