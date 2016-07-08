import glob
import os

### This was an attempt to recreate the NTP/TPC subdirectory structure, but it
### wasn't worth the hassle.
#ifpath = '/ghi/fs01/belle2/bdata/group/detector/BEAST/data/NTP/TPC'
#ofpath = '/home/belle/mhedges/tpc_skimmer/data/'
#for subdirs in os.walk(ifpath):
#    print subdirs[0]
#    raw_input('looking into sub_subdir structure')
#    if subdirs[0].split('/')[-1] == 'TPC':
#        for subdir in subdirs[1]:
#            print 'subdir is ', subdir
#            newdir = str(ofpath) + str(subdir)
#            #os.system('mkdir %s' % newdir)
#            #os.system('ls -lhrt')
#            print newdir
#            raw_input('good so far?')
#    print raw_input('Going on to next subdir, should be a date')

ifpath = '/ghi/fs01/belle2/bdata/group/detector/BEAST/data/NTP/TPC'
ofpath = '/home/belle/mhedges/tpc_skimmer/data/'

good_run = 1
for subdir, dirs, files in os.walk(ifpath):
    for f in files:
        test = subdir.split('/')
        for i in test:
            if i  == 'badtime' or i  == 'old' or i == 'ENV' or i == 'tmp':
                good_run = 0
        if good_run == 0:
            continue
        ifile = os.path.join(subdir, f)
        print ifile
        names=f.split('/')
        infile_name=names[-1].split('.')
        ofile = str('data/') + str(infile_name[0]) + str('_skim') + str('.root')
        if os.path.isfile(ofile):
            continue
        os.system('''root -l -q -b 'skimmer.C(FileName="%s", OutputName="%s")' '''
                % (ifile, ofile))
        #raw_input('good so far?')
