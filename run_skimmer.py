import glob
import os

from subprocess import call

path = '/Volumes/KEK_USB/beast/phase1/data/TPC/'
rfiles = []
for f in glob.glob( os.path.join(path, '*.root')):
    ifile = f
    print f
    names=f.split('/')
    infile_name=names[-1].split('.')
    ofile = str(infile_name[0]) + str('_skim') + str('.root')
    os.system('''root -l -q -b 'skimmer.C(FileName="%s", OutputName="%s")' '''
            % (ifile, ofile))
    raw_input('good so far?')
