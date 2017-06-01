"""
Run this test once you have installed Dark Sage to ensure everything is working correctly.
This will fetch pre-made Dark Sage output, run the code on your machine, then compare the outputs.
This script should run straight out of the box with "python test.py"
"""

import os
import urllib
import filecmp
import subprocess

dir = 'test/'
if not os.path.exists(dir):
    os.makedirs(dir)

if not os.path.isfile(dir+'model_to_test_against_z2.239_0'):
    zip = urllib.urlretrieve('https://github.com/arhstevens/DarkSageTest/archive/master.zip')
    subprocess.call(['unzip', '-j', zip[0], '-d', dir])
    subprocess.call(['rm', zip[0]])

if os.path.isfile(dir+'model_z2.239_0'):
    subprocess.call(['rm', dir+'model_z2.239_0'])

subprocess.call(['./darksage', dir+'test.par'])

if filecmp.cmp('test/model_to_test_against_z2.239_0', 'test/model_z2.239_0'):
    print 'Success! Dark Sage output matches what is expected!'
else:
    print 'Uh oh! The Dark Sage output did not match what was expected!'
    print 'This can happen if test/test.py or any of the Dark Sage codebase was modified from the main repository.'
    print 'If you recently updated your local repository for Dark Sage, try deleting the `test/\' directory and running this again.'