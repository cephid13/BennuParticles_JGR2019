#!/usr/bin/env python

import subprocess
import string

exeName = '/Applications/MATLAB_R2018a.app/bin/matlab -nodesktop -nodisplay -nojvm -nosplash' 

line1 = 'addpath(\'~/Documents/MATLAB/mice/src/mice/\')'
line2 = 'addpath(\'~/Documents/MATLAB/mice/lib/\' )'
# vmag = 10;
# posCase = 4;
# srpLatD = 0;
# srpLongD = 15;
line7 = 'test_particles_thermal_UA2(vmag,posCase);'

vlist = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]

posCaseUse = 4;

runNum = 41;

for ii in range(-1,10):

    if ii == -1:
      inputName = 'run00th.m'
      outputName = 'out0035th' + str(runNum) + '.txt'
    else: 
      inputName = 'run' + str(ii) + 'th.m'
      outputName = 'out' + str(ii) + '35th' + str(runNum) + '.txt'
    
    with open(inputName, 'w') as f_in:
        f_in.write(line1 + '\n')
        f_in.write(line2 + '\n')
        f_in.write('vmag = %d;\n' % vlist[ii+1])
        f_in.write('posCase = %d;\n' % posCaseUse )
        f_in.write(line7 + '\n')

    subprocess.Popen(exeName + ' < ' + inputName + '>& ' + outputName + ' &',shell=True,stdin=None, stdout=None, stderr=None, close_fds=True)
