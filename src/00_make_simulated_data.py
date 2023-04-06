####################################################
#                      readme                      #  
####################################################

# make parameter files for simulation
# then run plink

####################################################
#                        env                       #  
####################################################

import system
import glob
import os

####################################################
#                     variables                    #  
####################################################

# outputs
param_path = 'params/'
data_path = 'simulated_data/'

# variables
ncases = ['250', '1000', '10000']
ncontrols = ['1000', '10000']
afs = [[0.01, 0.02], [0.05, 0.06], [0.25, 0.35], [0.4, 0.49]]
odds = {'dom' : [[1, 1/1.1], [1,1/1.25], [1, 1/1.5], [1, 1/2]],
        'add' : [[0.5, 1/1.1], [0.5, 1/1.25], [0.5, 1/1.5], [0.5, 1/2]],
        'rec' : [[1.1, 1.1], [1.25, 1.25], [1.5, 1.5], [2, 2]],
        'mult' : [[1.1, "mult"], [1.25, "mult"], [1.5, "mult"], [2, "mult"]]}
prevalence = [0.01, 0.1, 0.3, 0.5]

nsims = 1
####################################################
#                        main                      #  
####################################################

# make params
for key in odds.keys():
    vals = odds[key]
    for val in vals:
        valx = [round(x, 2) for x in val if x != 'mult']
        for af in afs:
            fw = open(param_path + str(key) + '_maf_' + '-'.join([str(x) for x in af]) + '_or_' + '-'.join(str(x) for x in valx) + '.param', 'w')
            fw.write('9900 null 0.01 0.49 1 1\n')
            fw.write('100 disease ' + ' '.join([str(x) for x in af]) + ' ' + ' '.join(str(x) for x in valx) + '\n')
            fw.close()


param_files = glob.glob(param_path + "*param")

# run plink
for pf in param_files:
    for cacount in ncases:
        for cocount in ncontrols:
            for prev in prevalence:
                for i in range(0,nsims): 
                    out_file = data_path + os.path.basename(pf) + '_' + str(cacount) + '-' + str(cocount) + '_' + str(prev) + '-' + str(i)
                    print(out_file)
                    os.system('./plink --simulate ' + pf +' --simulate-ncases ' + str(cacount) + ' --simulate-ncontrols ' + str(cocount) + ' --simulate-prevalence ' + str(prev) + ' --recodeA --out ' + out_file)


