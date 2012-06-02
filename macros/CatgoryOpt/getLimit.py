#! /usr/bin/python

import sys,os


def getJobs(dir):
    masses=[ ]
    getCom='ls -l '+dir+' | grep ".log"'
    (fin,fout)=os.popen4(getCom)
    for line in fout.readlines():
        masses.append(float(line.rstrip().rsplit('_')[1].split('.log')[0]))
    return masses

def getResult(dir,mass):
    if abs(mass) >= 1.:
        return -1.
    massName=str(mass).split('.')[-1]
    if (len(massName) < 3):
        massName=massName+'0'
    if (len(massName) < 3):
        massName=massName+'0'
    if (mass<0.) :
        massName='-.'+massName
    else:
        massName='.'+massName
    fileName='job_'+massName+'.log'
    getCom='grep RES '+dir+'/'+fileName
    (fin,fout)=os.popen4(getCom)
    res=fout.readlines()
    if (len(res) == 0):
        return -1.
    theRes=res[-1].rstrip().split()[8]
    if(theRes.find('gxx') > 0):
        return '0.'
    else:
        return theRes

dir=sys.argv[1]
funname=sys.argv[2]

dMasses=getJobs(dir)
dMasses.sort()

print "void "+funname+"(TGraph*& obsG) {"

print "const int nMassD = "+str(len(dMasses))+";"
print "Double_t obsMass[nMassD];"
print "Double_t obsMean[nMassD];"

counter=0
for idx,mass in enumerate(dMasses):
    res=getResult(dir,mass)
    if (str(res) == 'file'):
        continue
    if (res < 0.) :
        continue
    print "obsMass["+str(counter)+"] = "+str(mass)+";"
    print "obsMean["+str(counter)+"] = "+str(res)+";"
    counter=counter+1


print "realN="+str(counter)+";"
print "obsG = new TGraph(realN,obsMass,obsMean);"
print "return; }"
