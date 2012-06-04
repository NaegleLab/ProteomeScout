#!/usr/local/bin/python
from pathname import *
import cgi
import cgitb;cgitb.enable()
import pickle
import math
import sets

print "Content-Type: text/html\n"

form = cgi.FieldStorage()
files = []
clusters = []
for item in form.keys():
    files.append(form.getvalue(item,''))
    clusters.append(item[-1])


allSetsFeatures ={}
allKeys = sets.Set()
i = 0
for file in files:
    print "Cluster %s<BR>"%clusters[i]
    i+= 1
    f = open(entropyPath+file,'rb')
    dict = pickle.load(f)
    f.close()
    
    featureList = []
    for key in dict.keys():
        featureList.append((dict[key],key))
    featureList.sort()
    featureList.reverse()
    for item in featureList:
        print item[1]+' '+str(item[0])+"<BR>"
    print "<BR><BR>"
    ## allSetsFeatures[file]=dict
##     allKeys=allKeys.union(sets.Set(dict.keys()))
    

## entropyList = []
## for item in allKeys:
##     num = 0
##     entropy = 0
##     for dict in allSetsFeatures.values():
##         p = dict.get(item,[0,0])
##         p,num = p[0],max(num,p[1])
##         if p !=0:
##             entropy += -1*p*math.log(p)
##     entropyList.append((entropy,item,num))

## entropyList.sort()
## print "Experimental Entropy Feature<br>"
## print "Feature: entropy, where entropy is a measure of how separating a feature is among clusters. A lower entropy score means the feature is more discriminating across subsets<br>"
## for item in entropyList: print "<br>"+item[1]+' '+str(item[0])+' ' + str(item[2])
    
