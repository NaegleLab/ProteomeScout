#!/usr/local/bin/python
from pathname import *
import os,sys
import cgi
import cgitb
if displayPythonErrors:
    cgitb.enable()
import cStringIO

os.environ['MPLCONFIGDIR']='%s'%mpl
import matplotlib
import pylab
import math
import random
from matplotlib.font_manager import *
from matplotlib.text import *
## get form data
form = cgi.FieldStorage()
chart = form.getvalue("chart","line")
group = form.getvalue("group","run")
errors = form.getvalue("errors",'true')
transform = None

dataDictFile = form.getvalue("dict","")
import pickle
try:
    dataDict = pickle.load(open(dataDictFile,'r'))
except IOError: dataDict = {}




# calculate the number of graphs to draw
numgraphs = 1 # get the number of graphs to show, only one with all data


dataDict = {'average': {'xlabels': ['0', '1', '2', '4', '8', '16', '32'], 'yvector': [3, 9, 4, 2, 0, 0, 0], 'xvector': [0, 1, 2, 4, 8, 16, 32], 'yerrors': ['0.080499884', '0.055152418', '0.003439604', '0', '0.035917862', '0.000630166', '0.021569119'], 'title': 'average'}}

group = 'run' # don't allow for group by site if there is only 1 site
title=''



def normToMax(datavector):
    maximum = max(datavector)
    for i in range(len(datavector)):
        if datavector[i]==None: datavector[i]=0 #fix this
    return [item*1.0/maximum for item in datavector]
def transformStdDev(stddev,datavector,transform = None):
    if not transform: return stddev
    if transform == "normToMax":
        for i in range(len(stddev)):
            if stddev[i]==None: stddev[i]=0
        return [item*1.0/max(datavector) for item in stddev]
    else: return stddev
def transformData(datavector,transform = None):
    if not transform:
        return datavector
    if transform == "normToMax":
        return normToMax(datavector)
    else:
        return datavector

## how many rows, columns of graphs
numCols = 1
numRows = math.ceil(numgraphs*1.0/numCols)
##
#if len(msids)>1:
 #   width = 1.0/len(msids)

#else:
    #width = 1.0

def getWidth(data,numids):
    if data.get('xvector') != None:
        xrange = max(data['xvector'])
    else:
        xrange = numids
    return 0.4*numids*1./xrange
width = 0

# make dictionary of colors for msids
msColors = {}
runColors = {}


def allZeros(list):
    allz=True
    for item in list:
        if item !=0:
            allz=False
    return allz
top=0.4
if group != 'site':
    pylab.figure(figsize=(9,9))
    top = 0.9
else:
    pylab.figure(figsize=(10,numgraphs))
    if numgraphs>6:
        top=0.99
plotCount = 1


    
numCols = 1
numRows = math.ceil(numgraphs*1.0/numCols)
            
# if group is run or all
maxSoFar = 0

dictionary = {}
for run in dataDict.keys():
         drewPlot = False
         pylab.subplot(1,1,1)
         dictionary = dataDict[run]
         yerrs=dictionary['yerrors']
         yerrs=transformStdDev(yerrs,dictionary['yvector'],transform)
         yerr=[]
         for item in yerrs:
             item = float(item)
             if item ==None: ## need to deal with missing data point
                 yerr.append(0)
             else:
                 yerr.append(item)
         leny=len(yerr)
         yvector = transformData(dictionary.get('yvector',[]),transform)
         xvector = dictionary['xvector']
         #pylab.plot([1,2,3],[1,2,3])
         if xvector !=None:
             pylab.errorbar(xvector,yvector,yerr=yerr)
         else:
             x = range(len(yvector))
             if chart == 'line' and not allZeros(yvector):
                 pylab.errorbar(x,yvector,yerr=yerr)
             else:
                 if not allZeros(yvector):
                     pylab.bar([(item-(counter*width)*0.3) for item in x],yvector,width,yerr=yerr)
         pylab.title(run)
         ymax = max(yvector)
         if dictionary.get('xlabels',None)!=None:
             xlabels = dictionary['xlabels']
             pylab.xticks(range(len(xlabels)))
             if chart=='bar':
                 pylab.axis([-0.5,max(dictionary['xvector'])+width+1*width,0,ymax])
             else:
                 pylab.axis([-0.5,max(dictionary['xvector']),0,ymax])
                
         else:
             xlabels = tuple(dictionary.get('xlabels',[]))
             pylab.xticks(range(len(xlabels)))
             text='\n'
             for i in range(len(xlabels)):
                 text = text + str(i)+': '+xlabels[i]
                 if i != len(xlabels)-1: text += '\n'
             pylab.text(len(xlabels)-1,ymax,text,bbox=dict(facecolor='pink', alpha=0.7),zorder=4)
             #pylab.xticks(pylab.arange(len(yvector)),tuple(dictionary['xlabels']),fontsize='xx-small')
             pylab.axis([-0.5,max(x)+counter*width+1*width,0,ymax*1.3])
        
        
##             # make legend
##             #pylab.legend()
        
##         # increase plotCount
##         pylab.xlabel(xlabel)
##         if drewPlot:
##             plotCount += 1
##     pylab.subplots_adjust(left=0.125,right=0.70,hspace=0.4)
    
##     prop = FontProperties( size="xx-small" )

        
        
#pylab.plot([1,2,3],[1,2,3])
# save picture to cString.StringIO as a png then read to the file
f = cStringIO.StringIO()
pylab.savefig(f,format='png')
f.seek(0)
data = f.read()
print "Content-Type: image/png\nContent-Length: %d\n" % len(data)
print data




# arguments to this file will be pid, exp_id
# all calculation of which plots to make done in dynamics

# later will be option to do stuff to the data
# right now just returns the data as is

