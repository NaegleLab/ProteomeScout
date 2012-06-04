#!/usr/local/bin/python
from pathname import *
import os,sys
import cgi
import cgitb
if displayPythonErrors:
	cgitb.enable()
import cStringIO

sys.path.append(path+'load')
os.environ['MPLCONFIGDIR']='%s'%mpl
from prevFcns import *
import matplotlib
import pylab
import math
import random
from sets import Set
from matplotlib.font_manager import *
from matplotlib.text import *
## get form data
form = cgi.FieldStorage()
chart = form.getvalue("chart","line")
group = form.getvalue("group","run")
exp_id=int(form.getvalue("exp_id",52))
pid = int(form.getvalue("pid",9337))
msids = form.getvalue("msids",'')
errors = form.getvalue("errors",'true')
transform = form.getvalue("transform",None)
form = cgi.FieldStorage()
filename = form.getvalue("file","")
try:
	text = open(filename,'r').read()
except IOError: text = ''
text = text.replace('\r','\n')
dcols,rcol=getColumns(text)
types,labels =  getTypes(dcols,rcol,text)
testPep =  previewPeptide(text,pepColumn(text))
pepLines = getTestLines(testPep,text)
dataDict = getDataToGraph(pepLines,types,labels,rcol,dcols)
x = getType(text) # get type of graphs to make
if len(x)>1:
    if len(x[1])>1:
        if len(x[1])>0:
            if 'time' in x[1][0]:
                chart = 'line'
            else:
                chart = 'bar'
try:
    xlabel= x[1][0]
    if xlabel == 'average': xlabel = ''
except IndexError: xlabel = ''

runs = dataDict.keys()
if chart == "bar" and len(msids)>5:
    msids = msids[:5]

# calculate the number of graphs to draw
if group == 'run':
    numgraphs = len(x[0]) # get number of graphs to show, one for each run
elif group == 'site':
    numgraphs = len(msids) # get the number of graphs to show, one for each site
else:
    numgraphs = 1 # get the number of graphs to show, only one with all data


dataList = []
msids=[]
legend = []
dataList.append(dataDict)
print dataList
# data list has one entry for each phosphopep
if len(dataList)==1 and len(x[0])==1: group = 'run' # don't allow for group by site if there is only 1 site
title=''
if len(msids)==1:
    c.execute("""select acc_gene from protein join MS on MS.protein_id=protein.id where MS.id=%s""",(msids[0][0],))
    p1=c.fetchall()
    title=p1[0][0]


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
if len(dataList)>1: numCols = 2
else: numCols = 1
numRows = math.ceil(numgraphs*1.0/numCols)
##
#if len(msids)>1:
 #   width = 1.0/len(msids)

#else:
    #width = 1.0
width=0.5

# make dictionary of colors for msids
msColors = {}
runColors = {}
for item in msids:
    msColors[item[1]] =(random.random(),random.random(),random.random())
for item in runs:
    runColors[item] = (random.random(),random.random(),random.random())

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

if group == 'site': # make one graph per site with a line/bar for each run
    for d in dataList:
        drewPlot = False
        if numgraphs > 1:
            pylab.subplot(numRows,numCols, plotCount) # if there is only one, defaults to subplot(111)
        else:
            pylab.figure(figsize=(8,8))
            #pylab.subplot(1,1,1)
            
        counter = 0
        for run in runs:
            dictionary = d[run]
            if errors == 'false' or len(dictionary['yerrors'])==0:
                yerr = None
                leny = 0
            else:
                yerrs=dictionary['yerrors']
                yerr=transformStdDev(yerrs,dictionary['yvector'],transform)
               ##  yerr=[]
##                 for item in yerrs:
##                     if item ==None: ## need to deal with missing data point
##                                     ## if bar chart, can just be 0
##                                     ## else:
##                         yerr.append(0)
##                     else:
##                         yerr.append(item)
                leny=len(yerr)
            yvector = transformData(dictionary['yvector'],transform)
            #yvector = transformData(dictionary['yvector'])  # only plot if yvec isn't all 0s
            if dictionary['xvector']!=None:
                if chart == 'line' and not allZeros(yvector):
                    pylab.errorbar(dictionary['xvector'],yvector,yerr=yerr,label=run)
                    drewPlot = True
                else:
                    if not allZeros(yvector):
                        x=dictionary['xvector']
                        pylab.bar([(item+counter*width) for item in x],yvector,width,yerr=yerr,label=run,color=runColors[run])
                        drewPlot = True
            else:
                x = range(len(yvector))
                if chart == 'line' and not allZeros(yvector):
                    pylab.errorbar(x,yvector,yerr=yerr,label=run)
                    drewPlot = True
                else:
                    if not allZeros(yvector):
                        pylab.bar([(item+counter*width) for item in x],yvector,width,yerr=yerr,label=run,color=runColors[run])
                        drewPlot = True
            if not allZeros(yvector): # only increase counter if we graphed something!
                counter =counter+ 1
        # get title and xlabel, just take it from the first item in the run    
        pylab.title(title)#(msids[dataList.index(d)][1])
        if dictionary['xvector']!=None:
            #xlabels = dataList[0].values()[0].get('xvector',pylab.arange(8)),tuple(dataList[0].values()[0].get('xlabels',''))
            #pylab.xticks(range(len(xlabels)))
            pylab.xticks(dataList[0].values()[0].get('xvector',pylab.arange(8)),tuple(dataList[0].values()[0].get('xlabels','')),fontsize='x-small')
            
        else:
            pylab.xticks(pylab.arange(len(yvector)),tuple(dictionary['xlabels']),fontsize='xx-small')
            #pylab.axis([-0.5,max(x),0,max(yvector)*2+1])
        if drewPlot: plotCount += 1
        pylab.xlabel(xlabel)
        ## deal with missing data points:
        ## dictionary['xvector'],
       ##  y = dictionary.get('yvector',[])
##         x = dictionary.get('xvector',[0])
##         if len(y)>0:
##             for i in range(len(y)):
##                 if dictionary['yvector'][i]==None:
##                     if i-1 in range(len(x)):
##                         leftCoord = x[i-1]
##                     else:
##                         leftCoord = 0
##                     if i+1 in range(len(x)):
##                         rightCoord = x[i+1]
##                     else:
##                         rightCoord = max(x)
##                     topCoord = max(y)
##                     pylab.text(leftCoord,0,'',bbox=dict(facecolor='white', alpha=0))
                    
        
    if numgraphs>1:pylab.subplots_adjust(left=0.125,right=0.90,hspace=0.5,top=top)
    else:
        pylab.subplots_adjust(hspace=0.2,left=0.2,right=0.9,top=0.7 ,bottom=0.15)
    pylab.figure.figsize=(6,15)
    


    
if len(runs)>1: numCols = 2
else: numCols = 1
numRows = math.ceil(numgraphs*1.0/numCols)
            
# if group is run or all
maxSoFar = 0
if group == 'run' or group == 'all':
    dictionary = {}
    for run in runs:
        drewPlot = False
        if numgraphs > 1:
            pylab.subplot(numRows,numCols, plotCount) # if there is only one, defaults to subplot(111)
        # do everything
        counter = 0
        for d in dataList:
            dictionary = d[run]
            if errors == 'false' or len(dictionary.get('yerrors',0))==0:
                yerr = []
                leny = 0
            else:
                yerrs=dictionary['yerrors']
                yerrs = [float(item) for item in yerrs]
                yerr=transformStdDev(yerrs,dictionary['yvector'],transform)
               ##  yerr=[]
##                 for item in yerrs:
##                     if item ==None: ## need to deal with missing data point
##                         yerr.append(0)
##                     else:
##                         yerr.append(item)
                leny=len(yerr)
            for i in range(len(yerr)):
                if yerr[i]==None: yerr[i]=0
            yvector = transformData(dictionary.get('yvector',[]),transform)
            for i in range(len(yvector)):
                if yvector[i] == None: yvector[i] = 0 # interpolate data here
            maxSoFar = max(maxSoFar,max(yvector))
            if dictionary.get('xvector',None)!=None:
                if chart == 'line' and not allZeros(yvector):
                    pylab.errorbar(dictionary['xvector'],yvector,yerr=yerr)
                    drewPlot = True
                else:
                    if not allZeros(yvector):
                        x=dictionary['xvector']
                        pylab.bar([((item)+int(counter*width)) for item in x],yvector,width,yerr=yerr)
                        drewPlot = True
            else:
                x = range(len(yvector))
                if chart == 'line' and not allZeros(yvector):
                    pylab.errorbar(x,yvector,yerr=yerr)
                    drewPlot = True
                else:
                    if not allZeros(yvector):
                        pylab.bar([(item-int(counter*width)*0.3) for item in x],yvector,width,yerr=yerr)
                        drewPlot = True
            if not allZeros(yvector): # only increment if we actually drew something!
                counter =counter+ 1
                                # get title and xlabel, just take it from the first item in the run
        if title != '' and len(dataList)>0:
           pylab.title(title + ' - '+dataList[0][run].get('title',''))
        elif len(dataList)>0:
            pylab.title(dataList[0][run].get('title',''))
        if transform == "normToMax":
            ymax = 1.01
        else:
            ymax = maxSoFar
        if dictionary.get('xvector',None)!=None:
            #xlabels = dataList[0].values()[0].get('xvector',pylab.arange(8)),tuple(dataList[0].values()[0].get('xlabels',''))
            #pylab.xticks(range(len(xlabels)))
            pylab.xticks(dataList[0].values()[0].get('xvector',pylab.arange(8)),tuple(dataList[0].values()[0].get('xlabels','')))
            if chart=='bar':
                pylab.axis([-0.5,max(dictionary['xvector'])+counter*width+1*width,0,ymax])
            else:
                pylab.axis([-0.5,max(dictionary['xvector']),0,ymax*1.3])
                
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
        
        
            # make legend
            #pylab.legend()
        
        # increase plotCount
        pylab.xlabel(xlabel)
        if drewPlot:
            plotCount += 1
    if numgraphs>1: pylab.subplots_adjust(left=0.125,right=0.9,hspace=0.4)
    
    prop = FontProperties( size="xx-small" )


    
        


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

