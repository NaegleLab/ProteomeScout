#!/usr/local/bin/python
from pathname import *
import os,sys
import cgi
import cgitb
if displayPythonErrors:
	cgitb.enable()
import cStringIO

os.environ['MPLCONFIGDIR']='%s'%mpl
from dynamics import *
import matplotlib
import pylab
import math
import random
## get form data
form = cgi.FieldStorage()
# will take labels, fractions, title
def sortListByElement(list,index):
	
	tmp = [(item[index],item) for item in list]
	
	tmp.sort()
	
	list = [t[1] for t in tmp]
	return list

fractions='0.03125:0.03125:0.03125:0.03125:0.03125:0.03125:0.0625:0.0625:0.0625:0.0625:0.1875:0.1875:0.03125:0.5625:0.09375:0.03125:0.21875:0.125:0.1875:0.125:0.1875:0.03125:0.0625:0.0625:0.09375:0.1875:0.03125:0.03125:0.0625:0.03125:0.03125:0.1875:0.5625:0.1875:0.1875:0.03125:0.0625:0.0625:0.03125:0.03125:0.1875:0.09375:0.0625:0.15625:0.1875:0.1875:0.09375:0.4375:0.03125:0.03125:0.0625:0.1875:0.25:0.03125:0.03125:0.09375:0.1875:0.375:0.03125:0.03125:0.03125:0.1875:0.25:0.09375:0.03125:0.09375:0.0625:0.25:0.0625:0.0625:0.21875:0.03125:0.03125:0.1875:0.03125:0.03125:0.03125:0.03125:0.03125:0.03125:0.03125:0.1875:0.09375:0.03125:0.03125'
labels='GO:0008284:GO:0008104:GO:0016337:GO:0031575:GO:0045930:GO:0030032:GO:0030154:GO:0009887:GO:0018108:GO:0006928:GO:0007155:GO:0048008:GO:0007163:GO:0007165:GO:0035023:GO:0007169:GO:0008219:GO:0009653:GO:0042981:GO:0007275:GO:0051301:GO:0009719:GO:0007389:GO:0007254:GO:0001701:GO:0016049:GO:0030838:GO:0051209:GO:0001568:GO:0006919:GO:0042531:GO:0016043:GO:0050789:GO:0030335:GO:0050852:GO:0007049:GO:0048538:GO:0006886:GO:0019725:GO:0016477:GO:0007186:GO:0016042:GO:0060017:GO:0008285:GO:0006468:GO:0050853:GO:0006629:GO:0008283:GO:0030218:GO:0006917:GO:0006810:GO:0008286:GO:0006464:GO:0042102:GO:0006882:GO:0015031:GO:0001558:GO:0008150:GO:0042110:GO:0051249:GO:0030097:GO:0048011:GO:0007010:GO:0007265:GO:0050862:GO:0009056:GO:0007507:GO:0007015:GO:0006897:GO:0009952:GO:0007242:GO:0019538:GO:0009725:GO:0007229:GO:0042493:GO:0030866:GO:0001935:GO:0007154:GO:0050870:GO:0030217:GO:0006417:GO:0007173:GO:0009790:GO:0006412:GO:0007172'
title='GO Biologcial Processes'


fractions = form.getvalue('fractions',fractions)
labels = form.getvalue('labels',labels)
title = form.getvalue('title',title)
type = form.getvalue('type','GO')
fractions = fractions.split(':')
fractions = [float(item)*100 for item in fractions]
labels= labels.split(':')
subtext = ''#'\nPercentages are relative to the chosen subset. \nView the help page to read more about \nhow percentages are calculated.'
if type == 'GO':
    newlabels= []
    for i in range(0,len(labels),2):
        newlabels.append(labels[i]+':'+labels[i+1])
    labels = newlabels
    for i in range(len(labels)):
        if labels[i] == 'null:':
            labels[i] = 'None'
if type == 'pfam' or type =='dom' or type in ['scan','kinase','pelm','bind']:
    for i in range(len(labels)):
        if labels[i] == "~~~":
            labels[i] = 'None'
  
### order from largest to greatest
# first make a list of tuples of (value,label)
vals = [(fractions[i],labels[i]) for i in range(len(fractions))]
newvals = sortListByElement(vals,0)
#newvals.reverse()
fractions = [item[0] for item in newvals]
labels = [item[1] for item in newvals]
### end ordering


### get rid of super small entries.. for now only graph the first 10
index = len(fractions)-1
set = False
for i in range(len(fractions)):
	if fractions[i] < 3 and set == False and i > 2:
		set = True
		index = i
		
if len(fractions)>15:
    fractions = fractions[:index] + [sum(fractions[index:])]
    labels = labels[:index]+['Other']
if len(fractions)>10:
	fractions = fractions[:10] + [sum(fractions[10:])]
	labels = labels[:10]+['Other']
### end getting rid of super small enries
if len(labels)>2:
	pylab.figure(1, figsize=(14,14))
	pylab.pie(fractions,labels=labels,shadow=True, autopct="%0.2f",colors=('b','g', 'r', 'c', 'm', 'y', 'w'))


## make legend
import sys
from pathname import *
sys.path.append(path)
f = open(path+'dictionaries/GOdict.dict.pkl','rb')
import pickle
dict = pickle.load(f)
f.close()
text = ''

for item in labels:
	text=  item + ' - ' + dict.get(item,'')+'\n'+text
if 'GO:' in labels[0]:
	pylab.text(-0.8,-1.6,text,size=12)

##
if len(labels)>2:
	pylab.title(title)
	pylab.text(-1.5,1.3,subtext)

# save picture to cString.StringIO as a png then read to the file
f = cStringIO.StringIO()
pylab.savefig(f,format='png')
f.seek(0)
data = f.read()
print "Content-Type: image/png\nContent-Length: %d\n" % len(data)
print data


