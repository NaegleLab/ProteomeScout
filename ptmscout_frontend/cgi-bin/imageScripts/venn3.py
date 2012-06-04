#!/usr/local/bin/python
from pathname import *
import os,sys
import cgi
import cgitb
if displayPythonErrors:
    cgitb.enable()
import cStringIO

os.environ['MPLCONFIGDIR']='%s'%mpl
import sys
sys.path.append(path+'includes')
sys.path.append(path)
import pylab
from matplotlib.patches import Ellipse,Circle
form = cgi.FieldStorage()
intall = form.getvalue("intall",'0')
int12 = form.getvalue("int12",'0')
int13 = form.getvalue("int13",'0')
int23 = form.getvalue("int23",'0')
only1 = form.getvalue("only1",'0')
only2 = form.getvalue("only2",'0')
only3 = form.getvalue("only3",'0')
title1 = form.getvalue("title1",'Experiment 1')
title2 = form.getvalue("title2",'Experiment 2')
title3 = form.getvalue("title3",'Experiment 3')
name1 = form.getvalue("name1","")
name2 = form.getvalue("name2","")
name3 = form.getvalue("name3","")
exp1 = form.getvalue("exp1",0)
exp2 = form.getvalue("exp2",0)
intersect = form.getvalue('intersect',0)
ambiguous = form.getvalue('ambiguous',0)
name = form.getvalue('title','')

comparison = int(form.getvalue('number',1))
type = 'site'

ells = [Circle((8,12),radius=4),Circle((12,12),radius=4),
        Circle((10,8),radius=4)]

colors={}

colors['site']=[[ 0.77650468 , 0.19206367 , 0.44102354],[ 0.95240583,  0.03944834 , 0.96916472],"red"]

fig = pylab.figure()
ax=fig.add_subplot(111,aspect='equal')

titles = {1:'Phospho.Elm',2:'Phosphosite',3: 'Literature Only Phosphosite'}
title = name2+'                 '+name
counter = 0
for e in ells:
    ax.add_artist(e)
    x = colors[type][counter]
    e.set_clip_box(ax.bbox)
    y= pylab.rand()+0.1
    e.set_alpha((counter+0.6)/5.0)
    e.set_facecolor(x)
    counter+=1




ax.set_xlim(4, 16)
ax.set_ylim(0, 16)

pylab.axis('off')
pylab.text(5.5,12,only1)#exp1only
pylab.text(12.5,12,only2)#exp2only
pylab.text(9.5,7,only3)#exp3only
pylab.text(9.5,10,intall)#intall
pylab.text(9.5,12.5,int12)#int12
pylab.text(7.3,9.4,int13)#int13
pylab.text(11.7,9.4,int23)#int23
pylab.text(3,16.1,title1)
pylab.text(11.5,16.1,title2)
pylab.text(6,3,title3)
#pylab.title(form.getvalue("name2"),fontsize="xx-small")
pylab.text(-2,-1,"Experiment 1: %s\nExperiment 2: %s\nExperiment 3: %s\n"%(name1, name2, name3),fontsize="x-small")
# save picture to cString.StringIO as a png then read to the file
f = cStringIO.StringIO()
pylab.savefig(f,format='png')
f.seek(0)
data = f.read()
print "Content-Type: image/png\nContent-Length: %d\n" % len(data)
print data

