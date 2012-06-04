from pathname import *
import pylab
import cgi
import cgitb
if displayPythonErrors:
	 cgitb.enable()
import os
import sys

sys.path.append(path+'includes')
sys.path.append(path)
from matplotlib.patches import Ellipse, Circle
import cStringIO
os.environ['MPLCONFIGDIR']=mpl
import sys
from pathname import *

## get arguments
form=cgi.FieldStorage()
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
##
ells = [Circle((8,12),radius=4),Circle((12,12),radius=4),
        Circle((10,8),radius=4)]

fig=pylab.figure()
ax=fig.add_subplot(111,aspect='equal')
counter=0
for e in ells:
    ax.add_artist(e)
    e.set_clip_box(ax.bbox)
    e.set_alpha((counter+0.6)/5.0)
    e.set_facecolor('red')
    counter+=1
ax.set_xlim(4,16)
ax.set_ylim(0,20)
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
pylab.text(-2,-1,"Experiment 1: %s\nExperiment 2: %s\nExperiment 3: %s\n"%(name1, name2, name3),fontsize="small")

f=cStringIO.StringIO()
pylab.savefig(f,format='png')
f.seek(0)
data = f.read()
print "Content-Type: image/png\nContent-Length: %d\n" % len(data)
print data
