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
from pathname import *
sys.path.append(path+'includes')
sys.path.append(path)
import pylab
from matplotlib.patches import Ellipse,Circle
form = cgi.FieldStorage()
exp1 = form.getvalue("exp1",0)
exp2 = form.getvalue("exp2",0)
intersect = form.getvalue('intersect',0)
ambiguous = form.getvalue('ambiguous',0)
name = form.getvalue('title','')
name2 = form.getvalue('title2','')
comparison = int(form.getvalue('number',1))
type = 'site'


ells = [Circle((8,8),radius = 4),
        Circle((12,8),radius = 4),
        Ellipse((12,10.5),width=4,height=3, zorder=0)]
colors={}

colors['site']=[[ 0.77650468 , 0.19206367 , 0.44102354],[ 0.95240583,  0.03944834 , 0.96916472],"red"]

fig = pylab.figure()
ax=fig.add_subplot(111,aspect='equal')

titles = {1:'Phospho.Elm',2:'Phosphosite',3: 'Literature Only Phosphosite'}
title = name+'                 '+name2
counter = 0
for e in ells[0:2]:
    ax.add_artist(e)
    x = colors[type][counter]
    e.set_clip_box(ax.bbox)
    y= pylab.rand()+0.1
    e.set_alpha((counter+0.6)/5.0)
    e.set_facecolor(x)
    counter+=1

if int(ambiguous)>0:
    e = ells[2]
    ax.add_artist(e)
    x = colors[type][counter]
    e.set_clip_box(ax.bbox)
    e.set_alpha(0.5)
    e.set_facecolor(x)
    pylab.text(11,11,"Ambiguous")
    pylab.text(11.5,10.5,ambiguous)
pylab.text(5.5,8,exp2)
pylab.text(12.5,8,exp1)
pylab.text(9,8,intersect)

ax.set_xlim(4, 16)
ax.set_ylim(4, 12)
pylab.title(title,fontsize="x-small")
pylab.axis('off')
# save picture to cString.StringIO as a png then read to the file
f = cStringIO.StringIO()
pylab.savefig(f,format='png')
f.seek(0)
data = f.read()
print "Content-Type: image/png\nContent-Length: %d\n" % len(data)
print data
