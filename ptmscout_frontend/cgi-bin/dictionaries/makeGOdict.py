import pickle
import sys
try:
    GOfile = sys.argv[1]
except IndexError:
    print "Usage: python makeGODict.py <GOfile.obo>"
    sys.exit(0)
f = open(GOfile,'r')
lines = f.readlines()
f.close()
d= {}
for i in range(len(lines)):
    if lines[i][0:3]=="id:":
        key=lines[i].replace('id: ','').strip()
        val = lines[i+1].split(':')[1].strip()
        d[key]=val

f=open("GOdict.dict.pkl","w")
pickle.dump(d,f)
f.close()
