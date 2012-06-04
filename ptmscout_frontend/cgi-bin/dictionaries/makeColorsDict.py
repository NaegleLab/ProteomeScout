from pathname import *
import random

import sys
sys.path.append(path)
import MySQLdb
db=MySQLdb.connect(user="mgymrek",passwd="melissa_mysql",db="webdev")
c=db.cursor()
import pickle

f = open('%s/dictionaries/colors.dict.pkl'%path,'r')
d=pickle.load(f)
f.close()

query = """select label from domain"""
c.execute(query)
x=c.fetchall()
domains = [item[0] for item in x]

def tohex(r,g,b):
    hexchars = "0123456789ABCDEF"
    return "#" + hexchars[r / 16] + hexchars[r % 16] + hexchars[g / 16] + hexchars[g % 16] + hexchars[b / 16] + hexchars[b % 16]


for item in domains:
    if item not in d.keys():
        key = item
        val = (random.randint(0,255),random.randint(0,255), random.randint(0,255))
        val = tohex(val[0],val[1],val[2])
        d[key]=val
f = open('%s/dictionaries/colors.dict.pkl'%path,'w')
pickle.dump(d,f)
f.close()

