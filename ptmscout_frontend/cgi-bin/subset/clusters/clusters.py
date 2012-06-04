import os
from pathname import *
class ClusterSet:
    def __init__(self,cluster_file_path):
        self.inputfile = cluster_file_path
#        inputfile = clusterPath+'exp%s.txt' % str(expid)
        self.parseError = None
        try:
            f = open(self.inputfile,'r')
            #f = open('clusters/MOAA_v1_Jan21_2009.txt','r')
            text = f.read()
            if '\r' not in text: text = text.replace('\n','\r')
            lines = text.split('\r')
            f.close()
        except IOError:
            lines = []
            self.parseError = "Error parsing file"
        self.indices = self.getClusterLineIndices(lines)
        self.names = self.getClusterNames(lines)
        #print self.indices, self.names
        self.clusterSets = self.getClusterSets(lines) ## clusterSets[name] = cluster:name
        if self.clusterSets == {}: self.parseError = "No clusters"
        #print self.clusterSets
    def getClusterLineIndices(self,lines):
        indices = []
        try:
            x = lines[0].split('\r')[0]
            x = x.split('\t')
            for i in range(len(x)):
                if 'cluster' in x[i] or 'Cluster' in x[i]:
                    indices.append(i)
        except IndexError: self.parseError = "Incorrect file format"
        return indices
        
    def getClusterNames(self,lines):
        try:
            clusters = lines[0].split('\r')[0]
            clusters = clusters.split('\t')[1:]
            clusters = [item[8:] for item in clusters if ('cluster' in item or 'Cluster' in item) and item !='']
        except IndexError:
            self.parseError = "Incorrect file format"
            return range(self.indices)
        return clusters

    def getClusterSets(self,lines):
        clusterSets = {}
        ##initialize dictionary
        for item in self.names:
            clusterSets[item]={}
        if len(lines)==1:
            newlines = lines[0].split('\r')[1:]
        else: newlines = lines[1:]
        for item in newlines:
            item = item.split('\t')
            msid = item[0]
            
            for i in range(len(self.indices)):
                index = self.indices[i]
                if len(item)> index:
                    #print i,self.names, self.indices, len(item)
                    if clusterSets[self.names[i]].get(item[index].replace('\r',''),None)==None and item[index]!='':
                        clusterSets[self.names[i]][item[index].replace('\r','')] = [msid]
                    elif item[index]!='':
                        clusterSets[self.names[i]][item[index].replace('\r','')].append(msid)
        return clusterSets

    def getK(self):
        return len(self.clusterSets.keys())


        
def checkFile(file):
    if 'cluster' not in file and 'Cluster' not in file:
        return False
    
    else:
        return True

def hackedFileName(exp_id, file_key):
    return clusterPath+'exp_%d_%d.txt'%(int(exp_id),int(file_key))
        
## x = ClusterSet(12)
## print x.clusterSets
## print x.parseError
