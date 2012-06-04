def clusterJS(formdata,names): # have all but cluster 0, which will be set in
    # each of the tabs
    K = len(names)
    if False not in [item.isalnum() and not(item.isalpha()) for item in names]:
        names = [int(item) for item in names]
    names.sort()
    print"""
<script type = "text/javascript">"""
    print "var count = parent.document.getElementsByName('counter')[0].value-1;"
    #print "alert(count+2);"
    for i in range(1,K+1):
        #print "alert('sub'+'%s'+1)"%i
        if i != K:
            print "parent.set();"
        print "parent.document.getElementsByName('sub'+(%s+count))[0].contentWindow.document.location = 'search.cgi?%s&cluster0=%s';" %(i,formdata,names[i-1])
        
        print """parent.document.getElementById('font'+(%s+count)).innerHTML =parent.document.getElementById('font'+(count+%s)).innerHTML.replace("Subset","Cluster");"""% (i, i)
        
        print """var rep3 = parent.document.getElementById('font'+(%s+count)).innerHTML.split('Cluster')[0];"""% (i)
        print """var rep4 = parent.document.getElementById('font'+(%s+count)).innerHTML.split('Cluster')[1].split('</font>')[1];"""% (i)
        print """var rep5 = rep3 + "Cluster "+%s + "</font>"+rep4;"""%i
        
        print """parent.document.getElementById('font'+(%s+count)).innerHTML=rep5;"""%i
        
        #print """parent.document.getElementById('font'+(%s+count)).innerHTML=parent.document.getElementById('font'+(%s+count)).innerHTML.replace(rep,new)"""
        
        print "parent.document.getElementById('font'+(%s+count)).display = 'block';"%i
        print "parent.document.getElementsByName('sub'+(%s+count))[0].display = 'block';"%i
   
    
    print"""
</script>
"""
