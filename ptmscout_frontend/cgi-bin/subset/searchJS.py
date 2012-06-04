from pathname import *
import sys
sys.path.append(path+'includes')
from template import *
def pvalJS():
	print """
	<script type="text/javascript">
	function pValCorrect(table,radio,text,pvaltext)
	{
	
	// get threshold value
	threshold = parseFloat(document.getElementsByName(text)[0].value);

	// get correction type
	radio = document.getElementsByName(radio);
	type = "None";
	for (i =0; i< radio.length; i++)
	{
	   if (radio[i].checked)
	   {
              type = radio[i].value;
	   }
	}

	// go down the table, highlight until get to values above threshold
	table = document.getElementsByName(table)[0];
	m = table.rows.length;
	alpha = 0;
	beta = 0;
	numenriched = 0;
	if(type=="BF"){alpha = threshold/m;}
	else {alpha = -1}
	if(type=="RFDR"){
	beta = (threshold*(m+1))/(2*m);
	}
		else {beta=-1}
        
	if(type=="FDR"){
	    eta = -1
	   for (i=1; i < table.rows.length; i++)
	   {
	      // get numerical value of the row
	      num = table.rows[i].cells[3].innerHTML;
	      
	      num = num.replace('<b>','');
	      num = num.replace('<\/b>','');
	       num = num.replace('<b>','');
	      num = num.replace('<\/b>','');
	      
	      num = num.replace('<\/sup>','');
	      num = num.replace('<sup>','x');
	      num = num.replace('<font color="black">','');
	      num = num.replace('<font color="red">','');
	      num = num.replace('<font color="green">','');
	      num = num.replace('<\/font>','');
	      var base = parseFloat(num.split('x')[0]);
	      var power = parseFloat(num.split('x')[2]);
	      var pvalue = base*Math.pow(10,power);
              if (pvalue <((i+1)*threshold/(m))){eta = pvalue;}
	      
	   }
	}
	else {eta = -1}
       
	      for (i = 1; i < table.rows.length; i++)
	      {
	      
	      
	      
	      // get numerical value of the row
	      num = table.rows[i].cells[3].innerHTML;
	      
	      num = num.replace('<b>','');
	      num = num.replace('<\/b>','');
	       num = num.replace('<b>','');
	      num = num.replace('<\/b>','');
	      
	      num = num.replace('<\/sup>','');
	      num = num.replace('<sup>','x');
	      num = num.replace('<font color="black">','');
	      num = num.replace('<font color="red">','');
	      num = num.replace('<font color="green">','');
	      num = num.replace('<font color="orange">','');
	      num = num.replace('<\/font>','');
	      var base = parseFloat(num.split('x')[0]);
	      var power = parseFloat(num.split('x')[2]);
	      var pvalue = base*Math.pow(10,power);
              
	      // code for different types of correction here
	      thresh = threshold;
              
	      if (pvalue <= thresh)
	         
              {
	      numenriched ++;
                  if (num.search('<b>')==-1)
		    {
		 
		         for(j = 0; j < 4; j++)
			 {
			   
			    oldHTML = table.rows[i].cells[j].innerHTML;
		            oldHTML = '<b>'+oldHTML+'<\/b>';	    
			   table.rows[i].cells[j].innerHTML = oldHTML;
			 }
		    }
              }
              
              
              else
	      {
                 for(j=0; j<4;j++)
		 {
                    oldHTML = table.rows[i].cells[j].innerHTML;
		    oldHTML = oldHTML.replace('<b>','');
		    oldHTML = oldHTML.replace('<\/b>','');
		    oldHTML = oldHTML.replace('red','black');
		    oldHTML = oldHTML.replace('green','black');
		     oldHTML = oldHTML.replace('orange','black');
		    table.rows[i].cells[j].innerHTML = oldHTML;
		 }
	      }
	      if (pvalue <= alpha)
	         
              {
	           numenriched ++;
		         for(j = 0; j < 4; j++)
			 {
			   
			    oldHTML = table.rows[i].cells[j].innerHTML;
		            oldHTML = oldHTML.replace('black','red');
			    oldHTML = oldHTML.replace('green','red');
			    oldHTML = oldHTML.replace('orange','red');
			   table.rows[i].cells[j].innerHTML = oldHTML;
			 }
		    
              }
              
              
              else
	      {
                 for(j=0; j<4;j++)
		 {
                    oldHTML = table.rows[i].cells[j].innerHTML;
		    //oldHTML = oldHTML.replace('red','black');
		    table.rows[i].cells[j].innerHTML = oldHTML;
		 }
	      }
	      if (pvalue <= beta)
	         
              {
	          numenriched ++;
                 
		         for(j = 0; j < 4; j++)
			 {
			   
			    oldHTML = table.rows[i].cells[j].innerHTML;
		            oldHTML = oldHTML.replace('black','green');
			    oldHTML = oldHTML.replace('red','green');
			    oldHTML = oldHTML.replace('orange','green');
			   table.rows[i].cells[j].innerHTML = oldHTML;
			 }
		    
              }
              
              
              else
	      {
                 for(j=0; j<4;j++)
		 {
                    oldHTML = table.rows[i].cells[j].innerHTML;
		    //oldHTML = oldHTML.replace('green','black');
		    table.rows[i].cells[j].innerHTML = oldHTML;
		 }
	      }
	      if (pvalue <= eta)
	         
              {
                 numenriched ++;
		         for(j = 0; j < 4; j++)
			 {
			   
			    oldHTML = table.rows[i].cells[j].innerHTML;
		            oldHTML = oldHTML.replace('black','orange');
			    oldHTML = oldHTML.replace('red','orange');
			    oldHTML = oldHTML.replace('green','orange');
			   table.rows[i].cells[j].innerHTML = oldHTML;
			 }
		    
              }
              
              
              else
	      {
                 for(j=0; j<4;j++)
		 {
                    oldHTML = table.rows[i].cells[j].innerHTML;
		    //oldHTML = oldHTML.replace('green','black');
		    table.rows[i].cells[j].innerHTML = oldHTML;
		 }
	      }
	      
             
	      
              

	      }
	      if (numenriched == 0){
	      document.getElementsByName(pvaltext)[0].innerHTML = "No terms enriched.";
	      }
	      else{
	      document.getElementsByName(pvaltext)[0].innerHTML = "";
	      }
	   
	}
	</script>
	"""

def printJavaScript():
	print """
	<script type="text/javascript">
	function setParentHeight()
	{
	
	var alarm;
        var now = new Date();
        var sleeping = true;
        var startingMSeconds = now.getTime();
        while(sleeping){
            alarm = new Date();
            alarmMSeconds = alarm.getTime();
            if(alarmMSeconds - startinMSeconds > 200000)
            { sleeping = false;}
            }
	var height = document.height;
	   for (i = 1; i < 30; i++){
	   parent.document.getElementsByName('sub'+i)[0].style.height=height;
	   }
	    
        }
	function tog(o,button)
	{
	    
	    var e = document.getElementById(o);
	    e.style.display = e.style.display == 'block' ? 'none' : 'block';
	   
	    button = document.getElementsByName(button)[0]
	    if(button.name == 'GOButton')
	    {button.value = button.value == 'Show GO' ? "Hide GO":"Show GO" ;}
	    if(button.name == 'SCANButton')
	    {button.value = button.value == 'Show Predictions' ? "Hide Predictions":"Show Predictions" ;}
	    if(button.name == 'DOMButton')
	    {button.value = button.value == 'Show Domains' ? "Hide Domains":"Show Domains" ;}
	    if(button.name == 'PFAMButton')
	    {button.value = button.value == 'Show PFAM' ? "Hide PFAM":"Show PFAM" ;}
	    if(button.name == 'DYNButton')
	    {button.value = button.value == 'Show Data' ? "Hide Data":"Show Data" ;}

	    if(button.name == 'GO_BPButton'||button.name == 'GO_MFButton' ||button.name == 'GO_CCButton' ||button.name =='PFAMTableButton'||button.name == 'DOMTableButton'||button.name == 'SCANTableButton'||button.name == 'KINTableButton' || button.name == 'PELMTableButton')
	    {button.value = button.value == 'Show Table' ? "Hide Table":"Show Table";}
	    if(button.name == 'GO_BPPieButton'||button.name == 'GO_MFPieButton' ||button.name == 'GO_CCPieButton' || button.name == 'PFAMPieButton'||button.name == 'DOMPieButton' ||button.name == 'SCANPieButton'||button.name == 'KINPieButton' ||button.name == 'PELMPieButton')
	    {button.value = button.value == 'Show Graph' ? "Hide Graph":"Show Graph";}
	    
	    var height = document.height;
	   
	    parent.document.getElementsByName('sub1')[0].style.height=height;
	     parent.document.getElementsByName('sub2')[0].style.height=height;
	      parent.document.getElementsByName('sub3')[0].style.height=height;
             parent.document.getElementsByName('sub4')[0].style.height=height;
	      parent.document.getElementsByName('sub5')[0].style.height=height;
	       parent.document.getElementsByName('sub6')[0].style.height=height;
	       parent.document.getElementsByName('sub7')[0].style.height=height;
	     parent.document.getElementsByName('sub8')[0].style.height=height;
	      parent.document.getElementsByName('sub9')[0].style.height=height;
             parent.document.getElementsByName('sub10')[0].style.height=height;
	      parent.document.getElementsByName('sub11')[0].style.height=height;
	       parent.document.getElementsByName('sub12')[0].style.height=height;
	    
	}
	</script>
	"""
	print """
 	<script type="text/javascript" src="%schart/includes/excanvas.js"></script>
 	<script type="text/javascript" src="%schart/includes/chart.js"></script>
 	<script type="text/javascript" src="%schart/includes/canvaschartpainter.js"></script>
 	<link rel="stylesheet" type="text/css" href="%schart/includes/canvaschart.css" >"""%(jsPath,jsPath,jsPath,jsPath)
 	print """<link rel="stylesheet" type="text/css" href="%stest.css">"""%jsPath
 	print """<link rel="stylesheet" type="text/css" href="%ssubset_search.css">"""%jsPath
	
	print """
	<script type="text/javascript">
	function toggle(o) {
	var e = document.getElementById(o);
	e.style.display = e.style.display == 'block' ? 'none' : 'block';
	}
	function toggleGO(d1,d2,d3,form)
	{

	var c = document.getElementById(d1);
	var b = document.getElementById(d2);
	var f = document.getElementById(d3);  
	c.style.display = (c.style.display == 'none') ? 'block' : 'none';
	b.style.display = (b.style.display == 'none') ? 'block' : 'none';
	f.style.display = (f.style.display == 'none') ? 'block' : 'none'; 
	if (form.tog.value =="show GO terms")
	    {
	    form.tog.value = "hide GO terms";
	    }
        else
	    {
	    form.tog.value = "show GO terms";
	    }
        if (form.tog.value =="show PFAM terms")
	    {
	    form.tog.value = "hide PFAM terms";
	    }
        else
	    {
	    form.tog.value = "show PFAM terms";
	    }
        if (form.tog.value =="show Domain terms")
	    {
	    form.tog.value = "hide Domain terms";
	    }
        else
	    {
	    form.tog.value = "show Scansite terms";
	    }
        if (form.tog.value =="show Scansite terms")
	    {
	    form.tog.value = "hide Scansite terms";
	    }
        else
	    {
	    form.tog.value = "show Scansite terms";
	    }
        if (form.tog.value =="show Dynamics terms")
	    {
	    form.tog.value = "hide Dynamics terms";
	    }
        else
	    {
	    form.tog.value = "show Dynamics terms";
	    }
	}
	"""
	print '</script>'
	print """
	<script type="text/javascript">
	function showmenu(elmnt)
	{
	document.getElementById(elmnt).style.visibility="visible";
	}
	function hidemenu(elmnt)
	{
	document.getElementById(elmnt).style.visibility="hidden";
	}
	</script>
        """
