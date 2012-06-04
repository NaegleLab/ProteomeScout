#!/usr/local/bin/python
from pathname import *
print "Content-Type: text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">"""
import sys


sys.path.append(path+'includes')
from template import *
form = cgi.FieldStorage()
printHeader("PTMScout - About")
print "<BR>"*2
print """

<div style="text-align:left;width:65%">
<p><font color="#0000A0"><h2>Attribution Reference</h2></font>
Please cite the following publication when publishing results derived from use of PTMScout:
<ul>
<li>Naegle KM, Gymrek M, Joughin BA, Wagner JP, Welsch RE, Yaffe MB, Lauffenburger DA, White FM. <i>PTMScout: A web resource for analysis of high-throughput post-translational proteomic studies.</i> Molecular and Cellular Proteomics, 2010. 
</ul>
Please cite the following publication when publishing results derived from use of MCAM:
<ul>
<li>Naegle KM, Welsch RE, Yaffe MB, White FM, Lauffenburger DA.  <i>MCAM: Multiple Clustering Analysis Methodology for Deriving Hypotheses and Insights from High-Throughput Proteomic Datasets.</i>  PLoS Comput Biol 7(7), 2011.
</ul>
</p>
<BR>
<font color="#0000A0"><h2>License, Disclaimer and Terms of Use</h2></font>
<h3><a href="http://creativecommons.org/licenses/by/3.0/legalcode">Licensed under Creative Commons Attribution 3.0</a></h3>

<font color="#0000A0"><h3>Disclaimer</h3></font>
<font size="-1"><p>The PTMScout web site (the "Site") is an online, open-content, collaborative project that represents a voluntary association of individuals and groups working to develop a common, educational resource. Use of the site is completely voluntary; and anyone may access the data presented on the Site and/or upload data to the Site. The Site's content is not uniformly peer reviewed, and, therefore, all of the Site's information is provided without representations or warranties of any kind, express or implied, including, without limitation, warranties of merchantability, fitness for a particular purpose, or noninfringement. None of the administrators, contributors, or anyone else associated with the Site is responsible for the appearance of any inaccurate or libelous information, your use of the information contained on this Site, or modifications made to Site content. The Site's database is stored on a server in the Commonwealth of Massachusetts in the United States of America, and is therefore maintained in reference to the protections afforded under local and federal law. Laws in your country or jurisdiction may not protect or allow the same kind(s) of speech or distribution. The Site does not encourage the violation of any laws, and, therefore, cannot assume responsibility for any violations of such laws attributable to your linking to this domain or using, reproducing, or publishing material contained herein.
<br><br>
The Site's information is furnished freely, and you are specifically granted a limited license to copy anything from this Site pursuant to a <a href="http://creativecommons.org/licenses/by/3.0/legalcode">"Creative Commons License"</a>. Notwithstanding the benefits provided by this limited license, no agreement or contract is created between you and the Site\'s owners or users; the owners of the servers on which the Site is housed; the individual Site contributors; any Site administrators, system operators; or anyone else connected (in any manner) to this project or sister projects subject to your claims directly. Furthermore, this license neither creates nor implies any contractual relationship or liability on the part of M.I.T., its trustees, directors, officers, employees, agents, or associates. Any trademarks, service marks, collective marks, design rights or similar rights that are mentioned, used, or cited on the Site remain the property of their respective owners. Unless otherwise specifically provided, the Site is not endorsed by or affiliated with the holders of any such marks/rights. Any use of the Site's material is completely at your own risk. Thus, users should not rely on any information contained in the Site without independent verification. The Site is a work in progress, so please notify Site administrators of any errors discovered on the Site.
</p></font>

<font color="#0000A0"><h3><a href="http://ptmscout.mit.edu/cgi-bin/terms.cgi">Terms of Use</a></h3></font>

<hr size=10 color="blue">
<BR>
<h2><font color="#0000A0">Developers</h2></font>
<ul>
<li>Kristen Naegle
<li>Melissa Gymrek
</ul>
<BR>
<h2><font color="#0000A0">Contact Us</h2></font>
<ul>
<li><a href="mailto:ptmscout_admin@mit.edu">ptmscout_admin@mit.edu</a>
</ul>


<h2><font color="#0000A0">Download PTMScout</h2></font>
<ul>
<li>In order to download your own version of PTMScout software, for analysis of non-public data, please e-mail the administrator at <a href="mailto:ptmscout_admin@mit.edu">ptmscout_admin@mit.edu</a>.  The software is distributed under the <a href=http://creativecommons.org/licenses/by-nc-sa/3.0/>Creative Commons Attribution Noncommercial license.</a>
</ul>

<h2><font color="#0000A0">Download MCAM Matlab Software</h2></font>
<ul>
<li><a href="/MCAM_July_2011.zip">MCAM Matlab Code</a>
</ul>

<BR>
<font color="#0000A0"><h2>Support</h2></font>
This work was supported by:
<ul>
<li>NIH-U54-CA112967
<li>NIH-RO1-CA096504
</ul>
<BR>
<font color="#0000A0"><h2>Privacy</h2></font>
<ul>
<li>Cookies are used on this site only for the collection of anonymous information.  The goal of this information is to inform the developers of PTMScout of the usefulness of the tools provided and to aid in the pursuit of funding for maintenance and further development of PTMScout.
</ul>
</div>
"""
print "<BR>"*3

printFooter()
