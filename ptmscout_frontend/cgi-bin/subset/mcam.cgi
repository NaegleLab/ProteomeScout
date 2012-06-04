#!/usr/local/bin/python
from pathname import *
import sets
import cgi
import math

import cgitb
if displayPythonErrors:
    cgitb.enable()

print "Content-Type: text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">"""
import MySQLdb
import sys

sys.path.append(path+'includes')
sys.path.append(path)
sys.path.append(path+'clusters')
from template import *
import pickle
import os



from pathname import *
sys.path.append('clusters')

from clusters import *

# these are variables that come from the http post back (the user clicks submit on this page)
motif_cutoff=""
exp_id=""
file_key=""
scansite_cutoff=""
domain_cutoff=""
alpha_p=""
email=""
email_confirm=""
correction=""

# this is special to the form on this page -- it tells us if we should read the postback data from
# the form and try to validate and process the user's inputs, or if we should draw the page
# with the initial values.  An empty value means draw the page fresh.
action="" 
user_file_name = ""

# these are error messages specifc to each field, set if validation detects an error
# if you add one, you'll have to touch the validation function and the has_errors function
# We call the perl script if has_errors returns false
motif_cutoff_error=""
exp_id_error=""
file_key_error=""
scansite_cutoff_error=""
domain_cutoff_error=""
alpha_p_error=""
email_error=""
email_confirm_error=""
correction_error=""

def has_errors():
    return(motif_cutoff_error    !="" or
           exp_id_error          !="" or
           file_key_error        !="" or
           scansite_cutoff_error !="" or
           domain_cutoff_error   !="" or
           alpha_p_error         !="" or
           email_error           !="" or
           email_confirm_error   !="" or
           correction_error      !="")

# reads from the HTTP post into or local script variables
def read_postback_data():
    # ---
    # Get the form post information
    form = cgi.FieldStorage()

    # this form is posted two in two different cases
    # 1) The user is coming from select.cgi and has chosen to do MCAM
    # 2) The user has put in their mcam inputs
    #
    # The inputs to this page are:
    # - exp_id (integer)
    # - file_key (integer)
    # - email (string)
    # - email_confirm (string)
    # - motif_cutoff - value between 0 and 1, inclusive?
    # - scansite_cutoff - value between 0 and 5, inclusive?
    # - domain_cutoff - scientific notation
    # - alpha_p - value between 0 and 1, inclusive
    # - correction - either BH FDR or Bonferroni
    # - user_file_name - the name of the user's original file.  Do not ever use this to access a file on disk
    global exp_id, file_key, email, email_confirm, motif_cutoff, scansite_cutoff, domain_cutoff, alpha_p, correction, action, user_file_name
    
    exp_id          = form.getvalue("exp_id", "")
    file_key        = form.getvalue("file_key", "")
    email           = form.getvalue("email", "")
    email_confirm   = form.getvalue("email_confirm", "")
    motif_cutoff    = form.getvalue("motif_cutoff", "")
    scansite_cutoff = form.getvalue("scansite_cutoff", "")
    domain_cutoff   = form.getvalue("domain_cutoff", "")
    alpha_p         = form.getvalue("alpha_p", "")
    correction      = form.getvalue("correction", "")
    action          = form.getvalue("action", "")
    user_file_name  = form.getvalue("user_file_name", "")

def set_default_values():
    global motif_cutoff, scansite_cutoff, domain_cutoff, alpha_p, correction
    if motif_cutoff=="":
        motif_cutoff="0.01"

    if scansite_cutoff=="":
        scansite_cutoff="3"

    if alpha_p=="":
        alpha_p="0.05"

    if correction=="":
        correction = "BH"

    if domain_cutoff=="":
        domain_cutoff = "1e-5"
    
# Draws a tabl row with a validation error
def render_error(validation_error):
    if validation_error != "":
        print """
        <tr>
        <td>&nbsp;</td>
        <td colspan="2" class="validation_error">%s</td>
        </tr>"""%(str(validation_error))

# Draws a simple text field in our form table
def render_text_area(field_label, field_name, current_value, validation_error, field_size, help_text):
    print """
    <tr title="%s">
      <th><label for="%s">%s</label></th>
      <td><input type="text" name="%s" value="%s" size="%d"></td>
    </tr>"""%(str(help_text), str(field_name), str(field_label), str(field_name), str(current_value), int(field_size))
      
    render_error(validation_error)  
      
# Draws a radio button and a label, handling if it should be checked or not
def render_radio_button(button_label, field_name, value, checked):
    print """ <input type="radio" name="%s" value="%s" """%(str(field_name), str(value))
    if ( checked ):
        print " checked"
    print ">"
    print """<label for="%s">%s</label>"""%(str(field_name), button_label)

# Renders the table row in our form input for the correction radio buttons
def render_correction():
    help_text = ""
    print "<tr title='%s'><th>Correction</th><td>"%(str(help_text))
    render_radio_button("BH FDR", "correction", "BH", correction == "BH")
    render_radio_button("Bonferroni", "correction", "BF", correction != "BH")        
    print "</td>"
    print "</tr>"
    render_error(correction_error)
        
    
# renders the webpage for providing MCAM inputs
def render_form():
    if has_errors():
        print "<div class='error'>There were errors processing your request.</div>"
        
    print """
    <form method="POST" action="mcam.cgi">
    <input type="hidden" name="action" value="process" />
    <input type="hidden" name="exp_id" value="%s" />
    <input type="hidden" name="file_key" value="%s" />"""%(str(exp_id), str(file_key))
    print """<input type="hidden" name="user_file_name" value="%s" />"""%(str(user_file_name))
    
    print "<fieldset><legend>MCAM Inputs</legend><table>"

    render_text_area("Email", "email", email, email_error, 60, "Notification of completion will be sent to this email address.")
    render_text_area("Email Confirmation", "email_confirm", email_confirm, email_confirm_error, 60, "")
    render_text_area("Motif Cutoff", "motif_cutoff", motif_cutoff, motif_cutoff_error, 10, "This is some help")
    render_text_area("Scansite Cutoff", "scansite_cutoff", scansite_cutoff, scansite_cutoff_error, 10, "")
    render_text_area("Domain Cutoff", "domain_cutoff", domain_cutoff, domain_cutoff_error, 20, "")
    render_text_area("Alpha P", "alpha_p", alpha_p, alpha_p_error, 10, "")
    render_correction()
    
    print """
    </table>
    <br />
    <input type="submit" value="Submit">
    </fieldset>"""
    
    print "</form>"

# Validates a float falls in a certain range
def validate_float(value, field_name, min, max):
    if (value ==""):
        return (value, "%s is a required field"%(str(field_name)))

    try:
        x = float(value)
    except ValueError:
        return (value, "%s must be between %f and %f"%(str(field_name), min, max))
    
    # x != x tests for Not a Number.  the math.isnan would be better, but we don't have that because we use ghetto phyton
    
    if ( x != x or x < min or x > max ):
        return (value, "%s mucatch must be between %f and %f"%(str(field_name), min, max))

    return (x, "")

# validates all our form inputs before we call the perl script
def validate_inputs():
    # We'll modify these if they pass validation so they are numbers, not strings
    global motif_cutoff, scansite_cutoff, domain_cutoff, alpha_p, correction

    # these we'll set to signal general errors with the page
    global motif_cutoff_error, scansite_cutoff_error, domain_cutoff_error, alpha_p_error, correction_error
    global email_error, email_confirm_error
    if email=="":
        email_error = "Email is required"
    else:
        if email != email_confirm:
            email_confirm_error = "Email confirmation does not match"

    (motif_cutoff, motif_cutoff_error)       = validate_float(motif_cutoff, "Motif Cutoff", 0.0, 1.0)
    (scansite_cutoff, scansite_cutoff_error) = validate_float(scansite_cutoff, "Scansite Cutoff", 0.0, 5.0)
    (alpha_p, alpha_p_error)                 = validate_float(alpha_p, "Alpha P", 0.0, 1.0)

    
    (domain_cutoff, domain_cutoff_error)     = validate_float(domain_cutoff, "Domain Cutoff", 0.0, 1.0)
    if (domain_cutoff_error == ""):
        domain_cutoff = "%e"%(domain_cutoff)

    if ( correction != "BH" and correction != "BF"):
        (correction, correction_error) = ("BH", "Invalid correction")

# run the perl script, display suitable feedbac
def run_mcam():
    clusterFile = "%sexp_%s_%s.txt"%(clusterPath,exp_id, file_key)
    print "<div class='feedback'>"
    os.system("""perl -I %s %sprintLineToEnrichQueue.pl %s %s %s %f %f %s %f %s %s"""%(libPath,scriptPath,email, clusterFile, exp_id, motif_cutoff,scansite_cutoff, domain_cutoff, alpha_p, correction, PERL_DB))
    print "Your MCAM request has been uploaded.  Notification of completion will be sent to: <b>%s</b>"%(str(email))
    print "</div>"
    print "<fieldset style='margin-bottom:5em;margin-top:1em;'><legend>MCAM Parameters</legend>"
    print "<p>Your parameters were:</p>"
    print "<table style='margin-bottom:1em'>"
    print "<tr><th>Motif Cutoff</th><td>%f</td></tr>"%(float(motif_cutoff))
    print "<tr><th>Scansite Cutoff</th><td>%f</td></tr>"%(float(scansite_cutoff))
    print "<tr><th>Domain Cutoff</th><td>%e</td></tr>"%(float(domain_cutoff))
    print "<tr><th>Alpha P</th><td>%f</td></tr>"%(float(alpha_p))
    print "<tr><th>Correction</th><td>%s</td>"%(str(correction))
    print "</table>"
    
    print "</fieldset>"
    
# handles the user providing MCAM inputs
def process_post():
    validate_inputs();

    if (has_errors()):
        render_form()
    else:
        run_mcam()
    


# renders an error message that couldn't be recovered from
def render_fatal_error(msg):
    print "Error: %s"%(str(msg))
    printFooter();
    sys.exit(0)

def summarize_clusters(filepath):
    clusterSets = ClusterSet(filepath)
    number_of_clusters = len(clusterSets.clusterSets)
    if (clusterSets.parseError != None):
        print "<div class='error'>Error parsing cluster file.  See the help pages for advice on importing clusters.</div>"
    else:
        print "<fieldset style='margin-bottom:1em'><legend>Cluster Information</legend>"
        print "<table>"
        print "<tr><th>Filename</th><td>%s</td></tr>"%(str(user_file_name))
        print "<tr><th>Number of Clusters</th><td>%d</td></tr>"%(int(number_of_clusters))
        print "</table>"
        print "</fieldset>"


def main():
    read_postback_data()

        
    # no exp_id -> error
    if exp_id == "":
        render_fatal_error("No experiment selected")

    if file_key == "":
        render_fatal_error("No clustering file specified")


    # This should really be refaactored and shared between pages
    # If the connection is made global, this becomes really easy
    db = MySQLdb.connect(user= "%s"%user, passwd="%s"%mysql, db="%s"%database)
    c=db.cursor()
    query = """select id from experiment where export = 0"""
    c.execute(query)
    x=c.fetchall()
    exps = [str(item[0]) for item in x]
    if exp_id in exps:
        print "<div class='error'>Illegal experiment id. Subset selection not allowed for this experiment.</div>"
        printFooter()
        sys.exit(0)
        
    printHeader("MCAM", exp_id)
    printExpHeader(exp_id,c)
    
    # if file doesn't exist - render error
    filepath = hackedFileName(exp_id, file_key)
    if (not os.path.exists(filepath)):
        print "<div class='error'>Clusters file not found</div>"
        return

    summarize_clusters(filepath)
    
    if action =="":
        set_default_values()
        render_form()
    else:
        process_post()





            
main()
