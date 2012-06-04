### usage
### python parseErrors.py userEmail experimentid notes
###
from pathname import *
import sys
from  smtplib import SMTP,SMTPSenderRefused

sys.path.append(path)
from MimeWriter import MimeWriter
try:
	from cStringIO import StringIO
except:
	from StringIO import StringIO
def main():
    email = sys.argv[1]
    expid = sys.argv[2]
    notes = sys.argv[3]
    text ="""
    Dear user,

    The log from your data load can be found here: %sload/errors.cgi?expid=%s.

    Notes:
    %s
    """%(urlpath,expid,notes)
    tempfile = StringIO()
    mw = MimeWriter(tempfile)
    mw.addheader('to',email)
    mw.addheader('from',sysEmail)
    mw.addheader('subject','[PTMScout] Dataset Load')
    mw.startmultipartbody('mixed')
    sw=mw.nextpart()
    j = sw.startbody('text/plain')
    j.write(text)
    mw.lastpart()
    message = tempfile.getvalue()
    mailer = SMTP()
    mailer.connect("localhost")
    ### make text of email
    
    print text
    ### write email with link to html error page
    try:
		mailer.sendmail(sysEmail,email,message)
    except:
        print "Invalid sender email address"




if len(sys.argv)<3:
    print """Usage: python parseErrors.py userEmail experimentid notes"""
else:
    main()
