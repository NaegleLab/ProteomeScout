from pyramid_mailer import get_mailer
from pyramid_mailer.message import Message
import celery.utils.mail
from ptmscout.config import settings
import re
from email.mime.text import MIMEText

class CeleryHTMLMessage(celery.utils.mail.Message):
    def __init__(self, to=None, sender=None, subject=None, body=None,
            charset='us-ascii'):
        celery.utils.mail.Message.__init__(self, to, sender, subject, body, charset)

    def __str__(self):
        msg = MIMEText(self.body, 'html', self.charset)
        msg['Subject'] = self.subject
        msg['From'] = self.sender
        msg['To'] = ', '.join(self.to)
        return msg.as_string()




def send_automail_message(request, recipients, subject, message):
    mailer = get_mailer(request)

    message = Message(subject=subject, sender=settings.automailerEmail,
            recipients=recipients, html=message.replace("\n","<br />\n"))

    mailer.send_immediately(message)
   

def celery_send_mail(recipients, subject, message):
    message = message.replace("\n","<br />\n")
    
    mailer = celery.utils.mail.Mailer()
    message = CeleryHTMLMessage(to=recipients, sender=settings.automailerEmail, subject=subject, body=message)
    mailer.send(message)
    
    
def email_is_valid(email):
    email_regex = re.compile(settings.email_regex, re.I)
    matches = email_regex.match(email)
    
    domain = ""
    if matches != None:
        domain = matches.group(1)

    domain_match = False
    for pattern in settings.valid_domain_suffixes:
        success = re.match(pattern, domain) != None
        domain_match |= success

    return matches != None, domain_match
    
