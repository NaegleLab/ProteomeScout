from pyramid_mailer import get_mailer
from pyramid_mailer.message import Message
import celery.utils
from ptmscout.config import settings
import re

def send_automail_message(request, recipients, subject, message):
    mailer = get_mailer(request)

    message = Message(subject=subject, sender=settings.automailerEmail,
            recipients=recipients, html=message)

    mailer.send_immediately(message)
    
def celery_send_mail(recipients, subject, message):
    mailer = celery.utils.mail.Mailer()
    message = celery.utils.mail.Message(to=recipients, sender=settings.automailerEmail, subject=subject, body=message)
    mailer.send(message)
    
    
def email_is_valid(email):
    email_regex = re.compile(settings.email_regex, re.I)
    matches = email_regex.match(email)
    
    domain_suffix = ""
    if matches != None:
        domain_suffix = matches.group(1)
    
    return matches != None, domain_suffix.lower() in settings.valid_domain_suffixes
    
