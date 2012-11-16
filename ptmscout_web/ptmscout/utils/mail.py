from pyramid_mailer import get_mailer
from ptmscout.config.settings import automailerEmail
from pyramid_mailer.message import Message
import celery.utils
from ptmscout.config import settings

def send_automail_message(request, recipients, subject, message):
    mailer = get_mailer(request)

    message = Message(subject=subject,
          sender=automailerEmail,
          recipients=recipients,
          body=message)

    mailer.send_immediately(message)
    
def celery_send_mail(recipients, subject, message):
    mailer = celery.utils.mail.Mailer()
    message = celery.utils.mail.Message(to=recipients, sender=settings.automailerEmail, subject=subject, body=message)
    mailer.send(message)