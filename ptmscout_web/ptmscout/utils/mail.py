from pyramid_mailer import get_mailer
from ptmscout.config.settings import automailerEmail
from pyramid_mailer.message import Message

def send_automail_message(request, recipients, subject, message):
    mailer = get_mailer(request)

    message = Message(subject=subject,
          sender=automailerEmail,
          recipients=recipients,
          body=message)

    mailer.send_immediately(message)
