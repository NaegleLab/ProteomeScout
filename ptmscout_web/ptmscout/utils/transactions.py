import transaction
from sqlalchemy.exc import IntegrityError
from smtplib import SMTPRecipientsRefused

class InternalDatabaseError(Exception):
    pass

class MailerError(Exception):
    pass

def commit():
    try:
        transaction.commit()
    except IntegrityError:
        transaction.abort()
        raise InternalDatabaseError
    except SMTPRecipientsRefused:
        transaction.abort()
        raise MailerError