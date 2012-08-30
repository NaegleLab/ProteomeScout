from sqlalchemy.ext.declarative import declarative_base

from sqlalchemy.orm import (
    scoped_session,
    sessionmaker,
    )

from zope.sqlalchemy import ZopeTransactionExtension
from zope.sqlalchemy.datamanager import mark_changed
import transaction
from sqlalchemy.exc import IntegrityError

DBSession = scoped_session(sessionmaker(extension=ZopeTransactionExtension()))
Base = declarative_base()


class InternalDatabaseError(Exception):
    pass


def commit():
    try:
        transaction.commit()
    except IntegrityError:
        transaction.abort()
        raise InternalDatabaseError