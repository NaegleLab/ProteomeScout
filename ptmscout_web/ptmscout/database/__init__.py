from sqlalchemy.ext.declarative import declarative_base

from sqlalchemy.orm import (
    scoped_session,
    sessionmaker,
    )

from zope.sqlalchemy import ZopeTransactionExtension
from zope.sqlalchemy.datamanager import mark_changed

DBSession = scoped_session(sessionmaker(extension=ZopeTransactionExtension()))
Base = declarative_base()

from uniprot import *
from taxonomies import *

from jobs import *

from experiment import *
from permissions import *
from user import *

from modifications import *
from mutations import *
from gene_expression import *
from protein import *

from annotations import *