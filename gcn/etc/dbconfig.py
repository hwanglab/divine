#
# COPYRIGHT (C) 2002-2011 Rajgopal Srinivasan
#
"""
.. module:: dbconfig
    :platform: Unix, Windows, MacOSX
    :synopsis: Configuration for Database connections

.. moduleauthor:: Rajgopal Srinivasan (rajgopal.srinivasan@tcs.com); modified by changjin.hong@gmail.com

The module defines a dictionary *DBCONFIG* that provides the parameters needed
to connect to various databases that are used in the Genome Commons Navigator.

The dictionary has as key the name of the resource for which a db connection
exist (e.g resources are REFGENE, REFMRNA). The value for each key is also
a dictionary with the following keys

    - db: The database type (sqlite3, postgres, mysql, ...)

    - name: Name of the database. For sqlite3 this the fully qualified
            filename of the database

    - host: The name (or IP Number) of the host on which the database is
            deployed
    - port: The port through which the database connection is available

    - user: The login name for the database

    - password: The password corresponding to the `user` login.

The parameters `host`, `port`, `user` and `password` may all be ignored for
sqlite3 and similar databases
"""
import os
from gcn.config import lib_config

_LOCALDIR = '/opt/db'


def set_localdir():
    """Set the local database directory from the enviornmental variable
    GCN_DB_DIR"""
    local = lib_config.gcn_path('GCN_DB_DIR')
    if local is not None:
        global _LOCALDIR
        _LOCALDIR = local

set_localdir()

DBCONFIG = {
    'REFMRNA': {'db': 'sqlite3',
                'name': os.path.join(_LOCALDIR, 'refmrna'),
                'tables':[],
                'port': None,
                'host': None,
                'user': None,
                'password': None,
                },
    'REFGENE': {'db': 'sqlite3',
                'name': os.path.join(_LOCALDIR, 'refgene'),
                'tables':[],
                'port': None,
                'host': None,
                'user': None,
                'password': None,
                },
    'GENEONTOLOGY': {'db': 'sqlite3',
                'name': os.path.join(_LOCALDIR, 'geneontology'),
                'tables':[],
                'port': None,
                'host': None,
                'user': None,
                'password': None,
                },
    'MIRNA': {'db': 'sqlite3',
              'name': os.path.join(_LOCALDIR, 'mirna'),
              'tables':[],
              'port': None,
              'host': None,
              'user': None,
              'password': None,
              },
    'UTRDB': {'db': 'sqlite3',
              'name': os.path.join(_LOCALDIR, 'utrdb'),
              'tables':[],
              'port': None,
              'host': None,
              'user': None,
              'password': None,
              },
    'REGULOME': {'db': 'sqlite3',
              'name': os.path.join(_LOCALDIR, 'regulomedb'),
              'tables':[],
              'port': None,
              'host': None,
              'user': None,
              'password': None,
              },
    'KGDB': {'db': 'sqlite3',
             'name': os.path.join(_LOCALDIR, 'kgdb'),
             'tables':[],
             'port': None,
             'host': None,
             'user': None,
             'password': None,
             },
    'DBSNP': {'db': 'sqlite3',
             'name': os.path.join(_LOCALDIR, 'dbsnp'),
             'port': None,
             'host': None,
             'user': None,
             'password': None,
         },
    'ESP': {'db': 'sqlite3',
             'name': os.path.join(_LOCALDIR, 'esp'),
             'port': None,
             'host': None,
             'user': None,
             'password': None,
         },
    'EXAC': {'db': 'sqlite3',
             'name': os.path.join(_LOCALDIR, 'exac'),
             'port': None,
             'host': None,
             'user': None,
             'password': None,
         },
    'HGMDDB': {'db': 'sqlite3',
             'name': os.path.join(_LOCALDIR, 'hgmd'),
             'port': None,
             'host': None,
             'user': None,
             'password': None,
         },
    'CLINVITAE': {'db': 'sqlite3',
              'name': os.path.join(_LOCALDIR,'clinvitae'),
              'port': None,
              'host': None,
              'user': None,
              'password': None,
         },
    'COSMIC': {'db': 'sqlite3',
              'name': os.path.join(_LOCALDIR,'cosmic'),
              'port': None,
              'host': None,
              'user': None,
              'password': None,
         },
    'CLINVARDB': {'db': 'sqlite3',
                  'name': os.path.join(_LOCALDIR, 'clinvardb'),
                  'port': None,
                  'host': None,
                  'user': None,
                  'password': None,
                  },
    'CLNPHESNP': {'db': 'sqlite3',
                'name': os.path.join(_LOCALDIR, 'clnphesnpdb'),
                'port': None,
                'host': None,
                'user': None,
                'password': None,
                 },
    'OMIM': {'db': 'sqlite3',
            'name': os.path.join(_LOCALDIR, 'mimdb'),
            'port': None,
            'host': None,
            'user': None,
            'password': None,
            },
    'NSFP': {'db': 'sqlite3',
            'name': os.path.join(_LOCALDIR, 'nsfpdb'),
            'port': None,
            'host': None,
            'user': None,
            'password': None,
            },
    'SPLICE': {'db': 'sqlite3',
               'name': os.path.join(_LOCALDIR, 'splicedb'),
               'port': None,
               'host': None,
               'user': None,
               'password': None,
               },
    'GENCODE': {'db': 'sqlite3',
                'name': os.path.join(_LOCALDIR, 'gencode'),
                'port': None,
                'host': None,
                'user': None,
                'password': None,
                },
    'GENCODEMRNA': {'db': 'sqlite3',
                    'name': os.path.join(_LOCALDIR, 'gencodemrna'),
                    'port': None,
                    'host': None,
                    'user': None,
                    'password': None,
                    },
    'INTERPRO': {'db': 'sqlite3',
             'name': os.path.join(_LOCALDIR, 'interpro'),
             'port': None,
             'host': None,
             'user': None,
            'password': None,
            },
}


def configure_db(resource, db, name, port=None, host=None, user=None,
                 password=None):
    """Configure a database connection.

    Args:
        resource (str): Name of the resource (e.g REFMRNA) for which the
                        database information is being configured
        name (str):     Name/Filename of the database
        host (str):     Name of the host on which the database is deployed
        port (integer): Port through which the database is to be accessed
        user (str):     Login id to connect to the database
        password (str): Password for the `user` to connect to the database

    Returns:
        None. Configures the database connection for the resource.

    This method is useful when individual users want to use custom resources
    instead of the system defaults.

    """

    DBCONFIG[resource] = {'db': db, 'name': name, 'port': port,
                          'host': host, 'user': user, 'password': password}

