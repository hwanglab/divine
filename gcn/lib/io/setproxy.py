#
# COPYRIGHT (C) 2002-2010 Rajgopal Srinivasan
#
"""
.. module:: entrystream
    :platform: Unix, Windows, MacOSX
    :synopsis: Functions to setup access through an authenticated web proxy

.. moduleauthor:: Rajgopal Srinivasan (rajgopal.srinivasan@gmail.com); modified by changjin.hong@gmail.com


Helper module to make it easy to work with proxies.  Assuming you are
behind a proxy you should first call setup_proxy to set up the proxy
parameters (including username and password if you are behind a proxy
that requires authentication).  Once that has been setup you can just
use the functions in the urllib2 module to retrieve data.  This module
enables proxy support for both http and ftp.
"""
import urllib2

def hexpass(pwd):
    """Convert no alpha numeric characters in a password string to
    hexadecimal

    Args:
        pwd(str): The password as a string

    Returns:
        The converted string

    Examples

    >>> print(hexpass('matter'))
    'matter'

    >>> print(hexpass('m@tt3r'))
    "m%40tt3r'

    """

    return ''.join([c if c.isalnum() else _ashex(c) for c in pwd])

def _ashex(char):
    return '%' + hex(ord(char))[2:]

def setup_proxy(host, port, username='', password=''):
    """Enable proxy support for accessing an external network

    Args:
        host (str): Name of the proxy host to connect to
        port (int): Port number on which proxy is listening
        username (str): Name to login as (in case of authenticated proxy)
        password (str): plain text password for proxy authentication

    Returns:
        None

    Examples:

    For use with a proxy not requiring authentatication

    >>> setup_proxy('10.0.0.1', 3128)

    and for one requiring authentication

    >>> setup_proxy('192.168.1.2', 3128, 'me', 'mypassword')

    """

    pwd = hexpass(password)

    prefix = ''
    if username:
        prefix += username
    if pwd:
        prefix = prefix + ':' + pwd

    if prefix:
        prefix += '@'

    auth = 'http://%s%s:%s' % (prefix, host, port)

    proxy_handler = urllib2.ProxyHandler({'ftp': auth,
                                          'http': auth})

    opener = urllib2.build_opener(proxy_handler,
                                  urllib2.HTTPBasicAuthHandler(),
                                  urllib2.HTTPHandler,
                                  urllib2.HTTPSHandler,
                                  urllib2.FTPHandler)
    urllib2.install_opener(opener)

