import urllib2

uniid="P07327"

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

try:
#   response = urllib2.urlopen('https://example.com') 
    response = urllib2.urlopen('https://www.uniprot.org/uniprot/'+uniid+'.fasta') 
    print 'response headers: "%s"' % response.info()
except IOError, e:
    if hasattr(e, 'code'): # HTTPError
        print 'http error code: ', e.code
    elif hasattr(e, 'reason'): # URLError
        print "can't connect, reason: ", e.reason
    else:
        raise
