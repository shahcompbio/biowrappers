import pysftp
import shutil
import urllib2


def download_from_url(src, dst):
    in_fh = urllib2.urlopen(src)
    
    out_fh = open(dst, 'w')
    
    shutil.copyfileobj(in_fh, out_fh)


def download_from_sftp(dst, host, host_path, user, password):
    with pysftp.Connection(host, username=user, password=password) as sftp:
        sftp.get(host_path, localpath=dst)
