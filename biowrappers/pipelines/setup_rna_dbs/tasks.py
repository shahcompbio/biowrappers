import gzip
import os
import zipfile


def unzip(in_file, out_file):
    ext = os.path.splitext(in_file)

    if ext == 'gz':
        in_file_opener = gzip.GzipFile
    elif ext == 'zip':
        member = os.path.basename(out_file)
        in_file_opener = lambda x: zipfile.ZipFile(x).open(member)

    with in_file_opener(in_file) as in_fh, open(out_file, 'w') as out_fh:
        for line in in_fh:
            out_fh.write(line)
