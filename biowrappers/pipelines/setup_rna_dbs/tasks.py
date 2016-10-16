import gzip
import os
import zipfile


def file_type(file_name):
    magic_dict = {
        '\x1f\x8b\x08': 'gz',
        '\x50\x4b\x03\x04': 'zip'
    }
    max_len = max(len(x) for x in magic_dict)
    with open(file_name) as f:
        file_start = f.read(max_len)
    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            return filetype
    raise Exception('File type cannot be identified for {}'.format(file_name))


def unzip(in_file, out_file):
    ext = file_type(in_file)
    if ext == 'gz':
        in_file_opener = gzip.GzipFile
    elif ext == 'zip':
        member = os.path.basename(out_file)
        in_file_opener = lambda x: zipfile.ZipFile(x).open(member)
    with in_file_opener(in_file) as in_fh, open(out_file, 'w') as out_fh:
        for line in in_fh:
            out_fh.write(line)
