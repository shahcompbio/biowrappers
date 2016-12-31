import pypeliner


def gunzip(gzipped_file):
    pypeliner.commandline.execute('gunzip', gzipped_file)
