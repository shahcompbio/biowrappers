'''
Created on 3 Apr 2017

@author: Andrew Roth
'''
from collections import OrderedDict

import re
import vcf


class SnpEffParser(object):

    def __init__(self, file_name):
        self._reader = vcf.Reader(filename=file_name)

        self.fields = self._get_field_names()

        self._buffer = []

    def __iter__(self):
        while True:
            yield self.next()

    def next(self):
        while len(self._buffer) == 0:
            record = self._reader.next()

            if 'ANN' not in record.INFO:
                continue

            for row in self._parse_record(record):
                self._buffer.append(row)

        return self._buffer.pop(0)

    def _get_field_names(self):
        fields = []

        match = re.search(":(.*)", self._reader.infos['ANN'].desc).groups()[0].replace("'", "")

        for x in match.split('|'):
            fields.append(x.strip().lower())

        return fields

    def _parse_record(self, record):
        for annotation in record.INFO['ANN']:
            out_row = OrderedDict((
                ('chrom', record.CHROM),
                ('coord', record.POS),
                ('ref', record.REF),
                ('alt', ','.join([str(x) for x in record.ALT])),
            ))

            fields = annotation.split('|')

            for i, key in enumerate(self.fields):
                out_row[key] = fields[i]

            yield out_row


class ClassicSnpEffParser(object):

    def __init__(self, file_name):
        self._reader = vcf.Reader(filename=file_name)

        self.fields = self._get_field_names()

        self._buffer = []

        self._effect_matcher = re.compile(r'(.*)\(')

        self._fields_matcher = re.compile(r'\((.*)\)')

    def __iter__(self):
        while True:
            yield self.next()

    def next(self):
        while len(self._buffer) == 0:
            record = self._reader.next()

            if 'EFF' not in record.INFO:
                continue

            for row in self._parse_record(record):
                self._buffer.append(row)

        return self._buffer.pop(0)

    def _get_field_names(self):
        fields = []

        match = re.search(r'\((.*)\[', self._reader.infos['EFF'].desc)

        for x in match.groups()[0].split('|'):
            fields.append(x.strip().lower())

        return fields

    def _parse_record(self, record):
        for annotation in record.INFO['EFF']:
            effect = self._effect_matcher.search(annotation).groups()[0]

            out_row = OrderedDict((
                ('chrom', record.CHROM),
                ('coord', record.POS),
                ('ref', record.REF),
                ('alt', ','.join([str(x) for x in record.ALT])),
                ('effect', effect),
            ))

            fields = self._fields_matcher.search(annotation).groups()[0].split('|')

            for i, key in enumerate(self.fields):
                out_row[key] = fields[i]

            yield out_row
