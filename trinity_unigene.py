#!/usr/bin/python
"""
Extract longest transcript from trinity assembly and 
output as unigenes in fasta format in STDOUT or to file
"""
import time
import argparse
import sys


class Ref_Seq(object):
    def __init__(self):
        self._internal_dict = {}
    
    def __setitem__(self, header, sequence):
        self._internal_dict[header]=sequence
    
    def __getitem__(self, header):
        return self._internal_dict[header]
    
    def __iter__(self):
        return self._internal_dict.iteritems()
    
    def __len__(self):
        return len(self._internal_dict)

    def headers(self):
        return self._internal_dict.keys()

    def sequences(self):
        return self._internal_dict.values()    
    

    def load_ref_file(self, ref_file):
        """
        Load ref file in fasta format
    
        Product ref_dict --> header:sequence
        """

        start = time.clock()
        #ref_dict = {}
        ref_count = 0
        loaded_ref = open(ref_file, 'rU')
        full_len_seq = ''
        key = ''
        first_header = True
        for line in loaded_ref:
            if line[0] == '>' and full_len_seq == '':
                key = line.strip().split()[0].split('i')[0][:-1]
                ref_count += 1
                if first_header:
                    first_header = False
            elif line[0] == '>' and full_len_seq != '':
                if key not in self._internal_dict:
                    self._internal_dict.update({key: full_len_seq})
                elif len(full_len_seq) > len(self._internal_dict[key]):
                    self._internal_dict.update({key: full_len_seq})
                key = line.strip().split()[0].split('i')[0][:-1]
                full_len_seq = ''
                ref_count += 1
            elif line[0] == '' and full_len_seq != '':
                if key not in self._internal_dict:
                    self._internal_dict.update({key: full_len_seq})
                elif len(full_len_seq) > len(self._internal_dict[key]):
                    self._internal_dict.update({key: full_len_seq})
                key = line.strip().split()[0].split('i')[0][:-1]
                full_len_seq = ''
            elif line[0] == '':
                pass
            else:
                full_len_seq += line.strip().upper()
    
        if key not in self._internal_dict:
            self._internal_dict.update({key: full_len_seq})
        elif len(full_len_seq) > len(self._internal_dict[key]):
            self._internal_dict.update({key: full_len_seq})
    
        print '\n----{0} reference sequences loaded for alignment----'\
            .format(ref_count)
        if len(self._internal_dict) ==1:
            print "\n{0} length = {1} bp".format(ref_file.split('/')[-1],
                                             len(full_len_seq))
        print "\nReference sequence loading time = "\
+ str((time.clock() - start)) + " seconds\n"

def main():
    parser = argparse.ArgumentParser( description='Filters trinity transcripts to longest unigene')

    ## output file to be written
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to an input FASTA file' )
    parser.add_argument('-o', '--output', type=str, required=False, help='Output file to be created.  Default = STDOUT' )
    args = parser.parse_args()

    fout = sys.stdout
    if args.output is not None:
        fout = open(args.output, 'wt')
    trinity = Ref_Seq()
    trinity.load_ref_file(args.input)

    for key, value in trinity: 
        fout.write(key+"\n")
        lines = (value[n:n + 64] for n in range(0, len(value), 64))
        for i in lines:
            fout.write(i+"\n")


if __name__ == '__main__':
    main()
