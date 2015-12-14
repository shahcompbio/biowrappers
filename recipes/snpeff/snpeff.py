#!/usr/bin/env python

import os
import subprocess
import sys

install_dir = os.path.dirname(os.path.realpath(__file__))
    
def main(args, snpeff_args):
    
    jar_file = os.path.join(install_dir, 'snpEff.jar')
    
    os.environ['MALLOC_ARENA_MAX'] = '4'
    
    cmd = [
        'java',
        '-Xms{0}g'.format(args.memory),
        '-Xmx{0}g'.format(args.memory),
        '-XX:CompressedClassSpaceSize=64m',
        '-XX:MaxMetaspaceSize=128m',
        '-XX:+UseConcMarkSweepGC',
        '-jar',
        jar_file
    ]
    
    cmd.extend(snpeff_args)

    run_cmd(cmd)

def run_cmd(cmd):

    cmd = [str(x) for x in cmd]
    
    process = subprocess.Popen(
        cmd, 
        stdin=None, 
        stdout=None, 
        stderr=None, 
        shell=False
    )
    
    stdout, stderr = process.communicate()
    
    sys.exit(process.returncode)
  
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--memory', default=4, type=int)
        
    args, unknown = parser.parse_known_args()
    
    main(args, unknown)