import os
import subprocess
import sys

install_dir = os.path.dirname(os.path.realpath(__file__))

model_file = os.path.join(install_dir, 'normal_tumour_model.pickle')

def main(args):
    if args[0] == 'classify':
        args.extend(['--model_file', model_file])
    
    cmd = ['museq', ]
    
    cmd.extend(args)
    
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
    
    _, args = parser.parse_known_args()
    
    main(args)
