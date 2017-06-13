import subprocess
import os

FILE =  os.path.abspath(__file__)
DIR = "/".join(FILE.split("/")[:-1])

print "add this following line to your .bashrc (linux) or .base_profile (mac)"
print "export  PYTHONPATH=$PYTHONPATH:"+DIR

print "installing all required packages"
subprocess.call('sudo pip install -r requirements.txt', shell=True)
