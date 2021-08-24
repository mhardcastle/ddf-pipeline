# Basic rclone functionality wrapper which allows you to create an
# object and send multiple commands, capturing output

from __future__ import print_function
import subprocess
import os

def splitlines(s):
    l=s.decode().split('\n')
    if l[-1]=='':
        return l[:-1]
    else:
        return l

class RClone(object):
    def __init__(self, cfg, debug=False):
        try:
            self.command=os.environ['RCLONE_COMMAND']
        except KeyError:
            self.command='rclone'

        try:
            self.config_dir=os.environ['RCLONE_CONFIG_DIR']
        except KeyError:
            self.config_dir=None

        # if no config dir specified, full path should be used

        self.debug=debug
        if self.config_dir is not None:
            self.config_file=os.path.join(self.config_dir,cfg)
        else:
            self.config_file=cfg

        if not os.path.isfile(self.config_file):
            raise RuntimeError('Config file not found at '+self.config_file)
        self.remote=None

    def execute(self,command):
        # generic execution
        
        if isinstance(command,str):
            command=command.split()
        
        fullcommand=[self.command,'--multi-thread-streams','1','--config='+self.config_file]+command
        if self.debug:
            print('Running command',' '.join(fullcommand))
        proc=subprocess.Popen(fullcommand,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        (out, err) = proc.communicate()

        if err:
            print('Rclone command returned error!\n',err)
                
        return {
            "code": proc.returncode,
            "out": splitlines(out),
            "err": splitlines(err)
        }

    def execute_live(self,command):
        # version of the execute command that does *not* catch stdout so you can see what's happening
        if isinstance(command,str):
            command=command.split()
        
        fullcommand=[self.command,'--multi-thread-streams','1','--config='+self.config_file]+command
        if self.debug:
            print('Running command',' '.join(fullcommand))
        proc=subprocess.Popen(fullcommand)
        proc.wait()
        return {"code": proc.returncode, "err": None, "out": None }        
    
    def get_remote(self):
        # If there is only one remote covered by the config file, find out what it is and store in self.remote, else raise exception

        d=self.execute('listremotes')
        if d['code']!=0 or d['err'] or len(d['out'])>1:
            raise RuntimeError('Unable to find unique remote: result was '+repr(d))
        else:
            self.remote=d['out'][0]
            
    def get_dirs(self,remote=None):
        if remote is None:
            if self.remote is None:
                self.get_remote()
            remote=self.remote

        d=self.execute(['lsd',remote])
        return [l.split()[4] for l in d['out']]
    
            
            
