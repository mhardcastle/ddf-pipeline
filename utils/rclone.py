# Basic rclone functionality thin wrapper which allows you to create an
# object and send multiple commands, capturing output where necessary.

from __future__ import print_function
import subprocess
import os
import tempfile

def splitlines(s):
    l=s.decode().split('\n')
    if l[-1]=='':
        return l[:-1]
    else:
        return l

class RClone(object):
    def __init__(self, cfg, debug=False):
        # Set up environment variables or sensible defaults
        try:
            self.command=os.environ['RCLONE_COMMAND']
        except KeyError:
            self.command='rclone'

        try:
            self.ada_command=os.environ['ADA_COMMAND']
        except KeyError:
            self.ada_command='ada'

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
        '''
        generic execution with standard out and error caught so that
        they can be parsed. Command is a string or a list that can be passed to Popen. stdout and stderr are caught and returned as elements of a dictionary along with any return code.
        '''
        
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
        '''
        version of the execute command that does *not* catch stdout so you can see what's happening. Returns a dictionary of the same format as execute for consistency but as stdout and stderr are caught they are always None
        '''
        
        if isinstance(command,str):
            command=command.split()
        
        fullcommand=[self.command,'--multi-thread-streams','1','--config='+self.config_file]+command
        if self.debug:
            print('Running command',' '.join(fullcommand))
        proc=subprocess.Popen(fullcommand)
        proc.wait()
        return {"code": proc.returncode, "err": None, "out": None }        

    def copy(self,source,dest):
        '''
        simplifying wrapper function -- one of source and dest needs
        to contain a 'remote' specification, probably self.remote,
        for this to do anything useful. As with rclone copy this will
        work on a single file or a whole directory.
        '''
        
        return self.execute_live(['-P','copy',source,dest])

    def multicopy(self,sourcedir,files,dest):
        '''
        another wrapper function, this time copy named files from
        the source directory to the destination. Better to use this
        than looping over copy if e.g. you want to exploit
        multi-threading or stage more than one file at a time but do
        not want to copy a whole directory.
        '''
        
        with tempfile.NamedTemporaryFile(suffix='.txt',delete=False,mode='w') as outfile:
            filename=outfile.name
            outfile.writelines([f+'\n' for f in files])
        result=self.execute_live(['-P','--include-from',filename,'copy',sourcedir,dest])
        os.unlink(filename)
        return result
        
    def get_remote(self):
        '''
        If there is only one remote covered by the config file, find out what it is and store in self.remote, else raise exception
        '''

        d=self.execute('listremotes')
        if d['code']!=0 or d['err'] or len(d['out'])>1:
            raise RuntimeError('Unable to find unique remote: result was '+repr(d))
        else:
            self.remote=d['out'][0]
            
    def get_dirs(self,base='',remote=None):
        '''
        wrapper round rclone lsd that returns a list of directories either in the root of the remote or in a specified base directory. If no remote specified use the result of get_remote().
        '''
        if remote is None:
            if self.remote is None:
                self.get_remote()
            remote=self.remote

        d=self.execute(['lsd',remote+base])
        return [l.split()[4] for l in d['out']]
    
    def get_files(self,base='',remote=None, exclude_dirs=True):
        '''
        wrapper round rclone lsf that returns a list of files either in the root of the remote or in a specified base directory. If no remote specified use the result of get_remote().
        '''
        if remote is None:
            if self.remote is None:
                self.get_remote()
            remote=self.remote

        d=self.execute(['lsf',remote+base])
        return [l for l in d['out'] if not exclude_dirs or not l.endswith('/')]
    
    def get_checksum(self,filename):
        ''' Use ada to get the checksum. Filename is the remote filename. ada does not use the remote. ada config file should contain the API information. '''
        command=self.ada_command+' --tokenfile '+self.config_file+' --checksum %s'% filename
        if self.debug:
            print('Running '+command)
        t = os.popen(command).read()
        return t.split()[1].replace('ADLER32=','')
