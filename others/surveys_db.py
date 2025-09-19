from __future__ import print_function
from future.utils import itervalues
from builtins import object
import sshtunnel
import socket
import os
import datetime
import numpy as np
from time import sleep
try:
    import MySQLdb as mdb
    import MySQLdb.cursors as mdbcursors
    mdb_type='MySQLdb'
except ImportError:
    import pymysql as mdb
    import pymysql.cursors as mdbcursors
    mdb_type='pymysql'

# General surveys_db functions for use with the various surveys that
# use the central KSP database. Survey-specific functions live in the
# trees for the various surveys.
    
def update_status(name,status,time=None,workdir=None,av=None,survey=None):
    # utility function to just update the status of a field
    # name can be None (work it out from cwd), or string (field name)

    if name is None:
        # work it out
        id=get_id(workdir=workdir)
    else:
        id=name
        
    with SurveysDB(survey=survey) as sdb:
      idd=sdb.get_field(id)
      if idd is None:
          raise RuntimeError('Unable to find database entry for field "%s".' % id)
      idd['status']=status
      tag_field(sdb,idd,workdir=workdir)
      if time is not None and idd[time] is None:
        idd[time]=datetime.datetime.now()
      if av is not None:
        idd['archive_version']=av
      sdb.set_field(idd)

def tag_field(sdb,idd,workdir=None):
    # Add location and user tags
    idd['clustername']=get_cluster()
    if workdir is None:
        idd['location']=os.getcwd()
    else:
        idd['location']=workdir
    idd['username']=get_user()
    idd['nodename']=sdb.hostname
    
def get_id(workdir=None):
    if workdir is None:
        dir=os.getcwd()
    else:
        dir=workdir
    dname=os.path.basename(dir)
    return dname

def get_user():
    return os.getenv('USER')

def get_cluster():
    cluster=os.getenv('DDF_PIPELINE_CLUSTER')
    if cluster:
        return cluster
    else:
        return 'Unknown'

def use_database():
    return 'DDF_PIPELINE_DATABASE' in os.environ

class SurveysDB(object):
    ''' Provides low-level and high-level interfaces to the surveys database '''

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        self.close()

    def __init__(self,readonly=False,verbose=False,survey=None,retrynumber=10,sleeptime=60):
        if verbose: print('MySQL type is',mdb_type)
        if survey is None:
            survey='hba' # preserve old default behaviour
        # get the config file -- this must exist
        home=os.getenv("HOME")
        mysql_host=os.getenv('DDF_PIPELINE_MYSQLHOST')
        if not mysql_host:
            mysql_host='lofar-server.data'
        if verbose:
            print('MySQL host is',mysql_host)
        cfg=[l.rstrip() for l in open(home+'/.surveys').readlines()]
        self.password=cfg[0]
        try:
            self.ssh_user=cfg[1]
        except:
            self.ssh_user=None
        
        try:
            self.ssh_key=cfg[2]
        except:
            self.ssh_key="id_rsa"

        self.readonly=readonly
        self.verbose=verbose
        self.survey=survey
        if self.survey=='hba':
            self.database='surveys'
        elif self.survey=='lba':
            self.database='lba'
        else:
            raise NotImplementedError('Survey "%s" not known' % self.survey)
        
        # set up an ssh tunnel if not running locally
        self.usetunnel=False
        self.hostname=socket.gethostname()
        if self.hostname=='lofar-server':
            # Running on the database server itself
            if verbose:
                print('Using direct connection to localhost')
            if mdb_type=='MySQLdb':
                self.con=mdb.connect('127.0.0.1', 'survey_user', self.password, self.database, cursorclass=mdbcursors.DictCursor)
            else:
                self.con=mdb.connect(host='127.0.0.1', user='survey_user', password=self.password, database=self.database, cursorclass=mdbcursors.DictCursor)
        else:
            try:
                dummy=socket.gethostbyname(mysql_host)
            except socket.gaierror:
                if verbose:
                    print('Cannot find host',mysql_host,'will use tunnel')
                self.usetunnel=True

            if self.usetunnel:
                # Create an ssh tunnel and connect through that
                if verbose:
                    logger=sshtunnel.create_logger(loglevel=10)
                else:
                    logger=None
                # Reading the key ensures an error if it doesn't exist
                self.pkey=sshtunnel.SSHTunnelForwarder.read_private_key_file(home+'/.ssh/'+self.ssh_key)
                connected=False
                retry=0
                while not connected and retry<retrynumber:
                    self.tunnel=sshtunnel.SSHTunnelForwarder('lofar.herts.ac.uk',
                                                             ssh_username=self.ssh_user,
                                                             ssh_pkey=self.pkey,
                                                             remote_bind_address=('127.0.0.1',3306),
                                                             local_bind_address=('127.0.0.1',),
                                                             host_pkey_directories=[],
                                                             allow_agent=True,
                                                             logger=logger
                                                             )

                    try:
                        self.tunnel.start()
                    except sshtunnel.BaseSSHTunnelForwarderError as e:
                        print('ssh tunnel temporary error %s! Sleep %i seconds to retry' % (e,sleeptime))
                        retry+=1
                        sleep(sleeptime)
                        continue
                    localport=self.tunnel.local_bind_port
                    try:
                        if mdb_type=='MySQLdb':
                            self.con = mdb.connect('127.0.0.1', 'survey_user', self.password, self.database, port=localport, cursorclass=mdbcursors.DictCursor)
                        else:
                            self.con = mdb.connect(host='127.0.0.1', user='survey_user', password=self.password, database=self.database, port=localport, cursorclass=mdbcursors.DictCursor)
                        connected=True
                    except mdb.OperationalError as e:
                        print('Database temporary error %s! Sleep %i seconds to retry\n' % (e,sleeptime))
                        retry+=1
                        sleep(sleeptime)
            else:
                # Network connection not using ssh tunnel
                connected=False
                retry=0
                while not connected and retry<retrynumber:
                    try:
                        if mdb_type=='MySQLdb':
                            self.con = mdb.connect(mysql_host, 'survey_user', self.password, self.database, cursorclass=mdbcursors.DictCursor)
                        else:
                            self.con = mdb.connect(host=mysql_host, user='survey_user', password=self.password, database=self.database, cursorclass=mdbcursors.DictCursor)                
                        connected=True
                    except mdb.OperationalError as e:
                        print('Database temporary error! Sleep %i seconds to retry\n' % sleeptime,e)
                        retry+=1
                        sleep(sleeptime)
                if not connected:
                    raise RuntimeError("Cannot connect to database server after repeated retry")
        self.cur = self.con.cursor()

        # get the tables list for locking
        self.cur.execute('show tables')
        result=self.cur.fetchall()
        self.tables=[list(itervalues(d))[0] for d in result]
        
        if self.readonly:
            pass
            #can't use this feature on lofar's version of MariaDB
            #self.cur.execute('set session transaction read only')
        else:
            command='lock table '
            for table in self.tables:
                command+=table+' write, '
            command=command[:-2]
            self.cur.execute(command)
        self.closed=False

    def execute(self,*args):
        if self.verbose:
            print(args)
        self.cur.execute(*args)

    def close(self):
        # if 'closed' doesn't exist, then we are most likely being called through __del__ due to a failure in the init call. So skip the rest.
        if hasattr(self,'closed'):
            if not self.closed:
                if not self.readonly:
                    self.cur.execute('unlock tables')
                self.con.close()
                if self.usetunnel:
                    self.tunnel.stop()
                self.closed=True # prevent del from trying again
    
    def __del__(self):
        self.close()

    def check_table(self,table):
        if table not in self.tables:
            table+='s'
            if table not in self.tables:
                raise RuntimeError('Unknown table %s requested' % table)
        return table
        
    def db_get(self,table,id):
        if self.closed: raise RuntimeError('Attempting DB operation but instance is closed')
        table=self.check_table(table)
        self.execute('select * from '+table+' where id=%s',(id,))
        result=self.cur.fetchall()
        if len(result)==0:
            return None
        else:
            return result[0]

    def db_set(self,table,record):
        if self.closed: raise RuntimeError('Attempting DB operation but instance is closed')
        if self.readonly: raise RuntimeError('Write requested in read-only mode')
        table=self.check_table(table)
        id=record['id'];
        for k in record:
            if k=='id':
                continue
            if record[k] is not None:
                if isinstance(record[k],float) and np.isnan(record[k]):
                    record[k]=None # should work for NULL
                self.execute('update '+table+' set '+k+'=%s where id=%s',(record[k],id))

    def db_create(self,table,id):
        if self.closed: raise RuntimeError('Attempting DB operation but instance is closed')
        table=self.check_table(table)
        if self.readonly: raise RuntimeError('Create requested in read-only mode')
        self.execute('insert into '+table+'(id) values (%s)',(id,))
        return self.db_get(table,id)

    def db_delete(self,table,id):
        if self.closed: raise RuntimeError('Attempting DB operation but instance is closed')
        table=self.check_table(table)
        if self.readonly: raise RuntimeError('Delete requested in read-only mode')
        self.execute('delete from '+table+' where id=%s',(id,))
        return None
    
    def get_field(self,id):
        return self.db_get('fields',id)

    def set_field(self,sd):
        self.db_set('fields',sd)

    def create_field(self,id):
        return self.db_create('fields',id)

    def get_observation(self,id):
        return self.db_get('observations',id)
    
    def set_observation(self,sd):
        self.db_set('observations',sd)
        
    def create_observation(self,id):
        self.db_create('observations',id)

    def get_transient(self,id):
        return self.db_get('transients',id)

    def set_transient(self,sd):
        return self.db_set('transients',sd)
        
    def create_transient(self,id):
        return self.db_create('transients',id)

    def create_quality(self,id):
        self.cur.execute('delete from quality where id="%s"' % id)
        return self.db_create('quality',id)

    def get_quality(self,id):
        return self.db_get('quality',id)
        
    def set_quality(self,sd):
        return self.db_set('quality',sd)
    
    def get_reprocessing(self,id):
        return self.db_get('reprocessing',id)
    
    def set_reprocessing(self,sd):
        self.db_set('reprocessing',sd)
        
    def create_reprocessing(self,id):
        self.db_create('reprocessing',id)

    def get_ffr(self,id,operation):
        table='full_field_reprocessing'
        self.execute('select * from '+table+' where id=%s and operation=%s',(id,operation))
        result=self.cur.fetchall()
        if len(result)==0:
            return None
        else:
            return result[0]

    def set_ffr(self,record):
        table='full_field_reprocessing'
        if self.readonly: raise RuntimeError('Write requested in read-only mode')
        id=record['id']
        operation=record['operation']
        for k in record:
            if k=='id' or k=='operation':
                continue
            if record[k] is not None:
                self.execute('update '+table+' set '+k+'=%s where id=%s and operation=%s',(record[k],id,operation))


    def create_ffr(self,id,operation):
        table='full_field_reprocessing'
        if self.readonly: raise RuntimeError('Create requested in read-only mode')
        self.execute('insert into '+table+'(id,operation) values (%s,%s)',(id,operation))
        return self.get_ffr(id,operation)

    def delete_ffr(self,id,operation):
        table='full_field_reprocessing'
        if self.readonly: raise RuntimeError('Delete requested in read-only mode')
        self.execute('delete from '+table+' where id=%s and operation=%s',(id,operation))
        

if __name__=='__main__':
    with SurveysDB(verbose=True,survey='hba') as sdb:
        result=sdb.db_get('fields','P35Hetdex10')
    print(result)

