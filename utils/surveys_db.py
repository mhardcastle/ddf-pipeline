from __future__ import print_function
from builtins import object
import sshtunnel
import socket
import os
import datetime
from time import sleep
try:
    import MySQLdb as mdb
    import MySQLdb.cursors as mdbcursors
except ImportError:
    import pymysql as mdb
    import pymysql.cursors as mdbcursors


def get_next():
    # return the name of the top-priority field with appropriate status
    sdb=SurveysDB(readonly=True)
    sdb.cur.execute('select fields.id as id,sum(nsb*integration/232) as s,count(observations.id) as c,fields.required_integration as ri, fields.priority,fields.lotss_field from fields left join observations on (observations.field=fields.id) where fields.status="Not started" and (observations.status="Archived" or observations.status="DI_processed") and (gal_b>10 or gal_b<-10 or gal_b is NULL or fields.priority>9) group by fields.id having (s>0.95*ri or lotss_field=0) order by fields.priority desc,ra desc')
    results=sdb.cur.fetchall()
    sdb.close()
    if len(results)>0:
        return results[0]['id']
    else:
        return None

def get_next_selfcalibration():
    sdb=SurveysDB(readonly=True)
    sdb.cur.execute('select reprocessing.id,reprocessing.priority,reprocessing.fields,reprocessing.extract_status from reprocessing where reprocessing.selfcal_status like "%SREADY%" group by reprocessing.priority desc')
    results=sdb.cur.fetchall()
    sdb.close()
    if len(results)>0:
        return results[0]
    else:
        return None

def get_next_extraction():
    # return the name of the top-priority field with appropriate status
    sdb=SurveysDB(readonly=True)
    sdb.cur.execute('select reprocessing.id,reprocessing.priority,reprocessing.fields,reprocessing.extract_status from reprocessing where reprocessing.extract_status like "%EREADY%" group by reprocessing.priority desc')
    results=sdb.cur.fetchall()
    sdb.close()
    if len(results)>0:
        return results[0]
    else:
        return None

def update_status(name,status,time=None,workdir=None,av=None):
    # utility function to just update the status of a field
    # name can be None (work it out from cwd), or string (field name)

    if name is None:
        # work it out
        id=get_id(workdir=workdir)
    else:
        id=name
        
    with SurveysDB() as sdb:
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

    def __init__(self,readonly=False,verbose=False):

        # get the config file -- this must exist
        home=os.getenv("HOME")
        mysql_host=os.getenv('DDF_PIPELINE_MYSQLHOST')
        if not mysql_host:
            mysql_host='lofar-server.data'
        if verbose:
            print('MySQL host is',mysql_host)
        cfg=open(home+'/.surveys').readlines()
        self.password=cfg[0].rstrip()
        try:
            self.ssh_user=cfg[1].rstrip()
        except:
            self.ssh_user=None
        
        try:
            self.ssh_key=cfg[2].rstrip()
        except:
            self.ssh_key="id_rsa"

        # read only use
        self.readonly=readonly
        self.verbose=verbose

        self.tables=['fields','observations','quality','transients','reprocessing']
        
        # set up an ssh tunnel if not running locally
        self.usetunnel=False
        self.hostname=socket.gethostname()
        if self.hostname=='lofar-server':
            if verbose:
                print('Using direct connection to localhost')
            self.con = mdb.connect('127.0.0.1', 'survey_user', self.password, 'surveys',cursorclass=mdbcursors.DictCursor)
        else:
            try:
                dummy=socket.gethostbyname(mysql_host)
            except socket.gaierror:
                if verbose:
                    print('Cannot find host',mysql_host,'will use tunnel')
                self.usetunnel=True

            if self.usetunnel:
                self.tunnel=sshtunnel.SSHTunnelForwarder('lofar.herts.ac.uk',
                                                         ssh_username=self.ssh_user,
                                                         ssh_pkey=home+'/.ssh/%s'%self.ssh_key,
                                                         remote_bind_address=('127.0.0.1',3306),
                                                         local_bind_address=('127.0.0.1',))

                self.tunnel.start()
                localport=self.tunnel.local_bind_port
                self.con = mdb.connect('127.0.0.1', 'survey_user', self.password, 'surveys', port=localport, cursorclass=mdbcursors.DictCursor)
            else:
                connected=False
                retry=0
                while not connected and retry<10:
                    try:
                        self.con = mdb.connect(mysql_host, 'survey_user', self.password, 'surveys',cursorclass=mdbcursors.DictCursor)
                        connected=True
                    except mdb.OperationalError as e:
                        print('Database temporary error! Sleep to retry',e)
                        retry+=1
                        sleep(20)
                if not connected:
                    raise RuntimeError("Cannot connect to database server")
        self.cur = self.con.cursor()
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
        table=self.check_table(table)
        self.execute('select * from '+table+' where id=%s',(id,))
        result=self.cur.fetchall()
        if len(result)==0:
            return None
        else:
            return result[0]

    def db_set(self,table,record):
        if self.readonly: raise RuntimeError('Write requested in read-only mode')
        table=self.check_table(table)
        id=record['id'];
        for k in record:
            if k=='id':
                continue
            if record[k] is not None:
                self.execute('update '+table+' set '+k+'=%s where id=%s',(record[k],id))

    def db_create(self,table,id):
        table=self.check_table(table)
        if self.readonly: raise RuntimeError('Create requested in read-only mode')
        self.execute('insert into '+table+'(id) values (%s)',(id,))
        return self.db_get(table,id)

    def db_delete(self,table,id):
        table=self.check_table(table)
        if self.readonly: raise RuntimeError('Create requested in read-only mode')
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


if __name__=='__main__':
    sdb=SurveysDB(verbose=True)
    result=sdb.db_get('fields','P35Hetdex10')
    #result['location']='Never Never Land'
    #sdb.set_id(result)
    print(result)
    sdb.close()

