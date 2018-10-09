import sshtunnel
import socket
import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors
import os
import datetime

def get_next():
    # return the name of the top-priority field with appropriate status
    sdb=SurveysDB(readonly=True)
    sdb.cur.execute('select fields.id as id,sum(nsb*integration/232) as s,count(observations.id) as c,fields.priority from fields left join observations on (observations.field=fields.id) where fields.status="Not started" and observations.status="DI_processed" group by fields.id having s>7 order by fields.priority desc,ra desc')
    results=sdb.cur.fetchall()
    sdb.close()
    return results[0]['id']

def update_status(name,status,time=None,workdir=None):
    # utility function to just update the status of an observation
    # name can be None (work it out from cwd), or string (field name)

    if name is None:
        # work it out
        id=get_id(workdir=workdir)
    else:
        id=name
        
    sdb=SurveysDB()
    idd=sdb.get_field(id)
    idd['status']=status
    tag_field(sdb,idd,workdir=workdir)
    if time is not None and idd[time] is None:
        idd[time]=datetime.datetime.now()
    sdb.set_field(idd)
    sdb.close()

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

    def __init__(self,readonly=False):

        # get the config file -- this must exist
        home=os.getenv("HOME")
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

        # set up an ssh tunnel if not running locally
        self.usetunnel=False
        self.hostname=socket.gethostname()
        if self.hostname=='lofar-server':
            self.con = mdb.connect('127.0.0.1', 'survey_user', self.password, 'surveys')
        else:
            try:
                dummy=socket.gethostbyname('lofar-server.data')
            except:
                self.usetunnel=True

            if self.usetunnel:
                self.tunnel=sshtunnel.SSHTunnelForwarder('lofar.herts.ac.uk',
                                                         ssh_username=self.ssh_user,
                                                         ssh_pkey=home+'/.ssh/%s'%self.ssh_key,
                                                         remote_bind_address=('127.0.0.1',3306),
                                                         local_bind_address=('127.0.0.1',))

                self.tunnel.start()
                localport=self.tunnel.local_bind_port
                self.con = mdb.connect('127.0.0.1', 'survey_user', self.password, 'surveys', port=localport)
            else:
                self.con = mdb.connect('lofar-server.data', 'survey_user', self.password, 'surveys')
        self.cur = self.con.cursor(cursorclass=mdbcursors.DictCursor)
        if self.readonly:
            pass
            #can't use this feature on lofar's version of MariaDB
            #self.cur.execute('set session transaction read only')
        else:
            self.cur.execute('lock table fields write, observations write, transients write')
        self.closed=False

    def close(self):
        if not self.closed:
            if not self.readonly:
                self.cur.execute('unlock tables')
            self.con.close()
            if self.usetunnel:
                self.tunnel.stop()
            self.closed=True # prevent del from trying again
    
    def __del__(self):
        self.close()
        
    def get_field(self,id):
        self.cur.execute('select * from fields where id=%s',(id,))
        result=self.cur.fetchall()
        if len(result)==0:
            return None
        else:
            return result[0]

    def set_field(self,sd):
        assert not self.readonly
        id=sd['id'];
        for k in sd:
            if k=='id':
                continue
            if sd[k] is not None:
                self.cur.execute('update fields set '+k+'=%s where id=%s',(sd[k],id))

    def create_field(self,id):
        self.cur.execute('insert into fields(id) values (%s)',(id,))
        return self.get_field(id)

    def get_observation(self,id):
        self.cur.execute('select * from observations where id=%s',(id,))
        result=self.cur.fetchall()
        if len(result)==0:
            return None
        else:
            return result[0]

    def set_observation(self,sd):
        assert not self.readonly
        id=sd['id'];
        for k in sd:
            if k=='id':
                continue
            if sd[k] is not None:
                self.cur.execute('update observations set '+k+'=%s where id=%s',(sd[k],id))

    def create_observation(self,id):
        self.cur.execute('insert into observations(id) values (%s)',(id,))
        return self.get_field(id)

    def get_transient(self,id):
        self.cur.execute('select * from transients where id=%s',(id,))
        result=self.cur.fetchall()
        if len(result)==0:
            return None
        else:
            return result[0]

    def set_transient(self,sd):
        assert not self.readonly
        id=sd['id'];
        for k in sd:
            if k=='id':
                continue
            if sd[k] is not None:
                self.cur.execute('update transients set '+k+'=%s where id=%s',(sd[k],id))

    def create_transient(self,id):
        self.cur.execute('insert into transients(id) values (%s)',(id,))
        return self.get_transient(id)

if __name__=='__main__':
    sdb=SurveysDB()
    result=sdb.get_field('P35Hetdex10')
    #result['location']='Never Never Land'
    #sdb.set_id(result)
    print result
    sdb.close()

