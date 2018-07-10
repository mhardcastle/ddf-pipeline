import sshtunnel
import socket
import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors
import os

def get_next():
    # return the name of the top-priority field with status 'Preprocessed'
    sdb=SurveysDB()
    sdb.cur.execute('select * from fields where status="Preprocessed" order by priority desc')
    results=sdb.cur.fetchall()
    sdb.close()
    return results[0]['id']

def update_status(name,status):
    # utility function to just update the status of an observation
    # name can be None (work it out from cwd), string (strip L) or int (use int)

    if name is None:
        # work it out
        id=get_id()
    else:
        if isinstance(name,str):
            id=int(name[1:])
        else:
            id=int(name)

    sdb=SurveysDB()
    idd=sdb.get_id(id)
    idd['status']=status
    tag_idd(sdb,idd)
    sdb.set_id(idd)
    sdb.close()

def tag_idd(sdb,idd):
    # Add location and user tags
    idd['clustername']=get_cluster()
    idd['location']=os.getcwd()
    idd['username']=get_user()
    idd['nodename']=sdb.hostname
    
def get_id():
    dir=os.getcwd()
    dname=dir.split('/')[-1]
    return int(dname[1:])

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

    def __init__(self,localport=33306):

        # get the config file -- this must exist
        home=os.getenv("HOME")
        cfg=open(home+'/.surveys').readlines()
        self.password=cfg[0].rstrip()
        try:
            self.ssh_user=cfg[1].rstrip()
        except:
            self.ssh_user=None
        
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
                                                         ssh_pkey=home+'/.ssh/id_rsa',
                                                         remote_bind_address=('127.0.0.1',3306),
                                                         local_bind_address=('127.0.0.1',localport))
                self.tunnel.start()
                self.con = mdb.connect('127.0.0.1', 'survey_user', self.password, 'surveys', port=localport)
            else:
                self.con = mdb.connect('lofar-server.data', 'survey_user', self.password, 'surveys')
        self.cur = self.con.cursor(cursorclass=mdbcursors.DictCursor)
        self.cur.execute('lock table fields write')
        self.closed=False

    def close(self):
        if not self.closed:
            self.cur.execute('unlock tables')
            self.con.close()
            if self.usetunnel:
                self.tunnel.stop()
            self.closed=True # prevent del from trying again
    
    def __del__(self):
        self.close()
        

    def get_id(self,id):
        self.cur.execute('select * from fields where id=%s',(id,))
        result=self.cur.fetchall()
        if len(result)==0:
            return None
        else:
            return result[0]

    def set_id(self,sd):
        id=sd['id'];
        for k in sd:
            if k=='id':
                continue
            if sd[k] is not None:
                self.cur.execute('update fields set '+k+'=%s where id=%s',(sd[k],id))

    def create_id(self,id):
        self.cur.execute('insert into fields(id) values (%s)',(id,))
        return self.get_id(id)

if __name__=='__main__':
    sdb=SurveysDB()
    result=sdb.get_id(647109)
    #result['location']='Never Never Land'
    #sdb.set_id(result)
    print result
    sdb.close()

