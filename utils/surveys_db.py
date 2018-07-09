import sshtunnel
import socket
import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors
import os

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
    result=sdb.get_id(123456)
    #result['location']='Never Never Land'
    #sdb.set_id(result)
    print result
    sdb.close()
