import os
import fcntl
import termios
import struct

def ioctl_GWINSZ(fd):
    try:
        cr = struct.unpack('hh',
                           fcntl.ioctl(fd, termios.TIOCGWINSZ, '1234'))
        return cr
    except:
        pass

def get_terminal_size_linux():
    ''' From https://gist.github.com/jtriley/1108174 '''
    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    if not cr:
        try:
            cr = (os.environ['LINES'], os.environ['COLUMNS'])
        except:
            return None,None
    return int(cr[1]), int(cr[0])
