
import chempy
from chempy import io

path = chempy.path + 'fragments/'

def get(name):
    return io.pkl.fromFile(path+name+'.pkl')

