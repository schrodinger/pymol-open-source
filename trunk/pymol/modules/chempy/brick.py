
from Numeric import *

class Brick:
   
   def __init__(self):
      self.valid = None

   def setup_from_min_max(self,mn,mx,grid,buffer=0.0):
      self.origin = [
         mn[0]-buffer,
         mn[1]-buffer,
         mn[2]-buffer
         ]
      self.range = [
         (mx[0]-mn[0])+2*buffer,
         (mx[1]-mn[1])+2*buffer,
         (mx[2]-mn[2])+2*buffer
         ]
      self.dim = [
         1 + int(self.range[0]/grid[0]),
         1 + int(self.range[1]/grid[1]),
         1 + int(self.range[2]/grid[2])
         ]
      self.grid = copy.deepcopy(grid)
      self.lvl = zeros(self.dim,Float)
      
      

