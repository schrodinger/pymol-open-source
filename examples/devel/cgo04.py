from pymol.cgo import *

cgo = [

   BEGIN, LINES,
   COLOR,   1,1,1,
   VERTEX, 0,0,0,
   COLOR,   0,0,1,
   VERTEX, 0,0,1,
   COLOR,   1,1,1,
   VERTEX, 0,0,0,
   COLOR,   0,1,0,
   VERTEX, 0,1,0,
   COLOR,   1,1,1,
   VERTEX, 0,0,0,
   COLOR,   1,0,0,
   VERTEX, 1,0,0,
   END,

# nitrogen
   COLOR, 0.25, 0.25, 1.00,
   SPHERE, 0.25, 0.25, 1.00, 0.12,

# oxygen
   COLOR, 1.00, 0.20, 0.20,
   SPHERE, 1.00, 0.20, 0.20, 0.12,

# sulfur
   COLOR, 1.00, 0.50, 0.00,
   SPHERE, 1.00, 0.50, 0.00, 0.12,
]

set1 = [

# carbon (GREEN)
   COLOR, 0.20, 1.00, 0.20,
   SPHERE, 0.20, 1.00, 0.20, 0.12,

# cyan (CYAN)
   COLOR, 0.00, 1.00, 1.00,   
   SPHERE, 0.00, 1.00, 1.00, 0.10,  

# magenta (MAGENTA)
   COLOR, 1.00, 0.20, 0.80,   
   SPHERE, 1.00, 0.20, 0.80, 0.10,  

# yellow (YELLOW)
   COLOR, 1.00, 1.00, 0.00,   
   SPHERE, 1.00, 1.00, 0.00, 0.10,

# salmon (SALMON)
   COLOR, 1.00, 0.60, 0.60,
   SPHERE, 1.00, 0.60, 0.60, 0.10,

# gray90 (WHITE)
   COLOR, 0.90, 0.90, 0.90,
   SPHERE, 0.90, 0.90, 0.90, 0.10,

# slate (BLUE)
   COLOR, 0.50, 0.50, 1.00,
   SPHERE, 0.50, 0.50, 1.00, 0.10, 

# brightorange (ORANGE)
   COLOR, 1.00, 0.70, 0.20,
   SPHERE, 1.00, 0.70, 0.20, 0.10,
]

set2 = [

# smudge (GREEN) *NEW
   COLOR, 0.55, 0.7, 0.40,
   SPHERE, 0.55, 0.7, 0.40, 0.06,

# deepteal (CYAN) 
   COLOR, 0.1, 0.60, 0.60,
   SPHERE, 0.1, 0.60, 0.60, 0.06,

# hotpink (MAGENTA) 
   COLOR, 1.00, 0.00, 0.50,
   SPHERE, 1.00, 0.00, 0.50, 0.06,

# yelloworange (ORANGE) *NEW
   COLOR, 1.00, 0.87, 0.37,
   SPHERE, 1.00, 0.87, 0.37, 0.04,

# violetpurple (SALMON) *NEW
   COLOR, 0.55, 0.25, 0.60,
   SPHERE, 0.55, 0.25, 0.60, 0.06,

# grey70 (WHITE)
   COLOR, 0.70, 0.70, 0.70,
   SPHERE, 0.70, 0.70, 0.70, 0.06,

# marine (BLUE)
   COLOR, 0.00, 0.50, 1.00,
   SPHERE, 0.00, 0.50, 1.00, 0.03,
   
# olive (ORANGE)
   COLOR, 0.77, 0.70, 0.0,   
   SPHERE, 0.77, 0.70, 0.0, 0.06,  

]


set3  = [
# lime (GREEN)
   COLOR, 0.50, 1.00, 0.50,
   SPHERE, 0.50, 1.00, 0.50, 0.04,

# aquamarine (CYAN)
   COLOR, 0.50, 1.00, 1.00,
   SPHERE, 0.50, 1.00, 1.00, 0.04,

# dirtyviolet (MAGENTA) *NEW
   COLOR, 0.70, 0.50, 0.70,   
   SPHERE, 0.70, 0.50, 0.70, 0.04,  

# wheat (YELLOW)
   COLOR, 0.99, 0.82, 0.65,
   SPHERE, 0.99, 0.82, 0.65, 0.06,

# deepsalmon (SALMON) *NEW
   COLOR, 1.00, 0.42, 0.42,
   SPHERE, 1.00, 0.42, 0.42, 0.04,

# lightpink (WHITE) *NEW
  COLOR, 1.00, 0.75, 0.87,
  SPHERE, 1.0, 0.75, 0.87, 0.04,

# teal (BLUE)
   COLOR, 0.00, 0.75, 0.75,
   SPHERE, 0.00, 0.75, 0.75, 0.04,

# limon (ORANGE) *NEW
   COLOR, 0.75, 1.00, 0.25,   
   SPHERE, 0.75, 1.00, 0.25, 0.08,

]

set4  = [

# bluegreen (GREEN) 

   COLOR, 0.0, 1.0, 0.50,
   SPHERE, 0.0, 1.0, 0.50, 0.08,

# skyblue (CYAN) *NEW
   COLOR, 0.20, 0.50, 0.80,
   SPHERE, 0.20, 0.50, 0.80, 0.08,

# warmpink (MAGENTA) *NEW
   COLOR, 0.85, 0.20, 0.50,
   SPHERE, 0.85, 0.20, 0.50, 0.08,

# paleyellow (YELLOW) 
   COLOR, 1.00, 1.00, 0.50,
   SPHERE, 1.00, 1.00, 0.50, 0.08,

# violet (SALMON)
   COLOR, 1.00, 0.50, 1.0,   
   SPHERE, 1.00, 0.50, 1.0, 0.03,  

# bluewhite (WHITE) *NEW
   COLOR, 0.85, 0.85, 1.0,
   SPHERE, 0.85, 0.85, 1.0, 0.08,

# greencyan (BLUE) *NEW
   COLOR, 0.25, 1.00, 0.75,   
   SPHERE, 0.25, 1.00, 0.75, 0.08,  

# sand (ORANGE) *NEW
   COLOR, 0.72, 0.55, 0.30,   
   SPHERE, 0.72, 0.55, 0.30, 0.08,

]

set5  = [
# forest (GREEN) *UPDATED
   COLOR, 0.20, 0.60, 0.20,
   SPHERE, 0.20, 0.60, 0.20, 0.03,

# lightteal (CYAN) *NEW
  COLOR, 0.40, 0.70, 0.70,
  SPHERE, 0.40, 0.70, 0.70, 0.03,

# darksalmon (MAGENTA) *NEW
   COLOR, 0.73, 0.55, 0.52,
   SPHERE, 0.73, 0.55, 0.52, 0.08,

# splitpea (YELLOW) *NEW
   COLOR, 0.52, 0.75, 0.00,   
   SPHERE, 0.52, 0.75, 0.00, 0.03,

# rasberry (SALMON) *NEW
   COLOR, 0.70, 0.30, 0.40,
   SPHERE, 0.70, 0.30, 0.40, 0.03,

# grey50 (WHITE)
   COLOR, 0.50, 0.50, 0.50,
   SPHERE, 0.50, 0.50, 0.50, 0.03,

# brown (ORANGE) *updated
   COLOR, 0.65, 0.32, 0.17,
   SPHERE, 0.65, 0.32, 0.17, 0.04,
]


cgo.extend(set1)
cgo.extend(set2)
cgo.extend(set3)
cgo.extend(set4)
cgo.extend(set5)

cgo.extend([
   COLOR,   1,1,1,
   SPHERE,   0,0,0,0.5,
])

cmd.load_cgo(cgo,"test")
cmd.set("cgo_dot_width",20)
