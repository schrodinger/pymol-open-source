#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific.
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information.
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-*
#-*
#-*
#Z* -------------------------------------------------------------------

normal = {
'HOH': {
   'atoms' : {
      'O'    : { 'symbol' : 'O' , } ,
      'H1'   : { 'symbol' : 'H' , } ,
      'H2'   : { 'symbol' : 'H' , } ,
      } ,
   'bonds' : {
      ( 'O'   , 'H1'   ) : { 'order' : 1 , } ,
      ( 'O'   , 'H2'   ) : { 'order' : 1 , } ,
      ( 'H1'  , 'H2'   ) : { 'order' : 1 , } ,
      } ,
   'aliases' : {
      'OH2'  : 'O'    ,
      } ,
   } ,
'WAT': {
   'atoms' : {
      'O'    : { 'symbol' : 'O' , } ,
      'H1'   : { 'symbol' : 'H' , } ,
      'H2'   : { 'symbol' : 'H' , } ,
      } ,
   'bonds' : {
      ( 'O'   , 'H1'   ) : { 'order' : 1 , } ,
      ( 'O'   , 'H2'   ) : { 'order' : 1 , } ,
      ( 'H1'  , 'H2'   ) : { 'order' : 1 , } ,
      } ,
   'aliases' : {
      'OH2'  : 'O'    ,
      } ,
   } ,
}
