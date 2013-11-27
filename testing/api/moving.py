
from pymol import cmd, testing, stored

class TestMoving(testing.PyMOLTestCase):

    def prep_movie(self):
        cmd.mset("1x60")
    
    def testAccept(self):
        cmd.accept
        self.skipTest("TODO")

    def testBackward(self):
        self.prep_movie()
        cmd.middle()
        cmd.backward()
        self.assertEquals(cmd.get_frame(),30)

        cmd.frame(1)
        cmd.backward()
        self.assertEquals(cmd.get_frame(),1)
        

    def testDecline(self):
        cmd.decline
        self.skipTest("TODO")

    def testEnding(self):
        self.prep_movie()
        cmd.ending()
        self.assertEquals(cmd.get_frame(),60)

    def testForward(self):
        self.prep_movie()
        cmd.forward()
        self.assertEquals(cmd.get_frame(),2)
        cmd.forward()
        self.assertEquals(cmd.get_frame(),3)
        cmd.frame(30)
        self.assertEquals(cmd.get_frame(),30)
        cmd.frame(60)
        cmd.forward()
        self.assertEquals(cmd.get_frame(),60)

        cmd.frame(60)
        cmd.forward()
        self.assertEquals(cmd.get_frame(),60)

    def testFrame(self):
        self.prep_movie()
        cmd.frame(30)
        self.assertEquals(cmd.get_frame(),30)

        cmd.frame(60)
        self.assertEquals(cmd.get_frame(),60)

        cmd.frame(65)
        self.assertEquals(cmd.get_frame(),60)

        cmd.frame(-5)
        self.assertEquals(cmd.get_frame(),1)

        cmd.frame(10)
        cmd.frame(1)
        self.assertEquals(cmd.get_frame(),1)


    def testGetFrame(self):
        self.prep_movie()
        cmd.frame(30)
        self.assertEquals(cmd.get_frame(),30)

        cmd.frame(60)
        self.assertEquals(cmd.get_frame(),60)

        cmd.frame(65)
        self.assertEquals(cmd.get_frame(),60)

        cmd.frame(-5)
        self.assertEquals(cmd.get_frame(),1)

        cmd.frame(10)
        cmd.frame(1)
        self.assertEquals(cmd.get_frame(),1)

    def testGetMoviePlaying(self):
        self.prep_movie()
        cmd.mplay()
        self.assertEquals(cmd.get_movie_playing(),1)

        cmd.mstop()
        self.assertEquals(cmd.get_movie_playing(),0)
        

    def testGetState(self):
        cmd.get_state
        self.skipTest("TODO")

    def testMadd(self):
        # test clearing
        self.prep_movie()
        self.assertEquals(cmd.count_frames(),60)

        # test states and frame mapping
        cmd.madd("1x10 1 -10 10 -1")
        self.assertEquals(cmd.count_frames(),90)
        
        cmd.frame(1)
        self.assertEquals(cmd.get_frame(),1)
        self.assertEquals(cmd.get("state"),"1")
        cmd.frame(10)
        self.assertEquals(cmd.get_frame(),10)
        self.assertEquals(cmd.get("state"),"1")
        cmd.frame(72)
        self.assertEquals(cmd.get_frame(),72)
        self.assertEquals(cmd.get("state"),"2")
        cmd.frame(82)
        self.assertEquals(cmd.get_frame(),82)
        self.assertEquals(cmd.get("state"),"9")
        cmd.frame(90)
        self.assertEquals(cmd.get_frame(),90)
        self.assertEquals(cmd.get("state"),"1")

    def testMappend(self):
        self.prep_movie()
        cmd.fragment("ala")
        cmd.mappend(5,"delete *")
        cmd.mappend(5,"frame 1")
        cmd.frame(5)
        self.assertEquals(len(cmd.get_names()),0)
        self.assertEquals(cmd.get_frame(),1)

        
    def testMclear(self):
        cmd.mclear
        self.skipTest("TODO")

    def testMcopy(self):
        cmd.mcopy
        self.skipTest("TODO")

    def testMdelete(self):
        cmd.mset("1x10 1 -10 10 -1")

        # 1-10 = state 1
        # 11-20 = states 1..10
        # 21-30 = states 10..1

        cmd.mdelete(count=10, frame=11)
        # simple delete
        self.assertEquals(cmd.count_frames(),20)

        # test state mapping
        cmd.frame(1)
        self.assertEquals(cmd.get("state"),"1")
        cmd.frame(11)
        self.assertEquals(cmd.get("state"),"10")
        cmd.frame(20)
        self.assertEquals(cmd.get("state"),"1")

    def testMdo(self):
        self.prep_movie()
        cmd.fragment("ala")
        cmd.mdo(5,"delete *")
        cmd.frame(5)
        self.assertEquals(len(cmd.get_names()),0)

    def testMdump(self):
        cmd.mdump
        self.skipTest("TODO: This prints to console, not sure how to test.")

    def testMiddle(self):
        self.prep_movie()
        cmd.middle()
        self.assertEquals(cmd.get_frame(),31)

    def testMinsert(self):
        self.prep_movie()
        cmd.minsert(20)
        self.assertEquals(cmd.count_frames(),80)

        # need a more sophisticated test for this
        

    def testMmatrix(self):
        cmd.mmatrix
        self.skipTest("TODO")

    def testMmove(self):
        cmd.mmove
        self.skipTest("TODO")

    def testMplay(self):
        self.prep_movie()
        cmd.mplay()
        self.assertEquals(cmd.get_movie_playing(),1)
        cmd.mstop()

    def testMpng(self):
        cmd.mpng
        self.skipTest("TODO")

    def testMset(self):
        # basic tet
        self.prep_movie()
        self.assertEquals(cmd.count_frames(),60)

        # test clearing
        cmd.mset()
        self.assertEquals(cmd.count_frames(),0)

        # test states and frame mapping
        cmd.mset("1x10 1 -10 10 -1")
        cmd.frame(1)
        self.assertEquals(cmd.get_frame(),1)
        self.assertEquals(cmd.get("state"),"1")
        cmd.frame(10)
        self.assertEquals(cmd.get_frame(),10)
        self.assertEquals(cmd.get("state"),"1")
        cmd.frame(12)
        self.assertEquals(cmd.get_frame(),12)
        self.assertEquals(cmd.get("state"),"2")
        cmd.frame(22)
        self.assertEquals(cmd.get_frame(),22)
        self.assertEquals(cmd.get("state"),"9")
        cmd.frame(30)
        self.assertEquals(cmd.get_frame(),30)
        self.assertEquals(cmd.get("state"),"1")

    def testMstop(self):
        self.prep_movie()
        cmd.mplay()
        cmd.mstop()
        self.assertEquals(cmd.get_movie_playing(),0)

    def testMtoggle(self):
        self.prep_movie()
        cmd.mstop()
        
        # toggle on
        cmd.mtoggle()
        self.assertEquals(cmd.get_movie_playing(),1)
        # toggle off
        cmd.mtoggle()
        self.assertEquals(cmd.get_movie_playing(),0)


    def testMview(self):
        cmd.mview
        self.skipTest("TODO")

    def testRewind(self):
        self.prep_movie()
        cmd.frame(30)
        self.assertEquals(cmd.get_frame(), 30)

        cmd.rewind()
        self.assertEquals(cmd.get_frame(), 1)

        
    def testSetFrame(self):

        self.prep_movie()

        # within
        cmd.set_frame(30)
        self.assertEquals(cmd.get_frame(),30)

        # extent < 
        cmd.set_frame(-1)
        self.assertEquals(cmd.get_frame(),1)

        # extent > 
        cmd.set_frame(300)
        self.assertEquals(cmd.get_frame(),60)

        
