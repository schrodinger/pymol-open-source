
# NOTE that in general, PyMOL API calls should not be made during the
# __call__ method of a callback object.  The PYMOL API is not designed
# to be re-entrant, and the callback object's __call__ method may in
# fact be called as far of an API fuction (refresh, ray, etc.).

# There are a few exceptions, which you may call from a __call__ method
# in order to get useful information:
#
#    cmd.get_frame()
#    cmd.get_state()

class Callback(object):
    def __call__(self):
        pass

    def get_extent(self):
        # should return [ [min_x, min_y, min_z], [ max_x, max_y, max_z ] ]
        return None
