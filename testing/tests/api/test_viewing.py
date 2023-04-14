from pymol import cmd


def test_clip_set_near_far():
    near_value = 10
    far_value = 20
    cmd.clip('near_set', near_value)
    cmd.clip('far_set', far_value)
    (near, far) = cmd.get_clip()
    assert near == near_value
    assert far == far_value


def test_clip_set_near_far_overlap():
    near_value = 10
    far_value = 8
    cmd.clip('near_set', near_value)
    cmd.clip('far_set', far_value)
    (near, far) = cmd.get_clip()
    # indeterminate? but adjusts itself for a distance of 1
    assert far - near == 1
