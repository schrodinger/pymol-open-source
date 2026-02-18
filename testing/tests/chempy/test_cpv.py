from unittest.mock import patch
import random

import pytest
import numpy as np

from chempy import cpv


@pytest.fixture
def vector_a() -> list[float]:
    return [3.0, 4.0, 0.0]


@pytest.fixture
def vector_b() -> list[float]:
    return [6.0, 8.0, 0.0]


@pytest.fixture
def matrix_a() -> list[list[float]]:
    return [[4, 8, 15], [16, 23, 42], [1984, 1, 10]]


@pytest.fixture
def matrix_b() -> list[list[float]]:
    return [[5, 1, 3], [1, 1, 1], [1, 2, 1]]


def test_get_null():
    assert cpv.get_null() == np.zeros(shape=3, dtype=np.float64).tolist()
    assert cpv.get_null() == np.array([0, 0, 0], dtype=np.float64).tolist()


def test_get_identity():
    assert cpv.get_identity() == np.identity(3, dtype=np.float64).tolist()
    assert (
        cpv.get_identity()
        == np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64).tolist()
    )


def test_distance_sq(vector_a, vector_b):

    assert (
        cpv.distance_sq(vector_a, vector_b)
        == np.linalg.norm(np.array(vector_a) - np.array(vector_b)) ** 2
    )
    assert cpv.distance_sq([3, 3, 3], [2, 2, 2]) == 3.0


def test_distance(vector_a, vector_b):
    assert cpv.distance(vector_a, vector_b) == np.linalg.norm(
        np.array(vector_a) - np.array(vector_b)
    )
    assert cpv.distance(vector_a, vector_b) == 5.0


def test_length(vector_a):
    assert cpv.length(vector_a) == np.linalg.norm(vector_a)
    assert cpv.length(vector_a) == 5.0


def test_random_displacement(vector_a):
    radius = 2

    with (
        patch.object(cpv, "_random_vector", return_value=[0.1, 0.1, 0.1]),
        patch("numpy.random.random", return_value=np.array([0.1, 0.1, 0.1])),
        patch("random.random", return_value=1),
    ):
        rand_np_vector = np.random.random(3)
        v_len = random.random() * radius / np.linalg.norm(rand_np_vector)
        assert (
            cpv.random_displacement(vector_a, radius)
            == (np.array(vector_a) + rand_np_vector * v_len).tolist()
        )
        assert (
            cpv.random_displacement(vector_a, radius)
            == np.array(
                [4.1547005383792515, 5.1547005383792515, 1.1547005383792512]
            ).tolist()
        )


def test_random_sphere(vector_a):
    radius = 2

    with (
        patch.object(cpv, "_random_vector", return_value=[0.1, 0.1, 0.1]),
        patch("numpy.random.random", return_value=np.array([0.1, 0.1, 0.1])),
    ):
        rand_np_vector = np.random.random(3)
        assert (
            cpv.random_sphere(vector_a, radius)
            == (
                np.array(vector_a)
                + rand_np_vector * (2 * radius / np.linalg.norm(rand_np_vector))
            ).tolist()
        )
        assert (
            cpv.random_sphere(vector_a, radius)
            == np.array(
                [
                    5.309401076758503,
                    6.309401076758503,
                    2.3094010767585025,
                ]
            ).tolist()
        )


def test_random_vector():
    with (
        patch.object(cpv, "_random_vector", return_value=[0.1, 0.1, 0.1]),
        patch("numpy.random.random", return_value=np.array([0.1, 0.1, 0.1])),
    ):
        assert cpv.random_vector() == (np.random.random(3) * 2).tolist()
        assert cpv.random_vector() == np.array([0.2, 0.2, 0.2]).tolist()


def test_add(vector_a, vector_b):
    assert (
        cpv.add(vector_a, vector_b)
        == (np.array(vector_a) + np.array(vector_b)).tolist()
    )
    assert cpv.add(vector_a, vector_b) == np.array([9.0, 12.0, 0.0]).tolist()


def test_sub(vector_a, vector_b):
    assert (
        cpv.sub(vector_a, vector_b)
        == (np.array(vector_a) - np.array(vector_b)).tolist()
    )
    assert cpv.sub(vector_a, vector_b) == (np.array([-3.0, -4.0, 0.0]).tolist())


def test_average(vector_a, vector_b):
    assert (
        cpv.average(vector_a, vector_b)
        == np.average([np.array(vector_a), np.array(vector_b)], axis=0).tolist()
    )
    assert cpv.average(vector_a, vector_b) == np.array([4.5, 6.0, 0.0]).tolist()


def test_scale(vector_a):
    assert cpv.scale(vector_a, 2) == (np.array(vector_a) * 2).tolist()
    assert cpv.scale(vector_a, 2) == np.array([6.0, 8.0, 0.0]).tolist()


def test_negate(vector_a):
    assert cpv.negate(vector_a) == (np.array(vector_a) * -1).tolist()
    assert cpv.negate(vector_a) == np.array([-3.0, -4.0, -0.0]).tolist()


def test_dot_product(vector_a, vector_b):
    assert cpv.dot_product(vector_a, vector_b) == np.dot(
        np.array(vector_a), np.array(vector_b)
    )
    assert cpv.dot_product(vector_a, vector_b) == 50.0


def test_cross_product(vector_a, vector_b):
    assert (
        cpv.cross_product(vector_a, vector_b)
        == np.cross(np.array(vector_a), np.array(vector_b)).tolist()
    )
    assert cpv.cross_product(vector_a, vector_b) == np.array([0.0, 0.0, 0.0]).tolist()


def test_transform(matrix_a, vector_a):
    assert (
        cpv.transform(matrix_a, vector_a)
        == np.dot(np.array(matrix_a), np.array(vector_a)).tolist()
    )
    assert cpv.transform(matrix_a, vector_a) == np.array([44.0, 140.0, 5956.0]).tolist()


def test_inverse_transform(matrix_a, vector_a):
    assert (
        cpv.inverse_transform(matrix_a, vector_a)
        == np.dot(np.array(matrix_a).T, np.array(vector_a)).tolist()
    )
    assert (
        cpv.inverse_transform(matrix_a, vector_a)
        == np.array([76.0, 116.0, 213.0]).tolist()
    )


def test_multiply(matrix_a, matrix_b):
    assert (
        cpv.multiply(matrix_a, matrix_b)
        == np.dot(np.array(matrix_a), np.array(matrix_b)).tolist()
    )
    assert (
        cpv.multiply(matrix_a, matrix_b)
        == np.array([[43, 42, 35], [145, 123, 113], [9931, 2005, 5963]]).tolist()
    )


def test_transpose(matrix_a):
    assert cpv.transpose(matrix_a) == np.array(matrix_a).T.tolist()
    assert (
        cpv.transpose(matrix_a)
        == np.array([[4, 16, 1984], [8, 23, 1], [15, 42, 10]]).tolist()
    )


def test_get_system2(vector_a, vector_b):
    def normalize(v: np.ndarray) -> np.ndarray:
        return (
            np.array(v) / vector_len
            if (vector_len := np.linalg.norm(np.array(v))) > cpv.RSMALL4
            else np.zeros(shape=3, dtype=np.float64)
        )

    z = normalize(np.cross(np.array(vector_a), np.array(vector_b)))
    y = normalize(np.cross(z, np.array(vector_a)))
    x = normalize(vector_a)
    assert cpv.get_system2(vector_a, vector_b) == np.array([x, y, z]).tolist()
    assert (
        cpv.get_system2(vector_a, vector_b)
        == np.array([[0.6, 0.8, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]).tolist()
    )


# ------------------------------------------------------------------------------
def test_scale_system(matrix_a):
    assert cpv.scale_system(matrix_a, 1984) == (np.array(matrix_a) * 1984).tolist()
    assert (
        cpv.scale_system(matrix_a, 1985)
        == np.array(
            [[7940, 15880, 29775], [31760, 45655, 83370], [3938240, 1985, 19850]]
        ).tolist()
    )


def test_transform_about_point(matrix_a, vector_a, vector_b):
    assert (
        cpv.transform_about_point(matrix_a, vector_a, vector_b)
        == (
            np.dot(np.array(matrix_a), np.array(vector_a) - np.array(vector_b))
            + np.array(vector_b)
        ).tolist()
    )
    assert (
        cpv.transform_about_point(matrix_a, vector_a, vector_b)
        == np.array([-38.0, -132.0, -5956.0]).tolist()
    )


def test_get_angle(vector_a, vector_b):
    np_result = (
        np.dot(np.array(vector_a), np.array(vector_b)) / denom
        if (
            denom := np.linalg.norm(np.array(vector_a))
            * np.linalg.norm(np.array(vector_b))
        )
        > 1e-10
        else 0.0
    )

    assert cpv.get_angle(vector_a, vector_b) == np.arccos(np_result)
    assert cpv.get_angle(vector_a, vector_b) == 0.0


def test_get_angle_formed_by(vector_a, vector_b):
    r1 = np.linalg.norm(np.array(vector_a) - np.array(vector_b))
    r2 = np.linalg.norm(np.array(vector_b) - np.array(vector_b))
    r3 = np.linalg.norm(np.array(vector_a) - np.array(vector_b))

    theta = (
        np.pi
        if (r1 + r2 - r3) < 1.0e-10
        else np.arccos((r1**2 + r2**2 - r3**2) / (2.0 * r1 * r2))
    )

    assert cpv.get_angle_formed_by(vector_a, vector_b, vector_b) == theta
    assert cpv.get_angle_formed_by(vector_a, vector_b, vector_b) == np.pi


def test_project(vector_a, vector_b):
    np_version_dot = np.dot(np.array(vector_a), np.array(vector_b))

    assert (
        cpv.project(vector_a, vector_b)
        == (np.array(vector_b) * np_version_dot).tolist()
    )
    assert cpv.project(vector_a, vector_b) == np.array([300.0, 400.0, 0.0]).tolist()


def test_remove_component(vector_a, vector_b):
    np_version_dot = np.dot(np.array(vector_a), np.array(vector_b))

    assert (
        cpv.remove_component(vector_a, vector_b)
        == (np.array(vector_a) - np.array(vector_b) * np_version_dot).tolist()
    )
    assert (
        cpv.remove_component(vector_a, vector_b)
        == np.array([-297.0, -396.0, 0.0]).tolist()
    )


def test_normalize(vector_a):
    np_result = (
        np.array(vector_a) / vector_len
        if (vector_len := np.linalg.norm(np.array(vector_a))) > cpv.RSMALL4
        else np.zeros(shape=3, dtype=np.float64)
    )

    assert cpv.normalize(vector_a) == np_result.tolist()
    assert cpv.normalize(vector_a) == np.array([0.6, 0.8, 0.0]).tolist()


def test_normalize_failsafe(vector_a):
    np_result = (
        np.array(vector_a) / vector_len
        if (vector_len := np.linalg.norm(np.array(vector_a))) > cpv.RSMALL4
        else np.array([1.0, 0.0, 0.0])
    )

    assert cpv.normalize_failsafe(vector_a) == np_result.tolist()
    assert cpv.normalize_failsafe(vector_a) == np.array([0.6, 0.8, 0.0]).tolist()


def test_rotation_matrix():
    angle = 30
    axis = [4, 8, 15]
    axis = np.array(axis)
    s = np.sin(angle)
    c = np.cos(angle)

    mag = np.linalg.norm(axis)

    if abs(mag) < cpv.RSMALL4:
        return np.eye(3)

    axis = axis / mag

    # Rodrigues' rotation formula
    R = (
        c * np.eye(3)
        + (1 - c) * np.outer(axis, axis)
        + s
        * np.array(
            [[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]]
        )
    )

    assert cpv.rotation_matrix(angle, [4, 8, 15]) == R.tolist()
    assert (
        cpv.rotation_matrix(angle, [4, 8, 15])
        == np.array(
            [
                [0.1986185869426616, 0.9373521674176533, -0.28621944580745806],
                [-0.759883619197343, 0.3317199981078943, 0.5590516327950812],
                [0.6189729737205398, 0.1064554230310823, 0.7781643147246124],
            ],
        ).tolist()
    )


def test_transform_array(matrix_a, matrix_b):

    assert (
        cpv.transform_array(matrix_a, matrix_b)
        == np.apply_along_axis(
            lambda vertex: np.dot(np.array(matrix_a), np.array(vertex)),
            1,
            np.array(matrix_b),
        ).tolist()
    )
    assert (
        cpv.transform_array(matrix_a, matrix_b)
        == np.array(
            [
                [73, 229, 9951],
                [27, 81, 1995],
                [35, 104, 1996],
            ]
        ).tolist()
    )


def test_translate_array(vector_a, matrix_a):
    assert (
        cpv.translate_array(vector_a, matrix_a)
        == np.apply_along_axis(
            lambda vertex: np.array(vector_a) + np.array(vertex), 1, np.array(matrix_a)
        ).tolist()
    )
    assert (
        cpv.translate_array(vector_a, matrix_a)
        == np.array(
            [[7.0, 12.0, 15.0], [19.0, 27.0, 42.0], [1987.0, 5.0, 10.0]]
        ).tolist()
    )


def test_fit(matrix_a, matrix_b):

    def fit(
        target_array: np.ndarray, source_array: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
        """fit(target_array, source_array) -> (t1, t2, rot_mtx, rmsd) [fit_result]

        Calculates the translation vectors and rotation matrix required
        to superimpose source_array onto target_array.  Original arrays are
        not modified.  NOTE: Currently assumes 3-dimensional coordinates
        """
        # Check dimensions of input arrays
        if target_array.shape != source_array.shape:
            raise ValueError("Error: arrays must be of same length for RMS fitting.")
        if target_array.shape[1] != 3 or source_array.shape[1] != 3:
            raise ValueError("Error: arrays must be dimension 3 for RMS fitting.")

        ndim = 3
        maxiter = 200
        tol = 0.001

        # Calculate translation vectors (center-of-mass)
        t1 = np.mean(target_array, axis=0)
        t2 = np.mean(source_array, axis=0)

        # Calculate correlation matrix
        tvec1 = target_array - t1
        tvec2 = source_array - t2
        corr_mtx = np.dot(tvec2.T, tvec1)

        # Initial rotation matrix (identity matrix)
        rot_mtx = np.eye(ndim)

        # Main iteration scheme (hardwired for 3X3 matrix)
        for iters in range(maxiter):
            iy = iters % ndim
            iz = (iters + 1) % ndim
            sig = corr_mtx[iz, iy] - corr_mtx[iy, iz]
            gam = corr_mtx[iy, iy] + corr_mtx[iz, iz]

            sg = np.sqrt(sig**2 + gam**2)

            if sg != 0.0 and abs(sig) > tol * abs(gam):
                sg = 1.0 / sg
                for i in range(ndim):

                    bb = gam * corr_mtx[iy, i] + sig * corr_mtx[iz, i]
                    cc = gam * corr_mtx[iz, i] - sig * corr_mtx[iy, i]
                    corr_mtx[iy, i] = bb * sg
                    corr_mtx[iz, i] = cc * sg

                    bb = gam * rot_mtx[iy, i] + sig * rot_mtx[iz, i]
                    cc = gam * rot_mtx[iz, i] - sig * rot_mtx[iy, i]
                    rot_mtx[iy, i] = bb * sg
                    rot_mtx[iz, i] = cc * sg

                continue

        # Calculate RMS deviation
        vt = np.dot(source_array - t2, rot_mtx.T) + t1
        rmsd = np.sqrt(np.mean(np.sum((target_array - vt) ** 2, axis=1)))

        return t1, t2, rot_mtx, rmsd

    t1_np, t2_np, rot_mtx_np, rmsd_np = fit(np.array(matrix_a), np.array(matrix_b))
    result = cpv.fit(matrix_a, matrix_b)

    if result is None:
        raise
    t1, t2, rot_mtx, rmsd = result
    assert t1 == t1_np.tolist()
    assert t2 == t2_np.tolist()
    for i in range(len(rot_mtx)):
        assert rot_mtx[i] == pytest.approx(rot_mtx_np[i].tolist(), 0.001)
    assert rmsd == pytest.approx(rmsd_np)


def test_fit_apply(matrix_a, matrix_b):
    result = cpv.fit(matrix_a, matrix_b)

    if result is None:
        raise
    t1, mt2, m = result[:3]

    origin_result = cpv.fit_apply(result, matrix_b)
    np_result = np.array(t1) + np.dot(np.array(m), (np.array(mt2) + np.array(matrix_b)))

    for i in range(len(origin_result)):
        assert origin_result[i] == pytest.approx(np_result[i].tolist(), 1)

    np.testing.assert_array_almost_equal(
        np.array(origin_result),
        np.array(
            [
                [661.0827157827758, 7.129705092201009, 17.79008759059663],
                [665.1572144529815, 8.019304893172958, 19.404730826919444],
                [665.5694270238597, 7.575498193306849, 18.609044701890262],
            ]
        ).tolist()
    )
    np.testing.assert_array_almost_equal(
        np_result,
        np.array(
            [
                [662.0405274477379, 8.365106579876256, 18.542585300247588],
                [665.1798712452187, 7.2261452638093076, 19.82182248600216],
                [662.2963554892253, 8.460848531966128, 18.41893138496441],
            ]
        ),
    )
