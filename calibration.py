import datetime
from pathlib import Path

import numpy as np
import scipy.io


class CameraCalibration(object):
    """Camera calibration saved in .mat file and method to assemble Projective (P) martrix.

    Notes:
        - Inspired by example code + notes from CiRC which are derived from Hartley and Zisserman (20030.)
        - Asssumes calibration saved in .mat file

    Args:
        calibration_file (str): Path to camera calibration file.

    Attributes:
        fname (str): Name of camera calibration file
        serial_number (int): Camera serial number
        camera_number (str): Camera number (e.g. 'c5')
        calibration_date (datetime): Date of camera calibration
        coordinate_system (str): Coordinate system used for calibration (e.g. 'xyz')
        beta (np.ndarray): Camera extrinsic calibration
            x (across shore), y (longshore), z (vertical), azimuth, tilt, roll
        lcp (dict): Lens Calibration Profile structure, the intrinsic camera calibration
        P (np.ndarray): Matrix containing intrinsic and extrinsic calibration
    """
    def __init__(self, calibration_file):
        calibration_file = Path(calibration_file)
        self.fname = calibration_file.name
        sn, cn, dc, cs, _ = self.fname.split('_')
        self.serial_number = int(sn)
        self.camera_number = cn
        self.calibration_date = datetime.datetime.strptime(dc, '%Y%m%d')
        self.coordinate_system = cs

        mat_data = scipy.io.loadmat(calibration_file)
        self.beta = mat_data['beta'][0]
        self.lcp = self._load_lcp(mat_data['lcp'])
        self.P = self._assembleP()

    def _load_lcp(self, lcp):
        """Return dict of lcp from lcp loaded from mat file"""
        NU = lcp[0, 0][0][0][0]
        NV = lcp[0, 0][1][0][0]
        c0U = lcp[0, 0][2][0][0]
        c0V = lcp[0, 0][3][0][0]
        fx = lcp[0, 0][4][0][0]
        fy = lcp[0, 0][5][0][0]
        d1 = lcp[0, 0][6][0][0]
        d2 = lcp[0, 0][7][0][0]
        d3 = lcp[0, 0][8][0][0]
        t1 = lcp[0, 0][9][0][0]
        t2 = lcp[0, 0][10][0][0]
        r = lcp[0, 0][11][0, :]
        caltech_fname = lcp[0, 0][12][0]
        fr = lcp[0, 0][13][0, :]
        x = lcp[0, 0][14][0, :]
        y = lcp[0, 0][15][0, :]
        dx = lcp[0, 0][16][:, :]
        dy = lcp[0, 0][17][:, :]

        return {
            'NU': NU,
            'NV': NV,
            'c0U': c0U,
            'c0V': c0V,
            'fx': fx,
            'fy': fy,
            'd1': d1,
            'd2': d2,
            'd3': d3,
            't1': t1,
            't2': t2,
            'r': r,
            'fr': fr,
            'caltech_fname': caltech_fname,
            'x': x,
            'y': y,
            'dx': dx,
            'dy': dy
        }

    def __repr__(self):
        msg = (
            f'serial_number: {self.serial_number}\n'
            f'camera_number: {self.camera_number}\n'
            f'calibration_date: {self.calibration_date}\n'
            f'coordinate_system: {self.coordinate_system}\n'
            f'beta: {self.beta}\n'
            f"sum of lcp r: {np.nansum(self.lcp['r'])}"
        )
        return msg

    def __str__(self):
        msg = (
            f'serial_number: {self.serial_number}, '
            f'camera_number: {self.camera_number}, '
            f'calibration_date: {self.calibration_date}, '
            f'coordinate_system: {self.coordinate_system}'
        )
        return msg

    def _assembleP(self):
        """Assembles and returns Projective (P) matrix from LCP and Beta values.

        Notes:
            - Derived from lcpBeta2P.m + CiRN notes
            - K converts angle away from the center of view into camera coordinates
            - R describes the 3D viewing direction of camera compared to world coordinates
            - beta[:3] camera location in world coordinates (x,y,z)
            - beta[3::] camera orientation (azimuth, tilt, roll)

        Returns:
            P (np.ndarray): Projective matrix
        """
        # K: intrinsic matrix, puts image in pixel units of the specific camera
        K = np.array([
            [self.lcp['fx'], 0,               self.lcp['c0U']],
            [0,              -self.lcp['fy'], self.lcp['c0V']],
            [0,              0,               1]
        ])
        # R: rotation matrix, puts image in camera orientation
        R = angle2R(
            self.beta[3],
            self.beta[4],
            self.beta[5]
        )
        # I: identify matrix augmented by camera center, puts image in camera coordinates
        IC = np.vstack((
            np.eye(3),
            -self.beta[:3]
        )).T
        KR = np.matmul(K, R)
        P = np.matmul(KR, IC)

        # Make the matrix homogenous, methods use homogenous coordinates for easier math
        # - normalize to make last element equal 1
        P = P/P[-1, -1]

        return P


def angle2R(azimuth, tilt, swing):
    """Assembles and returns a rotation matrix R from azimuth, tilt, and swing (roll)

    Notes:
        - derived from angles2R.m by Costal Imaging Research Network and Oregon State University
        - From p 612 of Wolf, 1983

    Arguments:
        azimuth (float): Azimuth
        tilt (float): Tilt
        swith (float): swing

    Returns:
        R (np.ndarray): Rotation matrix
    """
    a = azimuth
    t = tilt
    s = swing
    R = np.zeros((3, 3))

    R[0, 0] = np.cos(a) * np.cos(s) + np.sin(a) * np.cos(t) * np.sin(s)
    R[0, 1] = -np.cos(s) * np.sin(a) + np.sin(s) * np.cos(t) * np.cos(a)
    R[0, 2] = np.sin(s) * np.sin(t)

    R[1, 0] = -np.sin(s) * np.cos(a) + np.cos(s) * np.cos(t) * np.sin(a)
    R[1, 1] = np.sin(s) * np.sin(a) + np.cos(s) * np.cos(t) * np.cos(a)
    R[1, 2] = np.cos(s) * np.sin(t)

    R[2, 0] = np.sin(t) * np.sin(a)
    R[2, 1] = np.sin(t) * np.cos(a)
    R[2, 2] = -np.cos(t)

    return R
