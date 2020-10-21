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
        metadata
        extrinsics
        intrinsics
        local_origin

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
    def __init__(self, metadata, intrinsics, extrinsics, local_origin):
        self.fname = metadata['name']
        # this assumes a file naming convention we have not adopted
        self.serial_number = metadata['serial_number']
        self.camera_number = metadata['camera_number']
        self.calibration_date = metadata['calibration_date']
        # coordinate system is either 'xyz' or 'geo'
        self.coordinate_system = metadata['coordinate_system'].lower()

        self.lcp = intrinsics
        if self.coordinate_system == 'geo':
            self.world_extrinsics = extrinsics
            self.local_extrinsics = {}
        elif self.coordinate_system == 'xyz':
            self.world_extrinsics = {}
            self.local_extrinsics = extrinsics
        else:
            print('Invalid value of coordinate_system: ',metadata['coordinate_system'])

        self.local_origin = local_origin

        if self.coordinate_system == 'geo':
            # calculate local extrinsics
            self.beta, self.local_extrinsics = self._convert_beta()
        else:
            # already in local coordinates, but need to make beta array
            self.beta = beta = np.array([*self.local_extrinsics.values()], dtype='float64')

        self.P, self.R, self.IC = self._assembleP()

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

    def _convert_beta(self):
        """Changes world coordinates to local coordinate_system"""
        local_xo = self.local_origin['x']
        local_extrinsics = local_transform_extrinsics(\
            self.local_origin['x'],self.local_origin['y'],self.local_origin['angd'],\
            1,\
            self.world_extrinsics)
        beta = np.array([*local_extrinsics.values()], dtype='float64')
        return beta, local_extrinsics

    def _assembleP(self):
        """Assembles and returns Projective (P) matrix from LCP and Beta values.

        Notes:
            - Derived from lcpBeta2P.m + CiRN notes
            - K converts angle away from the center of view into camera coordinates
            - R describes the 3D viewing direction of camera compared to world coordinates
            - beta[:3] camera location in world coordinates (x,y,z)
            - beta[3::] camera orientation (azimuth, tilt, roll) in radians

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

        return P, R, IC


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
