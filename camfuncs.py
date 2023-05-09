import scipy.io
import numpy as np

"""
readUVd

readxyz

read_intrinsics

read_extrinsics

xyz2UVd

UVd2xyz

distortUV

undistortUV

make_matrices(intrinsics, extrinsics)

angle2R(extrinsics)

"""

def make_matrices(intrinsics, extrinsics):
    """Assembles and returns Projective (P) matrix from intrinsics and extrinsics.

    Notes:
        - Derived from lcpBeta2P.m + CiRN notes
        - K converts angle away from the center of view into camera coordinates
        - R describes the 3D viewing direction of camera compared to world coordinates
        - beta[:3] camera location in world coordinates (x,y,z)
        - beta[3::] camera orientation (azimuth, tilt, roll) in radians

    Input:
        intrinsics == lcp
        extrinsics == beta

    Returns:
        P (np.ndarray): Projective matrix
    """
    # K: intrinsic matrix, puts image in pixel units of the specific camera
    K = np.array([
        [intrinsics['fx'],        0,           intrinsics['c0U']],
        [0,                -intrinsics['fy'],  intrinsics['c0V']],
        [0,                       0,                 1          ]
    ])
    # R: rotation matrix, puts image in camera orientation
    R = angle2R(
        beta[3],
        beta[4],
        beta[5]
    )
    # I: identify matrix augmented by camera center, puts image in camera coordinates
    IC = np.vstack((
        np.eye(3),
        -beta[:3]
    )).T
    KR = np.matmul(K, R)
    P = np.matmul(KR, IC)

    # Make the matrix homogenous, methods use homogenous coordinates for easier math
    # - normalize to make last element equal 1
    P = P/P[-1, -1]

    return P, R, IC

def angle2R( extrinsics ):
    """Assembles and returns a rotation matrix R from azimuth, tilt, and swing (roll)

    Notes:
        - most recently extracted from CRS CoastCam
        - derived from angles2R.m by Costal Imaging Research Network and Oregon State University
        - From p 612 of Wolf, 1983

    Arguments:
        azimuth (float): azimuth (degre
        tilt (float): tilt
        swing (float): swing

    Returns:
        R (np.ndarray): Rotation matrix
    """
    dtr = np.pi/180.
    a = extrinsics['azimuth']*dtr
    t = extrinsics['tilt']*dtr
    s = extrinsics['swing']*dtr
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


intrinsics_path = 
extrinsics_path = 
gcpCoords = 

image_path = 

extrinsics_guess = 
extrsinsics_know_flag = [ 1 1 1 0 0 0];  # [ x y z azimuth tilt swing]

#  Azimuth is the horizontal direction the camera is pointing and positive CW 
#  from World Z Axis. 

#  Tilt is the up/down tilt of the camera. 0 is the camera looking nadir,
#  +90 is the camera looking at the horizon right side up. 180 is looking
#  up at the sky and so on.

#  Swing is the side to side tilt of the camera.  0 degrees is a horizontal 
#  flat camera. Looking from behind the camera, CCW rotation of the camera
#  would provide a positve swing.

#  Diagrams of these defintions are in Section 6 of the user manual. 

# load gcpxyz real-world coordinatesB

# load gcpUVd pixel coordinates

# assert they are the same length

xyz = np.array( (nid, 3) )
UVd = np.array( (nid, 2) )

for i, id in enumerate( gcpUVd.keys) :
   print(id)
   xyz[i, 0] = gcpxyz[id][1] #x, easting
   xyz[i, 1] = gcpxyz[id][0] #y, northing
   xyz[i, 2] = gcpxyz[id][1] #z, elevation
   UVd[i, 0] = gcpUVd[id][0] #pixel U
   UVd[i, 1] = gcpUVd[id][1] #pixel V

