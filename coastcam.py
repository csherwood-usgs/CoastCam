import imageio
import datetime
from pathlib import Path
import scipy.io
import numpy as np
from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator
from scipy.ndimage.morphology import distance_transform_edt

class TargetGrid(object):
    """Build grid onto which image(s) will be rectified

    CRS modified to make endpoints inclusive
    Notes:
        - Used to maps points in world coordinates to pixels
        - The limits should be specified in local coordinates using the same
          coordinates and units as camera calibrations.

    Args:
        xlims (ndarray) - min and max (inclusive) in the x-direction (e.g. [-50, 650])
        ylims (ndarray) - min and max (inclusive) in the y-direction (e.g. [0, 2501])
        dx (float) - resolution of grid in x direction (same units as camera calibration)
        dy (float) - resolution of grid in y direction (same units as camera calibration)
        z (float) - static value to estimate elevation at everypoint in the x, y grid

    Attributes:
        X (np.ndarray): Local grid coordinates in x-direction.
        Y (np.ndarray): Local grid coordinates in y-direction.
        Z (np.ndarray): Local grid coordinates in z-direction.
        xyz (np.ndarray): The grid where pixels are compiled from images for rectification.
    """
    def __init__(self,tgi,z=0.0):
        # tgi is shorthand for target_grid_info dict
        # x = np.arange(xlims[0], xlims[1]+dx, dx)
        # y = np.arange(ylims[0], ylims[1]+dx, dy)
        x = np.arange(tgi['xmin'], tgi['xmax']+tgi['dx'],tgi['dx'])
        y = np.arange(tgi['ymin'], tgi['ymax']+tgi['dy'],tgi['dy'])
        self.X, self.Y = np.meshgrid(x, y)
        self.Z = np.zeros_like(self.X) + z
        self.xyz = self._xyz_grid()

    def _xyz_grid(self):
        x = self.X.copy().T.flatten()
        y = self.Y.copy().T.flatten()
        z = self.Z.copy().T.flatten()
        return np.vstack((x, y, z)).T

class CameraCalibration(object):
    """Generate camera calibration object w/ transformation matrices

    Notes:
        - Inspired by example code + notes from CiRC which are derived from Hartley and Zisserman (20030.)
        - Asssumes calibration saved in .mat file

    Args:
        configuration dict (abbreviated as cfg)
        extrinsics
        intrinsics
        local_origin

    Attributes:
        coordinate_system (str): Coordinate system used for calibration (e.g. 'xyz')
        beta (np.ndarray): Camera extrinsic calibration
            x (across shore), y (longshore), z (vertical), azimuth, tilt, roll
        lcp (dict): Lens Calibration Profile structure, the intrinsic camera calibration
        P (np.ndarray): Matrix containing intrinsic and extrinsic calibration
    """
    def __init__(self, cfg, intrinsics, extrinsics):
        # coordinate system is either 'xyz' or 'geo'
        self.coordinate_system = cfg['coordinate_system'].lower()
        print('in CameraCalibration',intrinsics)
        self.lcp = intrinsics
        if self.coordinate_system == 'geo':
            self.world_extrinsics = extrinsics
            self.local_extrinsics = {}
        elif self.coordinate_system == 'xyz':
            self.world_extrinsics = {}
            self.local_extrinsics = extrinsics
        else:
            print('Invalid value of coordinate_system: ',cfg['coordinate_system'])

        self.local_origin = cfg['local_origin']

        if self.coordinate_system == 'geo':
            # calculate local extrinsics
            self.beta, self.local_extrinsics = self._convert_beta()
        else:
            # already in local coordinates, but need to make beta array
            self.beta = beta = np.array([*self.local_extrinsics.values()], dtype='float64')

        self.P, self.R, self.IC = assembleP( self.lcp, self.beta )

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
        """Changes world coordinates to local coordinate_system
        Returns beta (np.array) and local_extrinsics (dict), which contain the same info
        """
        local_xo = self.local_origin['x']
        local_extrinsics = local_transform_extrinsics(\
            self.local_origin['x'],self.local_origin['y'],self.local_origin['angd'],\
            1,\
            self.world_extrinsics)
        beta = np.array([*local_extrinsics.values()], dtype='float64')
        return beta, local_extrinsics

def assembleP(lcp, beta):
    """Assembles and returns Projective (P) matrix from LCP and Beta values.

    Notes:
        - Derived from lcpBeta2P.m + CiRN notes
        - K converts angle away from the center of view into camera coordinates
        - R describes the 3D viewing direction of camera compared to world coordinates
        - beta[:3] camera location in world coordinates (x,y,z)
        - beta[3::] camera orientation (azimuth, tilt, roll) in radians

    Input:
        intrinsics == lcp
        extrinsics.values == beta

    Returns:
        P (np.ndarray): Projective matrix
    """
    # K: intrinsic matrix, puts image in pixel units of the specific camera
    K = np.array([
        [lcp['fx'], 0,               lcp['c0U']],
        [0,              -lcp['fy'], lcp['c0V']],
        [0,              0,               1]
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

def apply_weights_to_pixels(K, W):
    """Return pixel intensities (K) weighted by W.

    Arguments:
        K (np.ndarray): Pixel intensity
        W (np.ndarray): Pixel weights used for merging images

    Returns:
        K_weighted(np.ndarray): Pixel intensity weighted for merging
    """
    W_nonan = W.copy()
    K_weighted = K*W_nonan[:, :, np.newaxis]
    return K_weighted

def local_transform_points( xo, yo, ang, flag, xin, yin):
    """
    Transforms between local World Coordinates and Geographical
    World Coordinates. Local refers to the rotated coordinate system where x
    is positive offshore and y is oriented alongshore. This function can go
    from Local to Geographical or in reverse.

    Based on localTranformPoints.m by Brittany Bruder, but in radians

    Input:
        xo and yo - location of local origin (0,0) in Geographical coordinates.
              Typically xo is E and yo is N coordinate.
        ang - relative angle between the local X axis and the Geo X axis,
              positive counter-clockwise from the Geo X.  Units are radians.

              Note: Regardless of transformation direction, xo, yo, and ang
              should stay the same.

        flag = 1 or 0 to indicate transform direction
              Geo-->local (1) or
              local-->Geo (0)

        xin - Local (X) or Geo (E) coord depending on transformation direction
        yin = Local (Y) or Geo (N) coord depending on transformation direction

    Returns:
        xout - Local (X) or Geo (E) coord depending on transformation direction
        yout - Local (Y) or Geo (N) coord depending on transformation direction
    """

    if flag == 1:
        # transform from world -> local
        # translate from origin
        easp = xin-xo
        norp = yin-yo

        #rotate
        xout = easp*np.cos(ang)+norp*np.sin(ang)
        yout = norp*np.cos(ang)-easp*np.sin(ang)

    if flag == 0:
        # rotate
        xout = xin*np.cos(ang)-yin*np.sin(ang)
        yout = yin*np.cos(ang)+xin*np.sin(ang)
        # translate
        xout = xout+xo
        yout = yout+yo

    return xout, yout

def local_transform_extrinsics(local_xo,local_yo,local_angd,flag,extrinsics_in):
    """
    Tranforms between Local World Coordinates and Geographical
    World Coordinates for the extrinsics vector. Local refers to the rotated
    coordinate system where X is positive offshore and y is oriented
    alongshore. The function can go from Local to Geographical and in
    reverse. Note, this only performs horizontal rotations/transformations.

    Based on localTranformExtrinsics.m by Brittany Bruder

    Input:

    local_xo and local_yo - local origin = Location of Local (0,0) in
        Geographical Coordinates. Typically first entry is E and second is N coordinate.

    local_angd -The relative angle between the new (local) X axis and old (Geo)
        X axis, positive counter-clockwise from the old (Geo) X.  Units are degrees.

    Note: Regardless of transformation direction, local_ang, local_xo, and local_yo
        should stay the same. z, tilt, and roll do not change.

    extrinsics_in - dict with Local or Geo x, y, z, a(zimuth), t(ilt), r(oll)
        By CIRN convention, extrinsic angles are in radians.

    flag = 1 or 0 to indicate whether you are going from
        Geo-->Local (1) or
        Local-->Geo (0)

    Output:

    extrinsics_out - dict with Local or Geo x, y, z, a(zimuth), t(ilt), r(oll)
    """
    local_angr = np.deg2rad(local_angd)
    extrinsics_out = extrinsics_in.copy()
    if flag == 1:
        # World to local
        extrinsics_out['x'], extrinsics_out['y'] = local_transform_points(local_xo,local_yo,local_angr,1,extrinsics_in['x'],extrinsics_in['y'])
        extrinsics_out['a'] = extrinsics_in['a']+local_angr

    if flag == 0:
        # local to world
        extrinsics_out['x'], extrinsics_out['y'] = local_transform_points(local_xo,local_yo,local_ang,0,extrinsics_in['x'],extrinsics_in['y'])
        extrinsics_out['a'] = extrinsics_in['a']-local_angr

    return extrinsics_out

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

def assemble_image_weights(K):
    """Return weight matrix W used for image merging.

    Notes:
        - Calculates Euclidean distance for each entry to nearest non-zero pixel value
        - edt: Euclidean Distance Transform

    Arguments:
        K (np.ndarray): Pixel intensity

    Returns:
        W (np.ndarray): Pixel weights used for merging images
    """
    # NaN in K indicates no pixel value at that location
    # edt finds euclidean distance from no value to closest value
    # so, find the nans, then invert so it works with the function
    W = distance_transform_edt(~np.isnan(K[:, :, 0]))

    # Not sure when this would happen, but included because it's in the MATLAB code
    if np.isinf(np.max(W)):
        W[:] = 1
    W = W / np.max(W)
    return W

def get_pixels(nx, ny, nc, DU, DV, image, interp_method='rgi'):
    """Return pixel values for each xyz point from the image

    Arguments:
        DU (np.ndarray): Pixel location in camera orientation and coordinate system
        DV (np.ndarray): Pixel location in cmaera orientation and coorindate system
        image (np.ndarray [nx,ny,nc]) with RGB values at U,V points
        interp_method (string):
            'rgi' - uses SciPy RegularGridInterpolator (linear, about 5x faster)
            'rbs' - use SciPy RectBivariateSpline (smoother?)

    Returns:
        K (np.ndarray): Pixel intensity for each point in the image
    """

    K = np.zeros((nx,ny,nc))

    # Having tested both interpolation routines, the rgi is about five times
    # faster, no visual difference, but that has not been checked quantitatively.
    if interp_method == 'rbs':
        for c, _ in enumerate(['r', 'b', 'g']):
            rbs = RectBivariateSpline(
                # use this range to match matlab exactly
                np.arange(1, image.shape[0] + 1),
                np.arange(1, image.shape[1] + 1),
                image[:, :, c],
                kx=1,
                ky=1
            )
            K[:, :, c] = rbs.ev(DV, DU)
    elif interp_method == 'rgi':
        for c, _ in enumerate(['r', 'b', 'g']):
            rgi = RegularGridInterpolator(
                (np.arange(0, image.shape[0]),
                 np.arange(0, image.shape[1])),
                image[:,:,c],
                method='linear',
                bounds_error=False,
                fill_value=np.nan)
            K[:, :, c] = rgi((DV,DU))
    else:
        #TODO - what is the proper way to handle this error?
        print('No valid interp method')

    # mask out values out of range like Matlab
    # avoid runtime nan comparison warning (DU, DV already have nans)
    with np.errstate(invalid='ignore'):
        mask_u = np.logical_or(
            DU <= 1,
            DU >= image.shape[1]
        )
        mask_v = np.logical_or(
            DV <= 1,
            DV >= image.shape[0]
        )
    mask = np.logical_or(
        mask_u,
        mask_v
    )
    K[mask,:] = np.nan
    return K

def find_distort_UV(target_grid, calibration):
    # get UV for pinhole camera
    xyz = np.vstack((
        target_grid.xyz.T,
        np.ones((len(target_grid.xyz),))
    ))
    UV = np.matmul(calibration.P, xyz)

    # make homogenous
    div = np.tile(UV[2, :], (3, 1))
    UV = UV / div

    # get and rename
    NU = calibration.lcp['NU']
    NV = calibration.lcp['NV']
    c0U = calibration.lcp['c0U']
    c0V = calibration.lcp['c0V']
    fx = calibration.lcp['fx']
    fy = calibration.lcp['fy']
    d1 = calibration.lcp['d1']
    d2 = calibration.lcp['d2']
    d3 = calibration.lcp['d3']
    t1 = calibration.lcp['t1']
    t2 = calibration.lcp['t2']
    u = UV[0, :]
    v = UV[1, :]

    # normalize distances
    x = (u - c0U) / fx
    y = (v - c0V) / fy
    # radial distortion
    r2 = x*x + y*y
    fr = 1. + d1*r2 + d2*r2*r2 + d3*r2*r2*r2
    # tangential distorion
    dx=2.*t1*x*y + t2*(r2+2.*x*x)
    dy=t1*(r2+2.*y*y) + 2.*t2*x*y
    # apply correction, answer in chip pixel units
    xd = x*fr + dx
    yd = y*fr + dy
    Ud = xd*fx+c0U
    Vd = yd*fy+c0V

    # Declare array for flagged values
    flag = np.ones_like(Ud)

    # find negative UV coordinates
    flag[np.where( Ud<0.)]=0.
    flag[np.where( Vd<0.)]=0.
    # find UVd coordinates greater than image size
    flag[np.where( Ud>=NU)]=0.
    flag[np.where( Vd>=NV)]=0.

    # Determine if Tangential Distortion is within Range
    #  Find Maximum possible tangential distortion at corners
    Um=np.array((0, 0, NU, NU))
    Vm=np.array((0, NV, NV, 0))

    # Normalization
    xm = (Um-c0U)/fx
    ym = (Vm-c0V)/fy
    r2m = xm*xm + ym*ym

    # Tangential Distortion
    dxm=2.*t1*xm*ym + t2*(r2m+2.*xm*xm)
    dym=t1*(r2m+2.*ym*ym) + 2.*t2*xm*ym

    # Find Values Larger than those at corners
    flag[np.where(np.abs(dy)>np.max(np.abs(dym)))]=0.
    flag[np.where(np.abs(dx)>np.max(np.abs(dxm)))]=0.

    DU = Ud.reshape(target_grid.X.shape, order='F')
    DV = Vd.reshape(target_grid.Y.shape, order='F')

    # find negative Zc values and add to flag
    UV = np.matmul(calibration.P, xyz)
    xyzC = np.matmul(calibration.R,np.matmul(calibration.IC,xyz))
    flag[np.where(xyzC[2,:]<=0.)]=0.
    flag = flag.reshape(target_grid.X.shape, order='F')

    # apply the flag to zero-out non-valid points
    return DU*flag, DV*flag, flag

def rectify_images(cfg, target_grid, image_files, intrinsic_cal_list, extrinsic_cal_list, fs, interp_method = 'rgi'):
    """Georectify and blend images from multiple cameras

    Arguments:
        target_grid (class):
        cfg (dict):
        image_files (list): list of paths to image files (one for each camera)
        intrinsic_cal_list (list): list of paths to internal calibrations (one for each camera)
        extrinsic_cal_list (list): list of paths to external calibrations (one for each camera)
        local_origin:
        fs: (object): fsspec file spec object for folder on S3 bucket. If none, normal file system will be used.
        interp_method: (string): either 'rbs' (rectilinear bicubic spline) or 'rgi' (regular grid interpolator: linear and faster)
        camera_calibration_files (list): List of calibrations for cameras used to get image_files.

    Returns:
        M (np.ndarray): Georectified images merged from supplied images.
    """
    local_origin = cfg['local_origin']
    nc = cfg['n_colors']
    # array for final pixel values
    M = np.tile(
        np.zeros_like(target_grid.X[:, :, np.newaxis]),
        (nc,)
    )
    # array for weights
    totalW = M.copy()

    nx = target_grid.X.shape[0]
    ny = target_grid.X.shape[1]

    for cur_idx, (image_file, intrinsic_cal, extrinsic_cal) in enumerate(zip(image_files, intrinsic_cal_list, extrinsic_cal_list)):
        #  print("loop",cur_idx,"calibrations:")
        #  print(intrinsic_cal, extrinsic_cal)
        # load camera calibration file and find pixel locations
        camera_calibration = CameraCalibration(cfg, intrinsic_cal, extrinsic_cal)
        U, V, flag = find_distort_UV(target_grid, camera_calibration)

        # load image and apply weights to pixels
        if fs:
            # using fsspec for S3 files
            with fs.open(image_file) as f:
                image = imageio.imread(f)
        else:
            # regular file system
            image = imageio.imread(image_file)

        K = get_pixels(nx, ny, nc, U, V, image, interp_method=interp_method)
        W = assemble_image_weights(K)
        K_weighted = apply_weights_to_pixels(K, W)

        # add up weights and pixel itensities
        totalW = totalW + W[:, :, np.newaxis]
        K_weighted[np.isnan(K_weighted)] = 0
        M = M + K_weighted

    # stop divide by 0 warnings
    with np.errstate(invalid='ignore'):
        M = M / totalW

    #TODO - is there any need to retain the NaNs, or is replacing by zero ok?
    M[np.isnan(M)]=0

    #TODO - don't need to return W, K or flag...they are from last image processed
    # return M.astype(np.uint8), W, K, flag

    return M.astype(np.uint8)
