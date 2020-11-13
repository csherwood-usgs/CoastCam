import imageio
import numpy as np
from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator
from scipy.ndimage.morphology import distance_transform_edt

from calibration_crs import CameraCalibration #CRS


class TargetGrid(object):
    """Grid generated to georectify image.

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
    def __init__(self, xlims, ylims, dx=1, dy=1, z=-0.91):
        x = np.arange(xlims[0], xlims[1]+dx, dx)
        y = np.arange(ylims[0], ylims[1]+dx, dy)
        self.X, self.Y = np.meshgrid(x, y)
        self.Z = np.zeros_like(self.X) + z
        self.xyz = self._xyz_grid()

    def _xyz_grid(self):
        x = self.X.copy().T.flatten()
        y = self.Y.copy().T.flatten()
        z = self.Z.copy().T.flatten()
        return np.vstack((x, y, z)).T


class Rectifier(object):
    """Georectifies an oblique image given RectifierGrid and ncolors.

    Note:
        Derived from example code for georectifying mini-Argus imagery.

    Attributes:
        target_grid (TargetGrid): Params and grid used to create georectified image.
        camera_calibration (CameraCalibration): CameraCalibration including intrinsic (LCP) and extrinsic (Beta) coefficients.
        ncolors (int): Number of colors in camera images.
        U (np.ndarray): horizontal image coordinates
        V (np.ndarray): vertical image coordinates (increasing down)
    """
    def __init__(self, target_grid, ncolors=3):
        self.target_grid = target_grid
        self.ncolors = ncolors

    def _find_distort_UV(self, calibration):
        # get UV for pinhole camera
        xyz = np.vstack((
            self.target_grid.xyz.T,
            np.ones((len(self.target_grid.xyz),))
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

        DU = Ud.reshape(self.target_grid.X.shape, order='F')
        DV = Vd.reshape(self.target_grid.Y.shape, order='F')

        # find negative Zc values and add to flag
        UV = np.matmul(calibration.P, xyz)
        xyzC = np.matmul(calibration.R,np.matmul(calibration.IC,xyz))
        flag[np.where(xyzC[2,:]<=0.)]=0.
        flag = flag.reshape(self.target_grid.X.shape, order='F')
        
        # apply the flag to zero-out non-valid points
        return DU*flag, DV*flag, flag

    def get_pixels(self, DU, DV, image, interp_method='rgi'):
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
        K = np.zeros((
            self.target_grid.X.shape[0],
            self.target_grid.X.shape[1],
            self.ncolors
        ))

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

    def assemble_image_weights(self, K):
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

    def apply_weights_to_pixels(self, K, W):
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

    def rectify_images(self, metadata, image_files, intrinsic_cal_list, extrinsic_cal_list, local_origin, fs=None, interp_method = 'rgi'):
        """Georectify and blend images from multiple cameras 

        Arguments:
            metadata (dict):
            image_files (list): List of image files
            intrinsic_cal_list (list): list of paths to internal calibrations (one for each camera)
            extrinsic_cal_list (list): list of paths to external calibrations (one for each camera)
            local_origin:
            fs: (object): fsspec file spec object for folder on S3 bucket. If none, normal file system will be used.
            interp_method: (string): either 'rbs' (rectilinear bicubic spline) or 'rgi' (regular grid interpolator: linear and faster)
            camera_calibration_files (list): List of calibrations for cameras used to get image_files.

        Returns:
            M (np.ndarray): Georectified images merged from supplied images.
        """
        # array for final pixel values
        M = np.tile(
            np.zeros_like(self.target_grid.X[:, :, np.newaxis]),
            (self.ncolors,)
        )
        # array for weights
        totalW = M.copy()

        for cur_idx, (image_file, intrinsic_cal, extrinsic_cal) in enumerate(zip(image_files, intrinsic_cal_list, extrinsic_cal_list)):
            #  print("loop",cur_idx,"calibrations:")
            #  print(intrinsic_cal, extrinsic_cal)
            # load camera calibration file and find pixel locations
            camera_calibration = CameraCalibration(metadata, intrinsic_cal, extrinsic_cal, local_origin)
            U, V, flag = self._find_distort_UV(camera_calibration)

            # load image and apply weights to pixels
            if fs:
                # using fsspec for S3 files
                with fs.open(image_file) as f:
                    image = imageio.imread(f)
            else:
                # regular file system
                image = imageio.imread(image_file)
                
            K = self.get_pixels(U, V, image, interp_method=interp_method)
            W = self.assemble_image_weights(K)
            K_weighted = self.apply_weights_to_pixels(K, W)

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
