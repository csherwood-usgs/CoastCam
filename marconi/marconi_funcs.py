import numpy as np
import os
import matplotlib.pyplot as plt
# from scipy import interpolate
from scipy.stats import linregress
from pyproj import Proj, transform
import yaml


def pcoord(x, y):
    """
    Convert x, y to polar coordinates r, az (geographic convention)
    r,az = pcoord(x, y)
    """
    r = np.sqrt(x**2 + y**2)
    az = np.degrees(np.arctan2(x, y))
    # az[where(az<0.)[0]] += 360.
    az = (az+360.)%360.
    return r, az


def xycoord(r, az):
    """
    Convert r, az [degrees, geographic convention] to rectangular coordinates
    x,y = xycoord(r, az)
    """
    x = r * np.sin(np.radians(az))
    y = r * np.cos(np.radians(az))
    return x, y


def UTM2local(eutm, nutm, eoff=420080.0, noff=4638320.0, rot=20.0):
    """
    Convert UTM NAD83 Zone 19N easting, northing to local Marconi alongshore, cross-shore coordinates
    xloc, yloc = UTM2local( eutm, nutm )
    Better to use values from the dict than defaults for translation/rotation values
    
    """
    [r, az] = pcoord(eutm-eoff, nutm-noff)
    az = az + rot
    [xloc, yloc] = xycoord(r,az)
    return xloc, yloc


def local2UTM(alongshore, across_shore, eoff=420080.0, noff=4638320.0, rot=20.):
    """Convert local coordinates to UTM
       Inverse of UTM2local()
       Better to use values from the dict than defaults for translation/rotation values

       Here is code for UTM2local:
          [r, az] = pcoord(eutm-eoff, nutm-noff)
          az = az + rot
          [xisl,yisl] = xycoord(r,az)
    """
    r, az = pcoord(alongshore, across_shore)
    az = az - rot
    eUTM, nUTM = xycoord(r, az)
    eUTM = eUTM + eoff
    nUTM = nUTM + noff
    return eUTM, nUTM


def LatLon2UTM(lat,lon,init_epsg='epsg:26918'):
    """
    Convert lat lon (WGS84) to UTM.
    Defaults to Zone 18N

    TODO: Update to Proj 6 and correct this syntax
    """
    inProj = Proj(init='epsg:4326')
    outProj = Proj(init=init_epsg)
    outx,outy = transform(inProj,outProj,lon,lat)
    return outx, outy


def UTM2LatLon(easting, northing, initepsg='epsg:26918'):
    """
    Convert UTM to lat, lon (WGS84)
    Defaults to Zone 18N

    TODO: Update to Proj 6 and correct this syntax
    """
    outProj = Proj(init='epsg:4326')
    inProj = Proj(init=initepsg)
    lon,lat = transform(inProj,outProj,easting,northing)
    return lon, lat


def rotrans(x, y, x0, y0, theta):
    '''
    2D rotation and translation of x, y
    Input:
        x, y - row vectors of original coordinates (must be same size)
        x0, y0 - Offset (location of x, y = (0,0) in new coordinate system)
        theta - Angle of rotation (degrees, CCW from x-axis == Cartesian coorinates)
    Returns:
        x_r, y_r - rotated, offset coordinates
    '''
    thetar = np.radians(theta)
    c, s = np.cos(thetar), np.sin(thetar)

    # homogenous rotation matrix
    Rh = np.array(((c, -s,  0.),
                   (s,  c,  0.),
                   (0., 0., 1.)))
    # homogenous translation matrix
    Th = np.array(((1., 0., x0),
                   (0., 1., y0),
                   (0., 0., 1.)))

    # homogenous input x,y
    xyh = np.vstack((x,y,np.ones_like(x)))

    # perform rotation and translation
    xyrh = np.matmul(np.matmul(Th,Rh), xyh)
    x_r = xyrh[0,:]
    y_r = xyrh[1,:]
    return x_r, y_r


def make_grid(name='marconi_local', e0=420080.0, n0=4638320.0,
    xs=0., xend=200., ys=-100., yend=300., dxdy=1., theta=10.):
    """
    Make a rectangular grid to interpolate elevations onto.

    where:
      e0 - UTM Easting of origin [m]
      n0 - UTM Northing of origin [m]
      xs - Start of alongshore axis [m]
      xend - End of alongshore axis [m]
      ys - Start of cross-shore axis [m]
      yend - End of cross-shore axis [m]
      dxdy - grid size (must be isotropic right now) [m]
      theta - rotation CCW from x-axis [deg]
    """
    xlen = (xend - xs)/dxdy
    ylen = (yend - ys)/dxdy
    nx = int((1./dxdy) * xlen)
    ny = int((1./dxdy) * ylen)
    print(nx, ny)

    xcoords = np.linspace(xs+0.5*dxdy,xend-0.5*dxdy,nx)
    ycoords = np.linspace(ys+0.5*dxdy,yend-0.5*dxdy,ny)

    # these will be the coordinates in rotated space
    xrot, yrot = np.meshgrid(xcoords, ycoords, sparse=False, indexing='xy')

    print('make_grid: Shape of xrot, yrot: ',np.shape(xrot),np.shape(yrot))
    shp = np.shape(xrot)
    xu, yu = rotrans(xrot.flatten(), yrot.flatten(), e0, n0, theta)
    xu=np.reshape(xu,shp)
    yu=np.reshape(yu,shp)
    # write the UTM coords of the corners to an ASCII file
    corners = np.asarray(  [[xu[0][0],yu[0][0]],
                           [xu[0][-1],yu[0][-1]],
                           [xu[-1][-1],yu[-1][-1]],
                           [xu[-1][0],yu[-1][0]],
                           [xu[0][0],yu[0][0]]])

    print('corners x, corners y (orig. coords)')
    print(corners)
    fn_corners = name+'.csv'
    print('Saving to '+fn_corners)
    np.savetxt(fn_corners, corners, delimiter=",")
    return nx, ny, xu, yu, xrot, yrot, xcoords, ycoords
