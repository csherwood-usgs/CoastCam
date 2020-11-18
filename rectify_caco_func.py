from pathlib import Path
import imageio
import fsspec
import numpy as np
import matplotlib.pyplot as plt
import datetime
from dateutil import tz
import pandas as pd

from coastcam_funcs import *
from calibration_crs import *
from rectifier_crs import *

rectify_caco_func(impaths,fs):
    """
    Rectify one pair of images for CACO01
    """
    # List of files...three for each camera. Calibration parameters are in .json format
    # These are the USGS image filename format
    extrinsic_cal_files = ['CACO01_C1_EOBest.json','CACO01_C2_EOBest.json']
    intrinsic_cal_files = ['CACO01_C1_IOBest.json','CACO01_C2_IOBest.json']

    # Dict providing the metadata that the Axiom code infers from the USACE filename format
    metadata= {'name': 'CACO-01', 'serial_number': 1, 'camera_number': 'C1', 'calibration_date': '2019-12-12', 'coordinate_system': 'geo'}
    # dict providing origin and orientation of the local grid
    local_origin = {'x': 410935.,'y':4655890., 'angd': 55.}

    # read cal files and make lists of cal dicts
    extrinsics_list = []
    for f in extrinsic_cal_files:
        extrinsics_list.append( json2dict(f) )
    intrinsics_list = []
    for f in intrinsic_cal_files:
        intrinsics_list.append( json2dict(f) )

    calibration = CameraCalibration(metadata,intrinsics_list[0],extrinsics_list[0],local_origin)

    xmin = 0.
    xmax = 500.
    ymin = 0.
    ymax = 700.
    dx = 1.
    dy = 1.
    z =  0.

    rectifier_grid = TargetGrid(
        [xmin, xmax],
        [ymin, ymax],
        dx,
        dy,
        z
    )

    rectifier = Rectifier(
        rectifier_grid
    )

    rectified_image = rectifier.rectify_images(metadata, impaths, intrinsics_list, \
                                    extrinsics_list, local_origin, fs=fs)
    return rectified_image