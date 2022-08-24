## CoastCam - Code for processing USGS CoastCam data

This resides under ../crs/src/CoastCam on my laptop and under ../src/CoastCam on my Workspaces instance.  

Main routines are in two files:  
* `coastcam.py` - Main calibration and rectification routines  
* `coastcam_funcs.py` - Various helper and utility routines  

### coastcam.py
#### Main routines
* `class TargetGrid(target_grid_info,z=0.0)` - Build grid onto which image(s) will be rectified
* `class CameraCalibration(cfg, intrinsics, extrinsics)` - Generate camera calibration object w/ transformation matrices
* `rectify_images(cfg, target_grid, image_files, intrinsic_cal_list, extrinsic_cal_list, fs=None, interp_method = 'rgi')` - Georectify and blend images from multiple cameras

#### Functions called by main routines
* `_convert_beta(self)` - Changes world coordinates to local coordinate_system
* `assembleP(lcp, beta)` - Assembles and returns Projective (P) matrix from LCP and Beta values
* `apply_weights_to_pixels(K, W)` - Return pixel intensities (K) weighted by W.
* `local_transform_points( xo, yo, ang, flag, xin, yin)` - Transforms points between local World Coordinates and Geographical (either direction)
* `local_transform_extrinsics(local_xo,local_yo,local_angd,flag,extrinsics_in)` - Tranforms extrinsic vector between Local World Coordinates and Geographical
* `angle2R(azimuth, tilt, swing)` - Assembles and returns a rotation matrix R from azimuth, tilt, and swing (roll)
* `assemble_image_weights(K)` - Return weight matrix W used for image merging.
* `get_pixels(nx, ny, nc, DU, DV, image, interp_method='rgi')` - Return pixel values for each xyz point from the image
* `find_distort_UV(target_grid, calibration)` - get UV for pinhole camera

### coastcam_funcs.py - Collection of utility routines
##### Image quality
* `estimate_sharpness(img)` - Estimate image sharpness and contrast
* `average_color(img)` - Calculate the average pixel intensity of an image
* `detect_blur_fft(img, size=60, vis=False)` - Use high-frequency content of image fft to determine blur
##### JSON
`def json2dict(jsonfile)` - Read a .json file into a dict
##### Date/time/filename
* `dts2unix(date_time_string, timezone='eastern')` - Return the unix epoch time and a datetime object with aware UTC time zone
* `unix2dts(unixnumber, timezone='eastern')` - Get local time from unix number
* `filetime2timestr(filepath, timezone='eastern')` - Return the local time and the Unix Epoch string from an image filename or path
* `timestr2filename(date_time_str, camera = 'c1', image_type = 'timex', timezone='eastern')` - Return a filename given a date_time_str and other info
### Subfolders   
#### `./awss` - Shell scripts with useful AWS CLI commands for examining S3
#### `./data` - Image files and calibration data for tests
#### `./misc` - Misc. code that probably should not be in this repo
#### `./original_code` - Code from Axiom provided by Brittany Bruder with minimum modifications to run with USACE examples
* `rectifyer.py` - Axiom code for rectification  
* `calibration.py` - Axiom code for reading .mat calib files and generating transformation matrices  
* `test_rectificatyion.ipynb` - Code to test Axiom routines with USACE data 
#### `./tests` - Examples and test code
#### `./util` - Utilities
* `write_IOEO_to_json.m` - Matlab script to convert calibrations in `.mat` format to JSON files.
* `parse_s3_inventory.ipynb` 

### Test notebooks to explore various processing
`test_read_json_cal_files.ipynb` - Code to read calibration files in proposed .json format.   
`test_bucket_read.ipynb` - Demos using fsspec to read/write files to S3 buckets.  
`test_R2Calc.ipynb` - Python version of Stockdon equation.  
`test_time_funcs.ipynb` - Demo of datetime routines and converting file names to times.  
`test_histogram_balance.ipynb` - Demo of using histogram balance to better match images.   

