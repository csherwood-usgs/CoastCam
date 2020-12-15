## CoastCam - Code for processing USGS CoastCam data

This resides under ../crs/src/CoastCam on my laptop and under ../src/CoastCam on my Workspaces instance.ls

### Code from Axiom provided by Brittany Bruder with minimum modifications to run with USACE examples

`rectifyer.py` - Axiom code for rectification  
`calibration.py` - Axiom code for reading .mat calib files and generating transformation matrices  
`test_rectificatyion.ipynb` - Code to test Axiom routines with USACE data  


### Code modified for USGS CoastCam data
There are a few differences between the USACE examples and the USGS data stored on CHS. USACE extrinsics are in local coordinates (e.g., alongshore/cross-shore; USGS extrinsics are in world coordinates (e.g., UTM northing/easting). USACE calibration data is stored in .mat files and includes some metadata. As of now, USGS calibration info is stored in two files (intrinsic and extrinsic) in JSON format that looks like the dict it is read into. The USGS routines need some additional metadata...right now, this is provided in the form of a dict that is hard-coded. (We are only processing one station).

`write_IOEO_to_json.m` - Matlab script to convert calibrations in `.mat` format to JSON files.  
`rectifier_crs.py` - Code to rectify images, USGS version.  
`calibration_crs.py` - Code to ingest intrinsic and extrinsic calibration data and provide rotation matrices.  
`coastcam_funcs.py` - Collection of useful functions.  
`rectify_caco01` - Working version for rectifying images on the S3 bucket  

#### Test notebooks to explore various processing

`test_read_json_cal_files.ipynb` - Code to read calibration files in proposed .json format.   
`test_bucket_read.ipynb` - Demos using fsspec to read/write files to S3 buckets.  
`test_R2Calc.ipynb` - Python version of Stockdon equation.  
`test_time_funcs.ipynb` - Demo of datetime routines and converting file names to times.  
`test_histogram_balance.ipynb` - Demo of using histogram balance to better match images.   

