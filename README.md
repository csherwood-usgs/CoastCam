## CoastCam - Code for processing USGS CoastCam data

This resides under ../crs/src/CoastCam on my laptop and under ../src/CoastCam on my Workspaces instance.ls

Routines are in two files:
`coastcam.py` - Main calibration and rectification routines
`coastcam_funcs.py` - Various helper and utility routines

### coastcam.py

### coastcam_funcs.py - Collection of utility routines

#### Image quality
`estimate_sharpness(img)` - Estimate image sharpness and contrast

`average_color(img)` - Calculate the average pixel intensity of an image

`detect_blur_fft(img, size=60, vis=False)` - Use high-frequency content of image fft to determine blur

#### JSON
`def json2dict(jsonfile)` - Read a .json file into a dict

#### Date/time/filename
`def dts2unix(date_time_string, timezone='eastern')` - Return the unix epoch time and a datetime object with aware UTC time zone

`unix2dts(unixnumber, timezone='eastern')` - Get local time from unix number

`filetime2timestr(filepath, timezone='eastern')` - Return the local time and the Unix Epoch string from an image filename or path

`timestr2filename(date_time_str, camera = 'c1', image_type = 'timex', timezone='eastern')` - Return a filename given a date_time_str and other info
 
### `./awss` - Shell scripts with useful AWS CLI commands for examining S3

### `./data` - Image files and calibration data for tests

### `./misc` - Misc. code that probably should not be in this repo

### `./original_code` - Code from Axiom provided by Brittany Bruder with minimum modifications to run with USACE examples

`rectifyer.py` - Axiom code for rectification  
`calibration.py` - Axiom code for reading .mat calib files and generating transformation matrices  
`test_rectificatyion.ipynb` - Code to test Axiom routines with USACE data 

### `./tests` - Examples and test code

### `./util` - Utilities

`write_IOEO_to_json.m` - Matlab script to convert calibrations in `.mat` format to JSON files.
`parse_s3_inventory.ipynb`

#### Test notebooks to explore various processing

`test_read_json_cal_files.ipynb` - Code to read calibration files in proposed .json format.   
`test_bucket_read.ipynb` - Demos using fsspec to read/write files to S3 buckets.  
`test_R2Calc.ipynb` - Python version of Stockdon equation.  
`test_time_funcs.ipynb` - Demo of datetime routines and converting file names to times.  
`test_histogram_balance.ipynb` - Demo of using histogram balance to better match images.   

