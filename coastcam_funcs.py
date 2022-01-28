#
import numpy as np
import datetime
import os
from dateutil import tz
from skimage import io
from skimage.color import rgb2gray

import json
import yaml
import matplotlib.pyplot as plt

# TODO - Pass image arrays, rather than path names to these functions
#     - ...and eliminate packages required to load images

def estimate_sharpness(img):
    """
    Estimate image sharpness and contrast
    Input:
        img - np.ndarray representing image
              img read as skimage.io.imread( 'imagefilename.jpg' )
    Returns:
        s,c - sharpness and contrast estimates

    https://stackoverflow.com/questions/6646371/detect-which-image-is-sharper
    """
    array = np.asarray(rgb2gray(img), dtype=np.int32)
    contrast = array.std()
    gy, gx = np.gradient(array)
    gnorm = np.sqrt(gx**2 + gy**2)
    sharpness = np.average(gnorm)
    return sharpness, contrast

def average_color(img):
    """ Calculate the average pixel intensity of an image
    Input:
        img - np.ndarray representing image
              img read as skimage.io.imread( 'imagefilename.jpg' )
    Returns:
        av, avall - av (np.array of average r, g, b values), avall average of r,g,b
    """
    av = img.mean(axis=0).mean(axis=0)
    avall = av.mean(axis=0)
    return av, avall

def detect_blur_fft(img, size=60, vis=False):
    """ Use high-frequency content of image fft to determine blur
    Input:
        img - np.ndarray representing image
              img read as skimage.io.imread( 'imagefilename.jpg' )
    Returns:
        av, avall - av (np.array of average r, g, b values), avall average of r,g,b
    From: https://www.pyimagesearch.com/2020/06/15/opencv-fast-fourier-transform-fft-for-blur-detection-in-images-and-video-streams/

    This is slower than estimate_sharpness and not any better as a metric for quality
    """

    # grab the dimensions of the image and use the dimensions to
    # derive the center (x, y)-coordinates
    imgbw = rgb2gray(img)
    (h, w) = imgbw.shape
    (cX, cY) = (int(w / 2.0), int(h / 2.0))
    fft = np.fft.fft2(img)
    fftShift = np.fft.fftshift(fft)
    # zero-out the center of the FFT shift (i.e., remove low
    # frequencies), apply the inverse shift such that the DC
    # component once again becomes the top-left, and then apply
    # the inverse FFT
    fftShift[cY - size:cY + size, cX - size:cX + size] = 0
    # check to see if we are visualizing our output
    if vis:
        # compute the magnitude spectrum of the transform
        np.seterr(all='ignore') # suppress log of zero error
        magnitude = 20 * np.log(np.abs(fftShift))
        # display the original input image
        (fig, ax) = plt.subplots(1, 2, )
        ax[0].imshow(img, cmap="gray")
        ax[0].set_title("Input")
        ax[0].set_xticks([])
        ax[0].set_yticks([])
        # display the magnitude image
        ax[1].imshow(magnitude, cmap="gray")
        ax[1].set_title("Magnitude Spectrum")
        ax[1].set_xticks([])
        ax[1].set_yticks([])
        # show our plots
        plt.show()

    fftShift = np.fft.ifftshift(fftShift)
    recon = np.fft.ifft2(fftShift)
    # compute the magnitude spectrum of the reconstructed image,
    # then compute the mean of the magnitude values
    magnitude = 20 * np.log(np.abs(recon))
    mean = np.mean(magnitude)
    # the image will be considered "blurry" if the mean value of the
    # magnitudes is less than the threshold value
    return mean

def json2dict(jsonfile):
    """ Import contents of a JSON file as a dict

    Args:
        jsonfile (str): json2dict file to read
    Returns:
        dict interpreted from JSON file
    """
    with open(jsonfile, "r") as data:
        dictname = json.loads(data.read())
    return dictname


def yaml2dict(yamlfile):
    """ Import contents of a YAML file as a dict

    Args:
        yamlfile (str): YAML file to read
    Returns:
        dict interpreted from YAML file
    """
    dictname = None
    with open(yamlfile, "r") as infile:
        try:
            dictname = yaml.safe_load(infile)
        except yaml.YAMLerror as exc:
            print(exc)
    return dictname


def dts2unix(date_time_string, timezone='eastern'):
    """
    Return the unix epoch time and a datetime object with aware UTC time zone

    Args:
        date_time_string in format YY-MM-DD hh:mm
        time zone for date_time string

    Returns:
        epoch number, datetime_object
    """
    if timezone.lower() == 'eastern':
        tzone = tz.gettz('America/New_York')
    elif timezone.lower() == 'utc':
        tzone = tz.gettz('UTC')

    date_time_obj = datetime.datetime.strptime(date_time_string, '%Y-%m-%d %H:%M').replace(tzinfo=tzone)
    ts = date_time_obj.timestamp()
    return int(ts), date_time_obj

def unix2dts(unixnumber, timezone='eastern'):
    """
    Get local time from unix number

    Input:
        unixnumber - string containing unix time (aka epoch)
    Returns:
        date_time_string, date_time_object in utc

    TODO: not sure why this returns the correct value without specifying that input time zone is eastern
    """
    if timezone.lower() == 'eastern':
        tzone = tz.gettz('America/New_York')
    elif timezone.lower() == 'utc':
        tzone = tz.gettz('UTC')

    # images other than "snaps" end in 1, 2,...but these are not part of the time stamp.
    # replace with zero
    ts = int( unixnumber[:-1]+'0')
    date_time_obj =  datetime.datetime.utcfromtimestamp(ts)
    date_time_str = date_time_obj.strftime('%Y-%m-%d %H:%M:%S')
    return date_time_str, date_time_obj

def filetime2timestr(filepath, timezone='eastern'):
    """
    Return the local time and the Unix Epoch string from an image filename or path
    Does not work with backslashes (e.g., Windows paths)
    """
    if timezone.lower() == 'eastern':
        tzone = tz.gettz('America/New_York')
    elif timezone.lower() == 'utc':
        tzone = tz.gettz('UTC')

    # remove path
    filename = os.path.split(os.path.normpath(filepath))[-1]
    # split on '.', take first on
    s = filename.split('.')[0]

    #TODO - Could check camera type and correct last digit, but it does not affect seconds

    date_time_str, date_time_obj = unix2dts(s)
    return date_time_str, s

def timestr2filename(date_time_str, camera = 'c1', image_type = 'timex', timezone='eastern'):
    """
    Return a filename given a date_time_str and other info
    """
    # filenames have extra digit added to time stamps - here is a dict listing them
    last_number = {'snap': 0, 'timex': 1, 'var': 2, 'bright': 3, 'dark': 4, 'rundark': 5}
    if timezone.lower() == 'eastern':
        tzone = tz.gettz('America/New_York')
    elif timezone.lower() == 'utc':
        tzone = tz.gettz('UTC')

    date_time_obj = datetime.datetime.strptime(date_time_str, '%Y-%m-%d %H:%M').replace(tzinfo=tzone)
    ts = int(date_time_obj.timestamp())+int(last_number[image_type])

    fn = str(ts)+'.'+camera+'.'+image_type+'.jpg'
    return fn
