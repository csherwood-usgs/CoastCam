from PIL import Image
import os
import glob

def shrink_all_jpgs(fac=2):
    '''
    Shrink all jpgs in the current folder by 1/fac

    Input:
        fac - factor to reduce image by (int)

    Returns:
        None - but writes a new, smaller image file named oldname.x2.jpg
    '''

    # get a list of all .jpg files in the current dir
    flist = map(os.path.basename, glob.glob('*.jpg'))

    for f in flist:
        # split name at dots
        n = f.split('.')
        # change last part
        n[-1]='x2.jpg'
        # reconstruct name
        tmp = list(map(str,n))
        new_name = '.'.join(tmp)

        with Image.open(f) as im:
            (width, height) = (im.width // int(fac), im.height // int(fac))
            # default uses bicubic interpolation when resizing
            im_resized = im.resize((width, height))
            im_resized.save(new_name,quality=95)
            print(new_name)

if __name__ == '__main__':
    shrink_all_jpgs()
