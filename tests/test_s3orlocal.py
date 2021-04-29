def s3orlocal(path):
    # determine if path resolves to s3 or s3 local
    isS3 = False
    return isS3

if __name__ == "__main__":
    # test the functions here
    path = 'junk'
    print('Testing isS3',s3orlocal(path))
else:
    print(__name__,' is loaded as a module.')
