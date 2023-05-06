# %%
# import the required library
import cv2

# define a function to display the coordinates of
# of the points clicked on the image

img_dir = 'D:/crs/proj/2019_CACO_CoastCam/CACO-02/2023-03-22_extrinsics/images/'
img_fname = '1679497200.Wed.Mar.22_15_00_00.GMT.2023.caco-02.c2.timex.jpg'

def click_event(event, x, y, flags, params):
   if event == cv2.EVENT_LBUTTONDOWN:
      print(f'({x},{y})')
      
      # put coordinates as text on the image
      cv2.putText(img, f'({x},{y})',(x,y),
      cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 0, 255), 2)
      
      # draw point on the image
      cv2.circle(img, (x,y), 3, (0,255,255), -1)
 
# read the input image
img = cv2.imread(img_dir+img_fname)

# create a window
cv2.namedWindow('Point Coordinates')

# bind the callback function to window
cv2.setMouseCallback('Point Coordinates', click_event)

# display the image
while True:
   cv2.imshow('Point Coordinates',img)
   k = cv2.waitKey(1) & 0xFF
   if k == 27:
      break
cv2.destroyAllWindows()


