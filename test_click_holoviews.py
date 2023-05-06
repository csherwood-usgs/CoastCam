# import the required library
import numpy as np
import holoviews as hv


img_dir = 'D:/crs/proj/2019_CACO_CoastCam/CACO-02/2023-03-22_extrinsics/images/'
img_fname = '1679497200.Wed.Mar.22_15_00_00.GMT.2023.caco-02.c2.timex.jpg'



from PIL import Image
from holoviews import opts, streams
from holoviews.plotting.links import DataLink
hv.extension('bokeh')

# Load image with pixel values as bounds
img_pil = Image.open(img_dir+img_fname)
bounds = (0, 0, img_pil.width, img_pil.height)
img = hv.RGB(np.array(img_pil), bounds=bounds)

# Create the PointDraw stream and link it to a Table
points = hv.Points(([], []))
point_stream = streams.PointDraw(data=points.columns(), source=points)
table = hv.Table(points, ['x', 'y'])
DataLink(points, table)
hv.render(img * points + table).opts(
    opts.RGB(data_aspect=1),
    opts.Layout(merge_tools=False),
    opts.Points(active_tools=['point_draw'], size=10, tools=['hover'],
                height=400, width=600),
    opts.Table(editable=True))
