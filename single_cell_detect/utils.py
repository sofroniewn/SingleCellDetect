###########################################################################
# Utilities for SingleCellDetect package
#
# Polar to Cartesian coordinate image transformation code lightly edited
# from https://github.com/NoahApthorpe/CellMagicWand
#
###########################################################################

from numpy import percentile, cos, sin, sqrt, arctan2, max, pad, meshgrid, \
    linspace, arange, round, clip, pi, zeros
from scipy.ndimage.interpolation import zoom

def norm(image, min_percentile=0, max_percentile=100):
    """
    Clip and rescale an image to be between [0, 1]

    Parameters
    ----------
    image : a 2d image
        The image to be rescaled.
    min : float
        Lower percentile value.
    max : float
        Upper percentile value.
    """

    tmp = image.astype('float').clip(percentile(image,min_percentile),percentile(image,max_percentile))
    tmp = tmp - tmp.min()
    return tmp/tmp.max()


def coord_polar_to_cart(r, theta, center):
    """
    Converts polar coordinates around center to Cartesian
    """
    x = r * cos(theta) + center[0]
    y = r * sin(theta) + center[1]
    return x, y


def coord_cart_to_polar(x, y, center):
    """
    Converts Cartesian coordinates to polar
    """
    r = sqrt((x-center[0])**2 + (y-center[1])**2)
    theta = arctan2((y-center[1]), (x-center[0]))
    return r, theta


def image_cart_to_polar(image, center, min_radius, max_radius, phase_width, zoom_factor=1):
    """
    Converts an image from cartesian to polar coordinates around center
    """

    # Upsample image
    if zoom_factor != 1:
        image = zoom(image, (zoom_factor, zoom_factor), order=4)
        center = (center[0]*zoom_factor + zoom_factor/2, center[1]*zoom_factor + zoom_factor/2)
        min_radius = min_radius * zoom_factor
        max_radius = max_radius * zoom_factor

    # pad if necessary
    max_x, max_y = image.shape[0], image.shape[1]
    pad_dist_x = max([(center[0] + max_radius) - max_x, -(center[0] - max_radius)])
    pad_dist_y = max([(center[1] + max_radius) - max_y, -(center[1] - max_radius)])
    pad_dist = int(max([0, pad_dist_x, pad_dist_y]))
    if pad_dist != 0:
        image = pad(image, pad_dist, 'constant')

    # coordinate conversion
    theta, r = meshgrid(linspace(0, 2*pi, phase_width),
                           arange(min_radius, max_radius))
    x, y = coord_polar_to_cart(r, theta, center)
    x, y = round(x), round(y)
    x, y = x.astype(int), y.astype(int)

    polar = image[x, y]
    polar.reshape((max_radius - min_radius, phase_width))

    return polar


def image_polar_to_cart(mask, center, min_radius, max_radius, output_shape, zoom_factor=1):
    """
    Converts a polar binary mask to Cartesian and places in an image of zeros
    """

    # Account for upsampling
    if zoom_factor != 1:
        center = (center[0]*zoom_factor + zoom_factor/2, center[1]*zoom_factor + zoom_factor/2)
        min_radius = min_radius * zoom_factor
        max_radius = max_radius * zoom_factor
        output_shape = map(lambda a: a * zoom_factor, output_shape)

    # new image
    image = zeros(output_shape)

    # coordinate conversion
    theta, r = meshgrid(linspace(0, 2*pi, mask.shape[1]),
                           arange(0, max_radius))
    x, y = coord_polar_to_cart(r, theta, center)
    x, y = round(x), round(y)
    x, y = x.astype(int), y.astype(int)

    x = clip(x, 0, image.shape[0]-1)
    y = clip(y, 0, image.shape[1]-1)
    ix,iy = meshgrid(arange(0,mask.shape[1]), arange(0,mask.shape[0]))
    image[x,y] = mask

    # downsample image
    if zoom_factor != 1:
        zf = 1/float(zoom_factor)
        image = zoom(image, (zf, zf), order=4)

    return image
