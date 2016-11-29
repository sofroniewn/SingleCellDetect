from numpy import zeros, pi
from skimage.morphology import watershed, disk, rectangle, dilation
from skimage.filters import sobel, sobel_h
from skimage.filters.rank import median
from .utils import norm, image_cart_to_polar, image_polar_to_cart

def watershed_edge(image, dilationSize=0, radial=True, filterSize=0):
    """
    Returns a single cell mask found using a marker controlled watershed on the
    edges of the image.

    Parameters
    ----------
    image : a 2d image
        The image in which to find the cell.
    dilationSize : int
        Amount of dilation applied to the central maker for the marker controlled
        segmentation.
    radial : bool
        If true the edge detection will be done radially by first transforming
        image to polar coordinates
    filterSize : int
        Amount of median filtering applied after edge detection before the watershed.
        If radial is false the filtering is with a disk, if radial is true the filtering
        is tangent to the radial direction.
    """

    if radial:
        edges = -sobel_r(image, filterSize)
    else:
        edges = sobel(image)
        edges = median(norm(edges), disk(filterSize))

    markers = zeros(image.shape)
    markers[markers.shape[0]/2,markers.shape[1]/2] = 1
    markers = dilation(markers, disk(dilationSize))

    markers[0, :] = 2
    markers[-1, :] = 2
    markers[:, 0] = 2
    markers[:, -1] = 2

    labels = watershed(edges, markers)

    return labels == 1

def sobel_r(image, filterSize):
    """
    Apply sobel filter radially to an image by first transforming it to polar coordinates
    then detect horizontal edges and then transform back to cartesian coordinates.

    Parameters
    ----------
    image : a 2d image
        The image for edge detection.
    filterSize : int
        Amount of median filtering applied after edge detection before the watershed.
        Filtering is done after the polar transform and is tangent to the radial direction.
    """

    phase_width = int(2 * pi * image.shape[1]/2 * 10)
    polar_image = image_cart_to_polar(image, [image.shape[0]/2, image.shape[1]/2], 0, image.shape[1]/2, phase_width=phase_width, zoom_factor=2)
    polar_image = median(norm(polar_image), rectangle(1,2*filterSize+1))
    edge_sobel = sobel_h(polar_image)
    return image_polar_to_cart(edge_sobel, [image.shape[0]/2, image.shape[1]/2], 0, image.shape[1]/2, image.shape, zoom_factor=2)
