# SingleCellDetect (WIP)

> Detect a single cell located in the center of a 2d image

This package provided several methods for detecting the boundaries of a single cell located at the center of an image and is useful when performing cell finding on biological images or movie. It may be used together with other packages that finds the locations of cell centers in a larger image.

### install

```bash
git clone https://github.com/sofroniewn/SingleCellDetect
git clone https://github.com/sofroniewn/CellMagicWand
```
Add both SingleCellDectect and CellMagicWand to your `PYTHONPATH`.

### example
```python
from single_cell_detect import watershed_edge as detect
from skimage.data import coins

image = coins()[:100,:100]
mask = detect(image, dilationSize=1, radial=True, filterSize=5)
```

### usage

#### `mask = watershed_edge(image, dilationSize=0, radial=True, filterSize=0)`

Returns a single cell mask found using a marker controlled watershed on the
edges of `image`. The `dilationSize`  controls the amount of dilation applied
to the central maker for the marker controlled segmentation. If `radial` is
true the edge detection will be done radially by first transforming image to
polar coordinates. The `filterSize` controls the amount of median filtering
applied after edge detection before the watershed. If radial is false the
filtering is with a disk, if radial is true the filtering is tangent to the
radial direction.


### license
MIT
