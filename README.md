# Grid Edge Classification
The purpose of this program is to create a list of grid edges that are classified to represent the effects of levees (or walls) in a raster grid flood inundation model. 

The method that is implemented in this code is described in a research
paper by Kahl et al. (2022) appearing in Advances in Water Resources
https://doi.org/10.1016/j.advwatres.2022.104287

The code was developed to improve the representation of levees in
dual-grid flood inundation models - which resolve topographic data on a
fine grid and perform solution updates on a relatively coarse grid. In
particular, the authors developed the code in support of the PRIMo model
described by Sanders and Schubert (2019)
 https://doi.org/10.1016/j.advwatres.2019.02.007

However, we anticipate that the method (and code) could support any raster
grid flood inundation model. That is, the code could be applied to come up
with a list of grid edges that, when connected, can approximately
represent a levee that would otherwise fall between neighboring grid
edges. It is left to the user to reformat the output for supporting other
applications.

This code relies on commands from the following matlab toolboxes
       Image Processing Toolbox
       Mapping Toolbox

Within the code, several additional inputs must be specified by the user:

    1) the location and input filenames for the DEM and levee shapefile.
    
    2) the upscale factor, which represents the ratio of the coarse grid
       size to the fine grid size. 
       
    3) the row and column indices of grid cell "plugs" that are needed to
       fully encircle interior channel areas with cells. That is, before
       running the "imfill.m" command, it is important to make sure there
       are not any openings/gaps in loop of cells around each channel.
       
    4) the row and column indices of cells that fall within "holes"
       See Fig. 2b in Kahl et al. (2022). This information supports the
       "imfill.m" routine (hole filling operation)
       
    5) the row and column indices of grid edges that are incorrectly
       classified as levees and, if implemented, would create a
       non-physical blockage effect within a channel. We term these edges
       as "weirs" since they would present a barrier similar to a weir
       if not removed. See Fig. 2c in Kahl et al. (2022).
       Note that there are separate lists for horizontal and vertical
       "weir" edges that neee to be removed.

The output of this program is an ascii file (*.levee) containing two 
lists of edges:

    1) A list of vertical edges defined by the row and cell index
    
    2) A list of horizontal edges defined by the row and cell index

The method is configured for a testcase at upscale factor of 10 and user inputs are declared. 

