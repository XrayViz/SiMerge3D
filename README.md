# SimpleMerge3D


A tool to make 3D volume using series of diffraction patterns. 

The patterns should correspond to a rotation series (an angle for each pattern has to be specified and the rotational axis for the whole measurement) or just a set of patterns with known Euler angles (for example simulation made with Moltrans).

The geometry of the detector is specified in CrystFEL format, but if you don't have one - the program would show you an example.

The program can read different file types (so far implemented .tif, .cbf, .raw, ...?) and saves the resulting 3D volume as a binary raw file or as a 3D tiff (ImageJ can read those).

Command line arguments:
...

**Useful tips**
- HDF5 support is currently excluded, but I plan to add it as an option during compilation. I just hate that library, but, unfortunately, there is no way around it :(
- _Scale_ parameter can be used to "down-sample" the resulting 3D array or, is the angular step was fine enough, even "up-sample". Usually dont-sampling is really needed - it the reciprocal space resolution in the experiment was higher than needed. In this case due to better statistics even weak signal, not visible in 2D, can be observed in 3D. 
- RAM required for the merging is easy to estimate: Dims3D_X*Dims3D_Y*Dims3D_Z * 8bytes (the float array + int array)
- _Interpolation_ option is useful when the angular step was too big during the measurements. But one has to be careful - it might take long and also it consumes a lot of memory. Due to the fact that it has to store im memory every "good" pixel from all patterns (I, x, y, z, pointer), it requires: num_good_pixels * num_patterns * 24 bytes. So having 800x800 pixels ROI and 200 patterns would result in almost 3Gb memory demand. And this is added to the size of the resulting array
- It is very important to know the geometry of your experiment. Especially the detector center and the sample-detector distance. **Wavelenght doesn't matter!** If you need help with optimization of the detector's geometry - contact me, I can advise several solutions that I've been using at different facilities around the world
- The program currently assumes a flat detector perpendicular to the beam, but it is easy to modify for more complex situation. At least to accomodate the detector at some angle one can replace the simplified geometry reading functions by the function from CrystFEL project - then the angle of each panel can be taken into account. 
