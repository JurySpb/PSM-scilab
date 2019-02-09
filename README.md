# PSM-scilab
### Point Shift Measurement (PSM) - Another way to find the shift of background screen points
The temperature of the gaseous medium is simply related to the refractive index through the well-known Lorenz-Lorentz equation or through Gladstone – Dale. Recently, a technique based on the use of a background screen - Background Oriented Schlieren [(BOS)](https://en.wikipedia.org/wiki/Background-oriented_schlieren_technique) has become widespread. The BOS method uses software particle image velocimetry [(PIV)](https://en.wikipedia.org/wiki/Particle_image_velocimetry). PIV calculates the shift caused by a change in the refractive index in the object under study. In some cases, the PIV gives erroneous data. The PSM program calculates the shift in a different way. The program works in the SciLab environment and uses the IPCV package. The program finds points on the image, determines their boundaries, then calculates the distance between the boundaries. 

Requires to work:
1. [Scilab](https://www.scilab.org/)
2. Image Processing and Computer Vision Toolbox [(IPCV)](https://atoms.scilab.org/toolboxes/IPCV/)
