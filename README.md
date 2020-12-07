# Semi-Autonomous-Method-for-Extracting-Nebular-Velocities-from-Images
Calculate the velocity of movement between corresponding regions (parsec scale distances from Earth) in two images (.fits) using RA and DEC coordinates, and generate a variety of graphics to represent said movement and velocity. [Python 2.7, SAOImage ds9 4.1]

This code will take a pair of pre-made SAOImage ds9 region files with the regions in each corresponding to points of interest that show movement between the images, calculate the velocity of said movement between each corresponding region, and generate four additional region files by mimicking the look of an actual SAOImage ds9 .reg file in .txt format, and changing the extention to .reg after saving the .txt file. 
The four additional region files comprise of:

 - Line regions between the pre-selected corresponding points of interest, labeled with numbers.
  
 - Line regions between the pre-selected corresponding ponts of interest, lebeled with the velocity of movement.
  
 - Vector regions between the pre-selected corresponding point of interest, labeled with numbers.
  
 - Scaled vector regions between the pre-selected corresponding points of interest, labeled with numbers.  
  
The generated vector region files require the degree angles for each vector to be entered manually. The degree angles can be found in the region descriptions of either generated line region file while open in SAOImage ds9. They must be entered in corresponding order. The numbering system is matched accross each generated region file with number labels.

Required packages: NumPy, os
