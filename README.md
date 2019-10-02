# flap_mdsplus
FLAP MDSPlus interface

This module provides interface between MDSplus data reading and the flap.DataObject. It read data and its coordinates from the given MDSplus tree and node combination.
A sample configuration file is provided for NSTX. Furthermore, a wrapper object is provided for EFIT data.

Comment for installing MDSPlus:
- MDSplus package should be installed in the OS from mdsplus.org / Software or from the main repo of the OS
http://www.mdsplus.org/index.php?title=Downloads&open=6913025428015808543&page=Software%2FDownloads
- Then the python package should be installed by python shell then:
http://www.mdsplus.org/index.php/Documentation:Users:MDSobjects:Python

The following install was not working on macOS:
conda install -c conda-forge mdsplus
