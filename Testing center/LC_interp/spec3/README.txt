
This tar file contains the files necessary to fit and/or plot the XRT spectra
yourself. The files all begin spec3 - which is the internal name given to this
spectum. The files spec3pc.pi and/or spec3wt.pi are those which were
fitted by the automatic processing. They have the BACKFILE, ANCRFILE and RESPFILE
keywords set, "bad 0-29" applied (as use of XRT data below 0.3 keV is not
recommended) and have had "group min 1" applied by grppha so the data can be
fitted with the C-statistic.

The files *back.pi, *.rmf and *.arf files are the background and ancilliary response
files respectively, which are needed for fitting the data. The *source.pi files
are the source spectrum before having any keywords set, or any binning applied.

If you use the spec3[wt|pc].pi files, Xspec will look for the back.pi and
RMF and ARF in your current working directory. 

If you make use of these data in any publication, please cite Evans et al. 2009
(MNRAS, 397, 1177) and include the following text in your 
acknowledgements:

This work made use of data supplied by the UK Swift Science Data Centre at the
University of Leicester.

