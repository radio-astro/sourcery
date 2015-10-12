# Testing

###On a command line:
> sourcery -i kat7restored.fits -p kat7psf.fits -pref "Test" -nisl=1 -npix=1

### Output Files 

1. A directory named 'Test_datenow' where datenow specfifies the date and hour at which you run the test.  
   Inside this diretory is 2 PNG files and 2 Source catalogues in Tigger format (name.lsm.html) as well as 
Fits image  Test_negative.fits.
2. A log file named 'logfile.log' which contains any relevant information which might be useful to a user.
