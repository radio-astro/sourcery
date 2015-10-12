# Testing

On a command line:
> sourcery -i kat7restored.fits -p kat7psf.fits -pref "Test" -nisl=1 -npix=1

### Output Files 

1. A directory named 'Test_datenow'.  
   Inside this diretory is Two PNG files and Two Source catalogues in Tigger format (name.lsm.html)
2. A log file named 'logfile.log' which contains any relevant information which might be useful to a user.
