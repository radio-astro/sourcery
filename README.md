# sourcery
Tools for creating high fidelity source catalogues from radio interferometric datasets.
It also provides a tool for selecting sources responsible for artifacts in an image.

## Requires

* pybdsm (part lofar software package)
* tigger
* numpy
* scipy
* astlib

You can find ubuntu packages for tigger, astlib, lofar on [radio-astro ppa](https://launchpad.net/~radio-astro/+archive/ubuntu/main). You may need to build the them source for other distributions, in that case find the source at the [ska-sa github repo](https://github.com/ska-sa)


## Install

### Pip
```
sudo pip install sourcery
```

### Direct build
```
git clone github.com/SpheMakh/sourcery
cd sourcery
python setup.py install --record sourcery_files.txt
```

## Unistall
### Pip
```
sudo pip uninstall sourcery
```

### Direct build
```
cd sourcery
cat sourcery_files.txt | sudo xargs rm
```


### Running command 
Using defaults:
```
    sourcery -i fitsimage -p psfimage  
```
For changing parameters:
    ```sourcery -h
    ```

If using config.json:  
   
1. Provide a -pref on the command line, e.g   
```
-prefix KAT7_DIRECTORY 

```
  
 2. You can either provide image input name inside the config file or on the comand line.
   In case you provide the input image inside the cofig.json only a single image can be executed per config file.
   But if providing the input image on the command line then more than one image can be specified and exceuted.
   One can use both the config file and the command line such as    
   
```
sourcery -i kat7.fits,kat8.fits -p kat7psf.fits -jc config.json

```

In which the parameter settings used are specified in config.json and the image to be executed are kat7.fits, kat8.fits etc.

If using the command line alone then sourcery -h should be your friend. 
Note more than one image can be specfied:    

```

sourcery -i kat7.fits,kat8.fits -p kat7psf.fits

```
The psf "kat7psf.fits" in this case will be used for both images, if each image has its own psf then provide this excplicitly as  
```
     sourcery -i kat7.fits,kat8.fits -p kat7psf.fits,kat8psf.fits
     
```
