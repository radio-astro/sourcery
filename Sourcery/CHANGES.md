##Current version: 0.3.1


##From version 0.3.0 to 0.3.1, the changes made were;

#### Masking is done as follows; first smoothing is performed at different scales and each of these are masked 
#### at a given threshold. The final mask is obtained by adding these masks and setting pixels > 1 to 1, this is followed
#### by multiplying the resulting mask and the original image.  
#### Whenever, no sources are not found by the source finder sourcery will give a warning and exit or continues with other images.
