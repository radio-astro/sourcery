##Current version: 0.2.9


##From version 0.2.8 to 0.2.9, the changes made were;

#### Shapes are corrected: sourcery able to assign pnt and Gau appropriately, and uses DC_Maj/Min as a model and Maj/Min for reliability computation.
#### The flux estimations are corrected. Uses masked image to create islands and the actual image for detection (model fitting)
#### Takes the radius of sources to exclude inside the reliability script (previously it was done inside main-- had problems)
#### Has an option to remove sources with correlation < 0.002 and reliability > 0.60
