# sourcery
Tools for creating high fidelity source catalogues from radio interferometric datasets.
It also provides tools to select sources in these images that are responsible for artifacts.

## Install

### Pip
```
pip install sourcery
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
pip uninstall sourcery
```

### Direct build
```
cd sourcery
cat sourcery_files.txt | sudo xargs rm
```
