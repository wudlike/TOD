# usage of the code

## 1 Smoothing.py

### 1.1 run smoothing.py

The purpose of this module is to save the smoothed maps (TT, EE, BB) for
cmb, thermaldust and synchrotron.

There are three options when running this module,

* -td : if True, then the thermaldust maps will be smoothed and saved at ../data/psm/smoothed_comps/thermaldust/
      default is False
* -sync : if True, the maps of synchrotron will be smoothed and saved at ../data/psm/smoothed_comps/synchrotron/
      default is False
* -cmb  : if True, the maps of cmb will be smoothed and saved at ../data/psm/smoothed_comps/cmb/
      default is False
The usage of smoothing.py as follows:
getting help:

```python
python smoothing -h
```

running code only smoothing and saving cmb maps:

```python
python smoothing -cmb True
```

### 1.2 saved data format

The format of the saved data is as follows:
