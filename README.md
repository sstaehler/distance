[![Build Status](https://travis-ci.org/sstaehler/distance.svg?branch=master)](https://travis-ci.org/sstaehler/distance)

# distance
Compute the distance of an earthquake given a seismic velocity model and a travel time difference between two phases (typically P and S).

## Installation
Using pip:
```shell script
$ git clone https://github.com/sstaehler/distance
$ cd distance
$ pip install -v -e .
```
## Example:
Contains a python script that can be called from shell:
```shell script
$ python -m distance.distance ./taup_file.nd 200 300 400
./taup_file.nd, S-P time: 200.0: NO SOLUTION FOUND!
./taup_file.nd, S-P time: 300.0, distance:  46.0
./taup_file.nd, S-P time: 400.0, distance:  62.9
```
as well as a Python function 
```python
from distance.distance import get_distance
from obspy.taup import TauPyModel
model = TauPyModel('iasp91')

dist = get_distance(tSmP=400, 
                    depth=50.)
print(dist)
```
