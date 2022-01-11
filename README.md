# Delsarte-K3-Picard-Numbers
## What is this?
In this repository I computed the Picard numbers of all smooth Delsarte K3 surfaces.
You can re-run the computations or just access the results.
This is an implementation of the algorithm explained in the article "An Explicit Algorithm for Computing the Picard Number of Certain Algebraic Surfaces" by Tetsuji Shioda.
The implementation computes the Picard number for each matrix encoding a smooth Delsarte K3 surface.
Many matrices encode the same K3 surface, for example the matrices:

```
A=[[4,0,0,0],[0,4,0,0],[0,0,4,0],[0,0,0,4]],
B=[[0,4,0,0],[4,0,0,0],[0,0,4,0],[0,0,0,4]]
```

The program runs separate computations even if matrices encode the surface.
Removing this could seriously speed up the computation. But it is already fairly fast, so I did not bother.

## Requirements

Python 3, Numpy 1.22

## Usage

Use the function `get_picard` in `delsarte_picard_computation.py` to compute the Picard number of a single smooth Delsarte K3 surface. For example:

```
import numpy as np
from delsarte_picard_computation import get_picard
A = np.array([[4,0,0,0],[0,4,0,0],[0,0,4,0],[0,0,0,4]])
print(get_picard(A))
```

Gives the output `20`, which is the Picard number of the K3 surface encoded by the matrix `A`.
This K3 surface is the Fermat quartic.

Run `delsarte_picard_computation.py` to compute the Picard numbers of all smooth Delsarte K3 surfaces.
The result will be written to a log file.

Check the file `__main__.log` for the results of such a run.
It contains 5761 K3 surfaces (some of them isomorphic, see above) with their Picard numbers.
For example, the line 

```2022-01-10 18:01:10,965 - shioda_picard_computation.py:182 - INFO - ID: 22536, rho=4, A=[array([0, 0, 0, 4]), array([1, 0, 3, 0]), array([0, 3, 1, 0]), array([3, 0, 0, 1])]```

Means that the Delsarte surface encoded by the matrix

```[[0, 0, 0, 4], [1, 0, 3, 0], [0, 3, 1, 0], [3, 0, 0, 1]]```

has Picard number 4, this is what the equation `rho=4` means.
The lowest occuring Picard number is 4, the highest occuring Picard number is 20.