import numpy as np

from einsteinpy.geodesic import Timelike
from einsteinpy.plotting.geodesic import GeodesicPlotter

position = [40., np.pi / 2, 0.]
momentum = [0., 0., 3.83405]
a = 0.
steps = 5500
delta = 1.

geod = Timelike(
    metric="Schwarzschild",
    metric_params=(a,),
    position=position,
    momentum=momentum,
    steps=steps,
    delta=delta,
    return_cartesian=True
)

gpl = GeodesicPlotter()
gpl.plot(geod)
gpl.show()
