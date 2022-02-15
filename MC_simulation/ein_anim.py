import numpy as np

from einsteinpy.geodesic import Timelike
from einsteinpy.plotting import StaticGeodesicPlotter


position = [4, np.pi / 3, 0.]
momentum = [0., 0.767851, 2.]
a = 0.99
steps = 400.
delta = 0.5

geod = Timelike(
    metric="Kerr",
    metric_params=(a,),
    position=position,
    momentum=momentum,
    steps=steps,
    delta=delta,
    return_cartesian=True
)

sgpl = StaticGeodesicPlotter()
sgpl.animate(geod, interval=1)
sgpl.show()

sgpl.ani.save('animation.gif', writer='imagemagick', fps=60)
