import astropy.units as u

from einsteinpy.rays import Shadow
from einsteinpy.plotting import ShadowPlotter
mass = 1 * u.kg
fov = 30 * u.km 
# What field of view is the user expecting
shadow = Shadow(mass=mass, fov=fov, n_rays=1000)
obj = ShadowPlotter(shadow=shadow, is_line_plot=True)
obj.plot()
obj.show()
obj = ShadowPlotter(shadow=shadow, is_line_plot=False)
obj.plot()
obj.show()
