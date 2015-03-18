from scipy.io.netcdf import netcdf_file
import scipy.interpolate, numpy, math

# Constants:
plate_scaling_factor    = 1.0 # Factor to scale plate velocities to RMS velocity of model (RMS_Model / RMS_Earth)
velocity_non_dim_factor = 1.0 # Factor to non-dimensionalise gplates velocities: often = d/kappa, but set here to 1.0
time_dim_factor         = 1.0 # Factor to dimensionalise model time: often = d^2/kappa, but set here to 1.0
deg2rad                 = math.pi/180.
gplates_stage_time      = 31536000000000.0 * (1./plate_scaling_factor)    # 1 Myr in seconds (scaled)
total_stages            = 200
start_time              = 6275663000000000.0

# Coordinate transformation function:
def xyz2spherical(X):
  # Returns r, phi, theta for a x,y,z coordinate. Note that phi is a weird latitude that starts at 0 
  # at the south pole and is pi at the north pole.
  r = math.sqrt(X[0]**2+X[1]**2+X[2]**2)
  phi = math.atan2(X[1], X[0]) + math.pi
  theta = math.pi-math.acos(X[2]/r)
  return r, phi, theta

class Gplates_Interpolator(object):

  def __init__(self):
    self.current_time = None
    self.stage_time   = None
    self.stage        = None
    self.lats         = None
    self.lons         = None
    self.coords       = None
    self.set_time(start_time)

  def set_time(self, t):
    if (self.current_time == t * time_dim_factor):
      return
    self.current_time = t * time_dim_factor
    self.stage_time   = self.current_time % gplates_stage_time
    print 'GPlates: Stage time (s):',self.stage_time,'of:',gplates_stage_time

    # Opens the GPlates data file appropriate for the current simulation time (stage):
    total_stage_time = gplates_stage_time * (1. / plate_scaling_factor)
    new_stage        = int(math.ceil(total_stages - (self.current_time / gplates_stage_time) ))

    if (self.stage == new_stage):
      return
    self.stage = new_stage

    # If appropriate, read GPlates file:
    if(self.stage != 200):
      print "GPlates *** Plate Stage Has Changed *** "
    print 'GPlates: Current dimensional time (s):',self.current_time
    print "GPlates *** Reading from file: velocity_%s.00Ma.nc *** " % str(self.stage) 

    nc = netcdf_file('velocity_%s.00Ma.nc' % str(self.stage), 'r')

    # As lats and lons do not change between files, only store these once:
    if self.lats is None:
      self.lats = (nc.variables['lat'][1:-1] + 90.)*deg2rad # do not include poles
    if self.lons is None:
      self.lons = (nc.variables['lon'][:-1] + 180.)*deg2rad # do not include 0 meridian twice
    if self.coords is None:
      self.coords = numpy.array([(lat,lon) for lat in self.lats for lon in self.lons])

    # Set up interpolators. Note that a nearest neighbour interpolation is used for Plate ID's to maintain unique values:
    self.intpx  = scipy.interpolate.RectSphereBivariateSpline(self.lats, self.lons, nc.variables['velocity'][0,:-1,1:-1].T)
    self.intpy  = scipy.interpolate.RectSphereBivariateSpline(self.lats, self.lons, nc.variables['velocity'][1,:-1,1:-1].T)
    self.intpz  = scipy.interpolate.RectSphereBivariateSpline(self.lats, self.lons, nc.variables['velocity'][2,:-1,1:-1].T)
    self.intpid = scipy.interpolate.NearestNDInterpolator(self.coords, nc.variables['plate_id'][:-1,1:-1].T.reshape(len(self.lats)*len(self.lons)))

  # Plate ID:
  def get_plate_ID(self,X):
    r, phi, theta = xyz2spherical(X)
    return self.intpid(theta,phi)

  # Nodal velocity vector:
  def get_velocities(self,X):
    r, phi, theta = xyz2spherical(X)
    return (self.intpx(theta,phi)[0,0] * velocity_non_dim_factor * plate_scaling_factor,
            self.intpy(theta,phi)[0,0] * velocity_non_dim_factor * plate_scaling_factor, 
            self.intpz(theta,phi)[0,0] * velocity_non_dim_factor * plate_scaling_factor)

  # X velocity:
  def get_velocity_x(self,X):
    r, phi, theta = xyz2spherical(X)
    return (self.intpx(theta,phi)[0,0] * velocity_non_dim_factor * plate_scaling_factor)

  # Y velocity:
  def get_velocity_y(self,X):
    r, phi, theta = xyz2spherical(X)
    return (self.intpy(theta,phi)[0,0] * velocity_non_dim_factor * plate_scaling_factor)

  # Z velocity:
  def get_velocity_z(self,X):
    r, phi, theta = xyz2spherical(X)
    return (self.intpz(theta,phi)[0,0] * velocity_non_dim_factor * plate_scaling_factor)

# Set up interpolation class:
gplates_interpolator=Gplates_Interpolator()





