import numpy as np
import astropy.units as u

def generate_inspiral_times(count, evolution_times, age=10*u.Gyr):
    """
        Create a random sample of insprial times *assuming constant star formation*
        
        Args:
            count           --> [int]             number of times to randomly sample
            evolution_times --> [array_like, Myr] t_evolve for each binary
            
            NOTE: the length of evolution times must be equal to count
            
        Returns:
            inspiral_times  --> [array_like, Gyr] time for which the binary will inspiral
            
            NOTE: many binaries will merge during this time (e.g. t_coalscence < t_inspiral)
    """
    
    if count == 1 and not isinstance(evolution_times.value, float) or len(evolution_times) != count:
        raise ValueError("Count must match the number of evolution times provided.")
    
    # randomise birth times uniformly across age
    # subtract maximum evolution time to ensure inspiral times are never negative
    formation_times = np.random.rand(count) * (age - max(evolution_times))

    # double compact object is formed t_form after t_birth
    DCO_times = formation_times + evolution_times

    # binary will inspiral for the rest of the age of the Milky Way
    inspiral_times = age - DCO_times
    
    return inspiral_times.to(u.Gyr)

def random_disk(count, scale_radius, scale_height):
    """ 
        Create a random sample of radii or heights using McMillan (2011)
        and the inverse CDF method.
        
        Args:
            count        --> [int]             number of lengths
            scale_radius --> [array_like, kpc] scale radius of the disk (Rd in McMillan)
            scale_height --> [array_like, kpc] scale height of the disk (Zd in McMillan)
        
        Returns:
            radii        --> [array_like, kpc] random radii in a disk
            heights      --> [array_like, kpc] random heights in a disk
            angles       --> [array_like, rad] random angles in a disk
    """
    # create a list of random radii
    u = np.random.rand(count)
    r = - scale_radius * np.log(1 - u)
    
    # create a list of random absolute heights
    u = np.random.rand(count)
    zabs = - scale_height * np.log(1 - u)
    
    # assign them randomly above or below the plane
    z = zabs * np.random.choice([-1, 1], count)
    
    # generate angles uniformly on a circle
    theta = 2 * np.pi * np.random.rand(count)
    
    # return the 3D position
    return (r[0], z[0], theta[0]) if count == 1 else (r, z, theta)

def simulate_mw_distances(count=1, components=[
                                    {"type": "disk", "scale_height": 0.3 * u.kpc, "scale_radius": 2.6 * u.kpc}
                                ]):
    """
        Create a random sample of positions in Milky Way and return distances to these points
        Uses McMillan (2011) for Milky Models
        
        Args:
            count      --> [int]              number of positions to simulate
            components --> [array_like, dict] array of Milky Way components
            
            Each component is a dictionary with information about what type of component it is
            and the scale of heights and radii.
            Example: {"type": "disk", "scale_height": 0.3 * u.kpc, "scale_radius": 2.6 * u.kpc}
            
            Allowed types are 'disk', 'bulge', and 'halo'.
            
        Returns:
            distances  --> [array_like, kpc]  array of distances to each random position
    """
    EARTH_TO_MW_CENTRE = 8.2 * u.kpc
    
    # single component Milky Way model
    if len(components) == 1 and components[0]["type"] == "disk":
        
        # generate random position (r, z, theta) in disk
        radii, heights, angles = random_disk(count, 2.6 * u.kpc, 0.3 * u.kpc)

        # calculate the distance to each coordinate
        distances = np.sqrt(np.square(heights) + np.square(radii) + np.square(EARTH_TO_MW_CENTRE) \
                            - 2 * radii * EARTH_TO_MW_CENTRE * np.cos(angles))
        return distances[0] if count == 1 else distances
    else:
        # TODO: implement multi-component galaxies
        raise NotImplementedError
        
        
def dedt(e, beta, c0):
    """
        Rate of change of eccentricity from Peters (1964) Eq. 5.13
        
        Args:
            e    --> [array_like, unitless] eccentricity
            beta --> [array_like, m^4 / s]  beta constant from Peters (1964), see Eq. 5.9
            c0   --> [array_like, m]        c0 constant from Peters (1964), see Eq. 5.11
    
        Returns:
            dedt --> [array_like, 1 / s]    rate of change of eccentricity
    """
    return - (19 / 12) * (beta / np.power(c0, 4)) * (np.power(e, -29/19) * np.power(1 - np.square(e), 3/2)) \
            / np.power(1 + (121/304) * np.square(e), 1181.0/2299.0)

def c0_peters(a0, e0):
    """
        Find constant c0 from Peters (1964) Eq. 5.11
        
        Args:
            a0 --> [array_like, AU]       initial semi-major axis
            e0 --> [array_like, unitless] initial eccentricity
            
        Returns:
            c0 --> [array_like, AU]       c0 constant from Eq. 5.11    
    """
    return (a0 * (1 - np.square(e0)) * np.power(e0,-12./19) * \
          np.power(1 + (121./304) * np.square(e0), -870./2299)).to(u.AU)

def evolve_eccentricity(e0, a0, m1, m2, t, nsteps=10000):
    """
        Evolve the eccentricity of a binary by numerically integrating using Euler's method
        
        Args:
            e0     --> [array_like, unitless] initial eccentricity
            a0     --> [array_like, AU]       initial semi major axis
            m1     --> [array_like, Msun]     primary star mass
            m2     --> [array_like, Msun]     secondary star mass
            t      --> [array_like, Gyr]      amount of time to inspiral
            nsteps --> [int]                  number of integration steps to take

        Returns:
            e      --> [array_like, unitless] final eccentricity
            c0     --> [array_like, m]        c0 constant from Peters (1964)
        """
    # calculate constants
    c0 = c0_peters(a0, e0).to(u.m).value
    beta = (64./5 * c.G**3 * m1 * m2 * (m1 + m2) / c.c**5).to(u.m**4 / u.s).value

    # initialise values
    e = np.array(e0)

    # find time step
    dt = np.divide(t.to(u.s).value, nsteps)

    # perform Euler's method
    for i in range(nsteps):
        e = e + dedt(e, beta, c0) * dt
    return e, c0 * u.m

def e_to_a(e, c0):
    """
        Convert eccentricity to semi major axis using Peters (1964) Eq. 5.11
        
        Args:
            e  --> [array_like, unitless] eccentricity
            c0 --> [array_like, AU]       c0 constant from Peters (1964), see Eq. 5.11
            
        Returns:
            a  --> [array_like, AU]       semi-major axis
    """
    a = c0 / (1.0 - np.square(e)) * np.power(e, 12./19) * np.power(1 + (121./304 * np.square(e)), 870./2299)
    return a.to(u.AU)