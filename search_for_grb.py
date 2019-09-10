""" 
Given a specific position and time and search radius, 
when is the last time a GRB happened in this location?
"""

from astropy.time import Time

# Search parameters, for testing
ra = 279.472820 
dec = 61.497984
erad = 1 # degree
time = Time(2458728.6798, format="jd")



