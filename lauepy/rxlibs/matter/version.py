##############################################################################
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""Definition of __version__ and __date__ for matter.
"""

__id__ = "$Id: version.py 2825 2009-03-09 04:33:12Z juhas $"

# obtain version information
# from pkg_resources import get_distribution


#__version__ = get_distribution('matter').version
#shortcut this for now
__version__ = 0.1

# we assume that tag_date was used and __version__ ends in YYYYMMDD
#__date__ = __version__[-8:-4] + '-' + \
#           __version__[-4:-2] + '-' + __version__[-2:]


# End of file
