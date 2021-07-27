""" Utilities for dealing with galaxy catalogs. """

class Catalog(object):
    def __init__(self, data, **kwargs):
        self.data = data
        self.kwargs = kwargs
