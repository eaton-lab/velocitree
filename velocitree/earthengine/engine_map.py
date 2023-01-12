#!/usr/bin/env python

"""
Subclass folium Map to enable earthengine rasters
"""

import ee
import folium


class EngineMap(folium.Map):
    """
    Subclass of folium Map for adding EE images. Here
    self is a folium.Map class instance, and we are adding 
    additional functions to this class type.
    """
    
    def __init__(self, location=None, zoom_start=None, **kwargs):
        # inherit from parent class
        super().__init__(location=location, zoom_start=zoom_start, **kwargs)
    
    
    def add_ee_vector(self, feature, **kwargs):
        """
        Add a google earth engine vector feature to a folium map.
        """
        # create vector feature
        feature = folium.GeoJson(
            data=feature.geometry().getInfo(),
            **kwargs,
        )
        
        # add vector to Map and return self to allow chaining
        self.add_child(feature)
        return self
    
    
    def add_ee_raster(self, image, vis_params=None, **kwargs):
        """
        Add a google earth engine raster layer to a folium map.
        """
        # get visparam kwrgs
        vis_params = {} if vis_params is None else vis_params

        # handle ImageCollections and Images
        if isinstance(image, ee.ImageCollection):
            image = image.mosaic()

        # get the JSON instructions to show image tiles
        map_id_dict = image.getMapId(vis_params)
        
        # get url of the raster tiles
        tiles = map_id_dict['tile_fetcher'].url_format
        
        # create a folium raster layer and add to self
        raster = folium.raster_layers.TileLayer(
            tiles=tiles,
            attr=kwargs.get("attr", "TODO-attr"),
            name=kwargs.get("name", "test"),
            overlay=True,
            control=True,
        )
        
        # add raster to Map and return self to allow chaining
        self.add_child(raster)
        return self
