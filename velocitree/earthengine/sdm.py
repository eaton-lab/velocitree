#!/usr/bin/env python

"""Generate species distribution models for all species in data.

- https://eaton-lab.org/hack-the-planet/tutorials/25.2-folium.html
- https://developers.google.com/earth-engine/apidocs/ee-classifier-amnhmaxent#colab-python
"""

import ee
# from engine_map import EngineMap


# make user login 
ee.Authenticate()
ee.Initialize()


class SpeciesDistributionModel:
    def __init__(self, species):
        self.species = species

    def get_filtered_sproc_points(self):
        pass

    def sample_training_data(self):
        # Create some sample species presence/absence training data.
        training_data = ee.FeatureCollection([
            # Species present points.
            ee.Feature(ee.Geometry.Point([-122.39567, 38.02740]), {'presence': 1}),
            ee.Feature(ee.Geometry.Point([-122.68560, 37.83690]), {'presence': 1}),
            
            # Species absent points.
            ee.Feature(ee.Geometry.Point([-122.59755, 37.92402]), {'presence': 0}),
            ee.Feature(ee.Geometry.Point([-122.47137, 37.99291]), {'presence': 0}),
            ee.Feature(ee.Geometry.Point([-122.52905, 37.85642]), {'presence': 0}),
            ee.Feature(ee.Geometry.Point([-122.03010, 37.66660]), {'presence': 0})
        ])

    # Import a Landsat 8 image and select the reflectance bands.
    image = (
        ee.Image('LANDSAT/LC08/C01/T1_SR/LC08_044034_20200606')
        .select(['B[0-9]*'])
    )

    # Sample the image at the location of the points.
    training = image.sampleRegions(**{
        'collection': training_data,
        'scale': 30
    })

    # Define and train a Maxent classifier from the image-sampled points.
    classifier = ee.Classifier.amnhMaxent().train(**{
        'features': training,
        'classProperty': 'presence',
        'inputProperties': image.bandNames()
    })

    # Classify the image using the Maxent classifier.
    image_classified = image.classify(classifier)    



if __name__ == "__main__":
    pass