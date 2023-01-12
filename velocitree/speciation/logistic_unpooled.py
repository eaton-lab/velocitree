#!/usr/bin/env python

"""
Subclass of BaseLogisticModel to fit a POOLED 
hierarchical regression model.
"""

import pymc3 as pm
from logistic import BaseLogisticModel


class LogisticUnpooled(BaseLogisticModel):
    def __init__(self, tree, spdata, ridata, model="linear"):
        super().__init__(tree, spdata, ridata, model)
        self.model_type = "logistic_pooled"
        self.setup_model()

    def setup_model(self):
        """
        Setup a logistic regression model.        
        """
        with self.model:
            # indexers
            sidx0 = pm.Data("spp_idx0", self.ridata.sidx0.values)
            sidx1 = pm.Data("spp_idx1", self.ridata.sidx1.values)

            # parameters and error
            psi_spp = pm.Normal('𝜓', mu=0, sigma=10, shape=self.tree.ntips)
            beta = pm.Normal('𝛽', mu=0., sigma=10., shape=1)
            
            # linear model prediction
            logit = pm.invlogit(
                self.function(
                    velocity=(beta + psi_spp[sidx0] + psi_spp[sidx1]),
                    distance=self.ridata.distnorm,
                    intercept=0,
                )
            )
            
            # data likelihood (normal distributed errors)
            pm.Bernoulli("y", p=logit, observed=self.ridata.RI)
  


if __name__ == "__main__":

    from generative import Standard

    # generate data using linear
    DATA = Standard(ntips=80, nclades=1, model="linear")
    SUB = DATA.sample_observation(500)

    # fit data using linear 
    MOD = LogisticUnpooled(DATA.tree, DATA.spdata, SUB, model="linear")
    print(MOD.model)
    MOD.sample()
