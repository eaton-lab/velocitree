#!/usr/bin/env python

"""
Subclass of BaseLogisticModel to fit a PARTIALLY-POOLED
hierarchical regression model.
"""

import pymc3 as pm
from velocitree.speciation.logistic import BaseLogisticModel


class LogisticPartPooled(BaseLogisticModel):
    def __init__(self, tree, spdata, ridata, model="linear"):
        super().__init__(tree, spdata, ridata, model)
        self.model_type = f"partpooled_{self.func}"
        self.setup_model()

    def setup_model(self):
        """
        Setup a logistic regression model.        
        """
        with self.model:
            # indexers
            sidx0 = pm.Data("spp_idx0", self.ridata.sidx0.values)
            sidx1 = pm.Data("spp_idx1", self.ridata.sidx1.values)
            gidx = pm.Data("gidx", self.spdata.gidx.values)

            # parameters and error
            psi_mean = pm.Normal(
                '𝜓_mean', mu=0., sigma=10., shape=len(self.clades))
            psi_std = pm.HalfNormal('𝜓_std', 5., shape=len(self.clades))
            psi_offset = pm.Normal(
                '𝜓_offset', mu=0, sigma=1., shape=self.tree.ntips)
            psi_spp = pm.Deterministic(
                '𝜓', psi_mean[gidx] + psi_std[gidx] * psi_offset)
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

    from velocitree.speciation.generative import RandomTree

    # generate data using linear
    MODEL = RandomTree(ntips=80, nclades=4, model="linear")
    DATA = MODEL.sample_observations(500)
    print("OBSERVATIONS:\n{}".format(DATA))

    # fit data using linear 
    FIT = LogisticPartPooled(MODEL.tree, MODEL.spdata, DATA, model="linear")
    print(FIT.model)
    FIT.sample()
    print("TRUE PARAMS:\n{}".format(MODEL.params))
