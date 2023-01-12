#!/usr/bin/env python

"""

https://www.statforbiology.com/nonlinearregression/usefulequations#asymptotic_regression_model
"""

import pandas as pd
import toytree
import pymc3 as pm
from velocitree.speciation.regression import linea_f, expon_f, quadr_f, logar_f, asymp_f

SPCOLUMNS = ["hybridizes", "speciesA", "speciesB"]


class BaseLogisticModel:
    """
    Base class for fitting a logistic regression model.
    """
    def __init__(self, tree, spdata, ridata, model="linear"):

        # store input args
        self.tree = tree
        self.spdata = spdata
        self.ridata = ridata
        self.func = model
        # self.check_inputs()

        # get the unique clade indices from spdata.gidx
        self.clades = [i for (i, j) in enumerate(self.spdata.gidx.unique())]

        # to be filled by a subclass 
        self.model = pm.Model()

        # default mcmcparams
        self.sample_kwargs = dict(
            tune=2000,
            draws=2000,
            target_accept=0.95,
            return_inferencedata=False,
            progressbar=True,
        )

        # select the regression model
        functions = {
            "linear": linea_f,
            "exponential": expon_f,
            "quadratic": quadr_f,
            "logarithmic": logar_f,
            "asymptotic_func": asymp_f,
        }
        self.function = functions[self.func]


    def check_inputs(self):
        """
        check tree, spdata and ridata formatting.
        """
        # check that tree is ultrametric
        self.tree = toytree.tree(self.tree)
        # assert self.tree.is_ultrametric(), "tree must be ultrametric"

        # check species assignments to groups, and match names in tree.
        self.spdata = pd.read_csv(self.spdata)
        assert [
            self.spdata.columns.isin(i) for i in SPCOLUMNS
        ], f"spdata column names must include {SPCOLUMNS}"

        # check ridata for names in tree.


    def sample(self, **kwargs):
        """
        Run MCMC sampler
        """
        self.sample_kwargs.update(kwargs)

        # sample posterior, skip burnin
        with self.model:
            trace = pm.sample(**self.sample_kwargs)[1000:]
            stats = pm.summary(trace)
            path = pm.save_trace(
                trace, 
                directory=f"velocitree_trace_{self.model_type}", 
                overwrite=True,
            )
        print(f'trace save to {path}')
        print(stats)        


if __name__ == "__main__":
    pass
