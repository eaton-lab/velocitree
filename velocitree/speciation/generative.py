#!/usr/env/bin python

"""
Generate data under a regression model with known parameters
for validating model setup and performing power analyses.
"""

from typing import List, Optional
import itertools
import toytree
import numpy as np
import pandas as pd
import scipy.special
from loguru import logger
from velocitree.speciation.regression import Models, MODEL_DICT
from velocitree.speciation.regression import (
    linea_f, expon_f, quadr_f, logar_f, asymp_f
)
# import toyplot
# import pymc3 as pm


class GenerativeBase:
    def __init__(self, tree, clade_idxs, model:Models, seed:Optional[int]):

        # set random seed
        self.rng = np.random.default_rng(seed=seed)

        # store user input
        self.tree = toytree.tree(tree).mod.node_scale_root_height(1.0)
        self.clade_idxs: List[int] = clade_idxs
        self.function = MODEL_DICT[model]

        # attrs to fill (type hints)
        self.ridata: pd.DataFrame = None
        self.subdata: pd.DataFrame = None
        self.group_sizes: dict = {}
        self.groups: list = []
        self.params: dict = {}

        self.get_groups()
        self.get_crosses()
        self.get_true_param_dists()
        self.get_species_data()
        self.get_velocities()
        self.get_logit()

    def get_groups(self):
        """
        Expand clade_idxs to the tip idxs
        """
        self.group_sizes = {
            i: len(self.tree.idx_dict[j]) 
            for (i, j) in enumerate(self.clade_idxs)
        }
        self.groups = np.concatenate([
            np.repeat(i, self.group_sizes[i]) 
            for i in self.group_sizes
        ])

    def get_crosses(self):
        """
        Get a dataframe of all possible pairs and a subsampled one 
        representing a set of observations.
        """
        # get all combinations of two tips
        tipsa, tipsb = zip(*itertools.combinations(range(self.tree.ntips), 2))

        # get tip idxs and genetic distances between them
        self.ridata = pd.DataFrame({
            "sidx0": tipsa,
            "sidx1": tipsb,
            "dist": [
                get_dist(self.tree, i, j) / 2. for (i, j) in zip(tipsa, tipsb)
            ],
        })

    def get_true_param_dists(self, random=False):
        """
        Sample a set of TRUE praam distributions from which data will
        be generated.
        """
        if random:
            self.params['𝛽'] = self.rng.uniform(3, 5)
            self.params['𝜓_means'] = [
                self.rng.uniform(0, 2) for i in self.clade_idxs
            ]
            self.params['𝜓_stds'] = [
                abs(self.rng.normal(0, 1)) for i in self.clade_idxs
            ]
        else:
            self.params['𝛽'] = 10
            self.params['𝜓_means'] = np.linspace(-0.5, 0.5, len(self.clade_idxs))
            self.params['𝜓_stds'] = np.linspace(0.05, 0.25, len(self.clade_idxs))

    def get_species_data(self):
        """
        Generate species-specific psi values by drawing from the
        distribution that describes them. This function varies in 
        the different model classes (pooled, unpooled, partpooled).
        """
        self.spdata = pd.DataFrame({
            "gidx": self.groups,
            "𝜓": self.rng.normal(
                [self.params['𝜓_means'][i] for i in self.groups],
                [self.params['𝜓_stds'][i] for i in self.groups], 
            )
        })

    def get_velocities(self):
        """
        Generate velocities of observed species crosses to use for 
        generating test data.
        """
        # add RI=0 for all within-species comparisons
        nulls = pd.DataFrame({
            'sidx0': range(self.tree.ntips),
            'sidx1': range(self.tree.ntips),
            'dist': 0.,
        })
        self.ridata = pd.concat([self.ridata, nulls], ignore_index=True)

        # sort RI by species indices
        self.ridata = (
            self.ridata
                .sort_values(by=['sidx0', 'sidx1'])
                .reset_index(drop=True)
            )

        # get normalized dists (dists will be approx -2 to 2, instead of 0-1).
        self.ridata["distnorm"] = (
            (self.ridata.dist - self.ridata.dist.mean()) / self.ridata.dist.mean())       

        # compute joint velocities
        self.ridata['joint_velo'] = (
            self.rng.normal(
                self.params['𝛽'] + 
                self.spdata['𝜓'][self.ridata.sidx0].values +
                self.spdata['𝜓'][self.ridata.sidx1].values
        ))

    def get_logit(self):
        """
        Calculate the log odds of hybridizing by implementing a 
        regression function between genetic distance and the joint
        velocity of RI for a pair of species. The options for model
        regressors are: "linear", "logarithmic", "asymptotic", 
        "exponential" and "quadratic".
        """
        # add log-odds of hybridizing
        values = self.function(
            velocity=self.ridata.joint_velo,
            distance=self.ridata.dist,            
            # distance=self.ridata.distnorm,
        )
        # try more normalizations
        values = values - values.mean() / values.std()
        self.ridata['expit'] = scipy.special.expit(values)

    def sample_observations(self, nsamples=None):
        """
        Randomly sample a set of observations with RI values sampled
        from the expit probabilities.
        """
        self.ridata["RI"] = np.random.binomial(1, self.ridata['expit'])
        if nsamples:
            return self.ridata.sample(nsamples).reset_index(drop=True).copy()
        return self.ridata.copy()


class UserTree(GenerativeBase):
    """
    Generate crossing dataset on a user supplied tree.
    """
    def __init__(self, tree, clade_idxs, model, seed):
        super().__init__(tree, clade_idxs, model, seed)


class RandomTree(GenerativeBase):
    """
    Generate a crossing dataset from a random tree.
    """
    def __init__(self, ntips=80, nclades=None, model="linear", seed=None):

        # store args
        self.ntips = ntips
        self.nclades = nclades

        # to be filled
        self.tree: toytree.Toytree.ToyTree = None
        self.clade_idxs: List[int] = None

        # init funcs
        self.get_bd_tree()
        self.get_clade_idxs()
        super().__init__(self.tree, self.clade_idxs, model, seed)


    def get_bd_tree(self):
        """
        generate a random bdtree for ntips and scale to 1.0
        """
        self.tree = toytree.rtree.bdtree(
            ntips=self.ntips,
            seed=666,
        )

    def get_clade_idxs(self):
        """
        Get objects for indexing individuals in groups/clades.
        """
        # generate group assignments. 
        node = self.tree.treenode
        self.clade_idxs = [node.idx]
        for _ in range(2, self.nclades + 1):
            sidxs = sorted(
                self.clade_idxs, 
                key=lambda x: len(self.tree.idx_dict[x].children)
            )
            nidx = sidxs[0]
            self.clade_idxs.remove(nidx)
            self.clade_idxs.extend([
                i.idx for i in self.tree.idx_dict[nidx].children
            ])



def get_dist(tree, idx0, idx1):
    """
    TODO: use toytree func.; 
    Returns the genetic distance between two nodes on a tree.
    """
    dist = tree.treenode.get_distance(
        tree.idx_dict[idx0], 
        tree.idx_dict[idx1],
    )
    return dist



if __name__ == "__main__":

    from viz import dist_RI_scatterplot, heatmap_tree_plot

    # generate a bdtree dataset
    TEST = RandomTree(ntips=80, nclades=4, model="linear")
    print(TEST.params)
    print(TEST.spdata.head())

    SUBSAMPLE = TEST.sample_observations(1000)
    print(SUBSAMPLE.head())

    # # visualization
    dist_RI_scatterplot(SUBSAMPLE)
    heatmap_tree_plot(TEST.tree, TEST.clade_idxs, SUBSAMPLE)
