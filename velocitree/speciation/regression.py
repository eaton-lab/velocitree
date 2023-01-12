#!/usr/bin/env python

"""
Regression functions for the relationship between velocity and gen. distance
"""

from enum import Enum
import numpy as np


def linea_f(velocity, distance, intercept=0):
    "linear function"
    return intercept + velocity * distance

def expon_f(velocity, distance, intercept=0):
    "exponential function"
    return intercept + np.exp(velocity * distance)

def quadr_f(velocity, distance, intercept=0):
    "quadratic function"
    return (intercept * distance) + (velocity * distance ** 2)

def logar_f(velocity, distance, intercept=1e-9):
    "logarithmic function"  
    return intercept + velocity * np.log(distance)

def asymp_f(velocity, distance, intercept=0):
    "asymptotic function"       
    return 1 - (1 - 0) * np.exp(-velocity * distance)


class Models(Enum):
    linear = "linear"
    exponential = "exponential"
    quadratic = "quadratic"
    logarithmic = "logarithmic"
    asymptotic = "asymptotic"


MODEL_DICT = {
    "linear": linea_f,
    "exponential": expon_f,
    "quadratic": quadr_f,
    "logarithmic": logar_f,
    "asymptotic_func": asymp_f,
}


if __name__ == "__main__":
    print(linea_f(velocity=2.0, distance=-2))
