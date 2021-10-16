from geobench import Sneddon
from numpy import linspace


def test_version():
    problem = Sneddon()
    # Compute stress along y axis with x=0.
    x = 0
    y = linspace(0, 20, num=41)
    sigma_x, sigma_y = problem.stress(x, y)

    # x = linspace(-l / 2, l / 2, 20)
    # w = problem.width(x)
