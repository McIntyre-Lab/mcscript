""" Various plotting function """

import numpy as np


def blandAltman(s1, s2, ax=None):
    """ Generate a Bland-Altman plot.

    Arguments:
        :type s1: numpy.array
        :param s1: An array of sample1 data.

        :type s2: numpy.array
        :param s2: An array of sample2 data.

        :param ax:
        :type ax: matplotlib.axes.AxesSubplot

    Returns:
        :returns: If avaiable returns a matplotlib.figure.Figure else adds plot
            to current axis.
        :rtype: matplotlib.figure.Figure

    """
    import matplotlib.pyplot as plt

    # Make sure s1 and s2 are numpy arrays
    s1 = np.asarray(s1)
    s2 = np.asarray(s2)

    # Calculate mean and difference
    mean = (s1 + s2) / 2
    diff = s1 - s2

    # make plot if not axis was provided
    if ax is None:
        fig, ax = plt.subplots(1, 1)

    ax.scatter(mean, diff)
    ax.axhline(0, color='r', ls='--', lw=2)
    ax.set_xlabel('Mean')
    ax.set_ylabel('Difference')

    # Try to return the figure if it exists
    try:
        return fig
    except:
        pass


if __name__ == '__main__':
    pass
