#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Small selection of classes and functions to make life easier

Created on Sat 12 14:40:58 2019
@author: Richard Boyne rmb115@ic.ac.uk
"""
import matplotlib.pyplot as plt
import sys
import numpy as np
from matplotlib.animation import FuncAnimation
from IPython.display import HTML


class animate_array():
    """
    Animate a 2D graph over time.
    """
    def __init__(self, array, x_points, times=None):
        """
        Parameters
        ----------
        array : 2d array
            Each row has the values at each x_node for that timestep
        x_points : 1d array
            The x_nodes that array corresponds to
        times : 1d array, optional
            The times corresponding to each row in array for text display.
            Defaults to displaying the row index.

        Optional Parameters
        -------------
        (Set these by changing the attributes manually)
        html : bool
            Pass the animation into a html wrapper (for jupyter notebooks)
        frame_interval : int
            The milisecond interval between frames
        frame_skip : int
            The frequency of frames to plot (good for saving animations)
        """

        # validate input
        assert array.ndim == 2, 'must be 2d'
        assert len(x_points) == len(array[0]), 'spacial dimension mismatch'

        # setup data
        self.arr = array
        self.N = len(array)
        self.x = x_points

        # set the times
        if times is None:
            self.pre_string = 'i='
            self.times = np.arange(self.N).astype(str)
        else:
            assert len(times) == self.N, 'The length of the time array is off'
            self.pre_string = 't='
            self.times = times.astype(str)

        # defualt options
        self.html = False
        self.frame_interval = 20
        self.frame_skip = 1

    def blank(self):
        "The blank animation frame"
        self.line.set_data([], [])
        self.text.set_text('')
        return self.line, self.text

    def update(self, i):
        "Plot the ith animation frame"
        self.line.set_data(self.x, self.arr[i, :])
        self.text.set_text(self.pre_string + self.times[i])
        return self.line, self.text

    def set_figure(self):
        "Set the figure seperately so it can be customised if wanted"
        # initialise figure
        self.fig, self.ax = plt.subplots()
        self.line, = self.ax.plot([], [], lw=3, label='Numerical')
        self.text = self.ax.text(0.05, 0.95, '', transform=self.ax.transAxes)

        # set axis limits
        x_r = (self.x[-1] - self.x[0]) * 0.05
        y_r = (self.arr[0, :].max() - self.arr[0, :].min()) * 0.05
        self.ax.set_xlim(self.x[0] - x_r, self.x[-1] + x_r)
        self.ax.set_ylim(self.arr[0, :].min() - y_r,
                         self.arr[0, :].max() + y_r)

    def animate(self):
        "Run the animation, call set_figure first to customise plot"
        if 'fig' not in vars(self):
            self.set_figure()

        # animate
        self.ani = FuncAnimation(self.fig,
                                 self.update,
                                 frames=range(self.N),
                                 interval=self.frame_interval,
                                 blit=True,
                                 init_func=self.blank)

        # html wrapper if wanted
        if self.html:
            plt.close(self.fig)
            self.HTML = HTML(self.ani.to_jshtml())
            return self.HTML

    def save(self, path):
        "Save the animation to a file in path (takes a few seconds to run)"
        # ensure the animation has been made
        if 'ani' not in vars(self):
            self.animate()

        # make the save
        with open(path + '.html', 'w') as f:
            f.write(self.ani.to_jshtml())
        f.close()


class multi_animate():
    """
    Animate a 2D graph over time.
    """
    def __init__(self, arrays, x_points, times=None, labels=None):
        """
        Parameters
        ----------
        array : list of 2d arrays
            each array has rows with values at each x_node for that timestep
        x_points : 1d array
            The x_nodes that array corresponds to
        times : 1d array, optional
            The times corresponding to each row in array for text display.
            Defaults to displaying the row index.
        labels : list
            The legend labels for each array given

        Optional Parameters
        -------------
        (Set these by changing the attributes manually)
        html : bool
            Pass the animation into a html wrapper (for jupyter notebooks)
        frame_interval : int
            The milisecond interval between frames
        frame_skip : int
            The frequency of frames to plot (good for saving animations)
        """

        # validate input
        if type(arrays) is not list and type(arrays) is not tuple:
            arrays = [arrays]
        assert arrays[0].shape[1] == len(x_points),\
            "array dimensions dont match"
        for arr in arrays:
            assert arr.ndim == 2, 'must be 2d arrays'
            assert arr.shape == arrays[0].shape, \
                'all arrays must have same shape'

        # setup data
        self.arrs = arrays
        self.N_arrs = len(arrays)
        self.N = len(arrays[0])
        self.x = x_points

        # set the times
        if times is None:
            self.pre_string = 'i='
            self.times = np.arange(self.N).astype(str)
        else:
            assert len(times) == self.N, 'The length of the time array os off'
            self.pre_string = 't='
            self.times = times.astype(str)

        # set the labels
        if labels is None:
            self.labels = [str('Array %i' % n) for n in range(self.N_arrs)]
        else:
            assert len(labels) == self.N_arrs, 'wrong number of labels'
            self.labels = labels

        # defualt options
        self.html = False
        self.frame_interval = 20
        self.frame_skip = 1
        self.legend = True
        self.figsize = [5, 5]
        self.titles = ['', '', '']
        
        # plotting format string
        if self.N_arrs > 1:
            self.fmts = ['-'] + ['--']*(self.N_arrs - 1)
        else:
            self.fmts = ['-']

        # set axis limits
        x_r = (self.x[-1] - self.x[0]) * 0.05
        y_max = max([arr[0, :].max() for arr in self.arrs])
        y_min = min([arr[0, :].min() for arr in self.arrs])
        y_r = (y_max - y_min) * 0.05
        self.xlims = (self.x[0] - x_r, self.x[-1] + x_r)
        self.ylims = (y_min - y_r, y_max + y_r)

    def blank(self):
        "The blank animation frame"
        for line in self.lines:
            line.set_data([], [])
        self.text.set_text('')
        return (self.text, *self.lines)

    def update(self, i):
        "Plot the ith animation frame"
        for line, arr in zip(self.lines, self.arrs):
            line.set_data(self.x, arr[i, :])
        self.text.set_text(self.pre_string + self.times[i])
        return (self.text, *self.lines)

    def set_figure(self):
        "Set the figure seperately so it can be customised if wanted"
        # initialise figure
        self.fig, self.ax = plt.subplots(figsize=self.figsize)
        self.lines = [self.ax.plot([], [], fmt, lw=3, label=l)[0]
                      for l,fmt in zip(self.labels, self.fmts)]
        self.text = self.ax.text(0.05, 0.95, '', transform=self.ax.transAxes)

        # set axis limits
        self.ax.set_xlim(*self.xlims)
        self.ax.set_ylim(*self.ylims)

        # set the titles
        self.ax.set(title=self.titles[0], xlabel=self.titles[1], ylabel=self.titles[2])
        
        # set legend
        if self.legend:
            self.ax.legend()

    def animate(self):
        "Run the animation, call set_figure first to customise plot"
        self.set_figure()

        # animate
        self.ani = FuncAnimation(self.fig,
                                 self.update,
                                 frames=range(self.N),
                                 interval=self.frame_interval,
                                 blit=True,
                                 init_func=self.blank)

        # html wrapper if wanted
        if self.html:
            plt.close(self.fig)
            self.HTML = HTML(self.ani.to_jshtml())
            return self.HTML

    def save(self, path):
        "Save the animation to a file in path (takes a few seconds to run)"

        # ensure the animation has been made
        if 'ani' not in vars(self):
            self.animate()

        # make the save
        with open(path + '.html', 'w') as f:
            f.write(self.ani.to_jshtml())
        f.close()


def percent_print(i, i_max, interval=1):
    """
    Print a progress bar or percentage value

    Parameters
    ----------
    i : int
        The current iteration number
    i_max : int
        The total number of iterations to complete
    interval : int, optional
        The frequency of updating the progress bar (set high for long
        simulations)

    Returns
    -------
    updated : bool
        Whether the update bar was updated
    """

    if i % interval == 0:

        # print percent
        # sys.stdout.write("\r %.1f %% Done" % (100*i/i_max))

        # print progress bar
        m = int(50 * i/i_max) + 1
        n = 50 - m
        sys.stdout.write("\rProgress |" + "#"*m + " "*n + "|")

        # update the string on screen
        sys.stdout.flush()
        return True

    return False


if __name__ == '__main__':
    x = np.linspace(0, 2*np.pi, 200)
    dts = np.linspace(0, 10, 1000)
    sin_arr = np.array([np.sin(x + dt) for dt in dts])
    cos_arr = np.array([np.cos(x + dt) for dt in dts])
    test = multi_animate([sin_arr, cos_arr], x, labels=['sin', 'cos'])
#    test = animate_array(sin_arr, x, times=np.round(dts, 2))
    test.ylims = [-2, 2] 
    test.animate()
    test.save('wasuuuuup')