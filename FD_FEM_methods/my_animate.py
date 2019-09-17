#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
<Description>

Created on Fri Nov 30 14:40:58 2018
@author: Richard Boyne rmb115@ic.ac.uk
"""
import matplotlib.pyplot as plt
import sys
from matplotlib.animation import FuncAnimation


class animate_array():
    def __init__(self, array, x_points, dt=1):
        """
        array is 2d with time on axis=0 and function values for x
        points on axis=1
        """
        # setup data
        self.arr = array
        self.N = len(array)
        self.x = x_points
        self.dt = dt

        # setup options
        self.save = ''
        self.interval = 20
        self.times = [dt*n for n in range(self.N)]

    def blank(self):
        self.line.set_data([], [])
        self.text.set_text('')
        return self.line, self.text

    def update(self, i):
        self.line.set_data(self.x, self.arr[i, :])
        self.text.set_text('t={0:.2f}'.format(self.times[i]))
        return self.line, self.text

    def animate(self):
        # initialise figure
        self.fig, self.ax = plt.subplots()
        self.line, = self.ax.plot([], [], lw=3, label='Numerical')
        self.text = self.ax.text(0.75, 0.9, '', transform=self.ax.transAxes)

        # set axis limits
        x_max = self.x[-1]
        x_min = self.x[0]
        y_max = self.arr[0, :].max()
        y_min = min(self.arr[0, :].min(), y_max * -0.01)
        self.ax.set_xlim([x_min - abs(x_min*0.1), x_max + abs(x_max*0.1)])
        self.ax.set_ylim([y_min - abs(y_min*0.1), y_max + abs(y_max*0.1)])

        # animate
        self.ani = FuncAnimation(self.fig,
                                 self.update,
                                 frames=range(self.N),
                                 interval=self.interval,
                                 blit=True,
                                 init_func=self.blank)
        if self.save != "":
            self.ani.save(self.save)


class dual_animate():
    def __init__(self, array1, array2, x_points, dt=1):
        """
        array is 2d with time on axis=0 and function values for x
        points on axis=1
        """
        # setup data
        self.arr1 = array1
        self.arr2 = array2
        self.N = len(array1)
        self.x = x_points
        self.dt = dt

        # setup options
        self.save = ''
        self.interval = 20
        self.times = [dt*n for n in range(self.N)]
        self.xlims = [self.x[0], self.x[-1]]
        self.ylims = [self.arr1[0, :].min(), self.arr1[0, :].max()]

    def blank(self):
        self.line1.set_data([], [])
        self.line2.set_data([], [])
        self.text.set_text('')
        return self.line1, self.line2, self.text

    def update(self, i):
        self.line1.set_data(self.x, self.arr1[i, :])
        self.line2.set_data(self.x, self.arr2[i, :])
        self.text.set_text('t={0:.2f}'.format(self.times[i]))
        return self.line1, self.line2, self.text

    def animate(self):
        # initialise figure
        self.fig, self.ax = plt.subplots()
        self.line1, = self.ax.plot([], [], lw=2, label='Numerical')
        self.line2, = self.ax.plot([], [], '--', lw=2, label='Analtical')
        self.text = self.ax.text(0.75, 0.9, '', transform=self.ax.transAxes)

        # set axis limits
        self.ax.set_xlim(*self.xlims)
        self.ax.set_ylim(*self.ylims)
        self.ax.legend(loc=2)

        # animate
        self.ani = FuncAnimation(self.fig,
                                 self.update,
                                 frames=range(self.N),
                                     interval=self.interval,
                                 blit=True,
                                 init_func=self.blank)
        if self.save != "":
            self.ani.save(self.save)


def percent_print(i, i_max, interval=20):
    if i % interval == 0:
        # print percent
        # sys.stdout.write("\r %.1f %% Done" % (100*i/i_max))

        # print progress bar
        m = int(50 * i/i_max) + 1
        n = 50 - m
        sys.stdout.write("\rProgess |" + "#"*m + " "*n + "|")

        # update the string on screen
        sys.stdout.flush()
        return


if __name__ == '__main__':
    import numpy as np
    x = np.linspace(0, 2*np.pi, 200)
    dts = np.linspace(0, 10, 1000)
    C = np.array([np.sin(x + dt) for dt in dts])
    test = animate_array(C, x)
    test.animate()
