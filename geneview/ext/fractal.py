# -*- coding: utf-8 -*-
from math import log
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt


def _iter_point(c, iter_num=100, escape_radius=2):
    z = c
    n = 1
    for i in range(1, iter_num):  # 最多迭代100次
        n = i
        if abs(z) > escape_radius:
            break  # 半径大于该值则认为逃逸
        z = z * z + c

    return n  # 迭代次数


def _smooth_iter_point(c, iter_num=20, escape_radius=10):
    z = c
    n = 1
    for i in range(1, iter_num):  # 最多迭代100次
        n = i
        if abs(z) > escape_radius:
            break  # 半径大于该值则认为逃逸
        z = z * z + c

    mu = n - log(log(abs(z), 2), 2) if abs(z) > 2.0 else n
    return mu  # 迭代次数


def _mandelbrot(cx, cy, d, n=200, iter_point_func=_iter_point,
                iter_num=100, escape_radius=2):
    """                                                                                  
    绘制点 (cx, cy) 附近正负d的范围的 Mandelbrot
    """
    x0, x1, y0, y1 = cx - d, cx + d, cy - d, cy + d
    y, x = np.ogrid[y0:y1:n * 1j, x0:x1:n * 1j]
    c = x + y * 1j

    m = np.frompyfunc(iter_point_func, 3, 1)(c, iter_num, escape_radius).astype(np.float)
    return m


def mandelbrot():
    plt.subplot(231)
    plt.imshow(_mandelbrot(-0.5, 0, 1.5),
               cmap=cm.Blues_r,
               extent=[-0.5 - 1.5, -0.5 + 1.5, -1.5, 1.5])
    plt.gca().set_axis_off()

    x, y = 0.27322626, 0.595153338
    for i in range(2, 7):
        plt.subplot(230 + i)
        _mandelbrot(x, y, 0.2 ** (i - 1))
        d = 0.2 ** (i - 1)
        plt.imshow(_mandelbrot(x, y, d),
                   cmap=cm.Blues_r, extent=[x - d, x + d, y - d, y + d])
        plt.gca().set_axis_off()

    plt.subplots_adjust(0.02, 0, 0.98, 1, 0.02, 0)


def mandelbrot_smooth():
    x, y, d = -0.5, 0, 1.5
    escape_radius, iter_num, N = 10, 20, 300

    plt.subplot(121)
    plt.imshow(_mandelbrot(x, y, d, N, _iter_point,
                           iter_num, escape_radius),
               cmap=cm.Blues_r, extent=[x - d, x + d, y - d, y + d])
    plt.gca().set_axis_off()

    plt.subplot(122)
    plt.imshow(_mandelbrot(x, y, d, N, _smooth_iter_point,
                           iter_num, escape_radius),
               cmap=cm.Blues_r, extent=[x - d, x + d, y - d, y + d])

    plt.subplots_adjust(0.02, 0, 0.98, 1, 0.02, 0)
    plt.gca().set_axis_off()


if __name__ == "__main__":
    # mandelbrot()
    mandelbrot_smooth()
    plt.show()
