# Modified demo file on animating matrix/imshow with matplotlib to represent forest-fire

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors
#from random import randint
fig = plt.figure()

#%%
from matplotlib.colors import LinearSegmentedColormap
colors = [(0, 0, 0), (0, 1, 0), (1, 1, 1)]  # Black -> Green -> White
n_bins = 3  # Discretizes the interpolation into bins
cmap_name = 'ForestFire'
# Create the colormap
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

#plt.figure(figsize = (15,15)) 

dimension = 100

#def f(x, y):
#    return np.sin(x) + np.cos(y)

#x = np.linspace(0, 2 * np.pi, 120)
#y = np.linspace(0, 2 * np.pi, 100).reshape(-1, 1)
randArr = np.random.randint(-1, 2, size = (dimension,dimension))


im = plt.imshow(randArr, cmap = cm , animated=True)
norm = matplotlib.colors.BoundaryNorm(np.arange(-1.5,2,1), cm.N)
plt.colorbar(ticks=np.linspace(-1,1,3))

def updatefig(*args):
    #global x, y
    #x += np.pi / 15.
    #y += np.pi / 20.
    global randArr
    randArr = np.random.randint(-1, 2, size = (dimension,dimension))
    im.set_array(randArr)
    return im,



ani = animation.FuncAnimation(fig, updatefig, frames = 100, interval=100, blit=True)

HTML(ani.to_html5_video())
#dpi = 300
#writer = animation.writers['ffmpeg'](fps=30)
#ani.save('test.mp4',writer=writer,dpi=dpi)