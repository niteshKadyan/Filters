from matplotlib.widgets import Button
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import random
from math import *

landmarks  = [[20.0, 20.0], [80.0, 80.0], [20.0, 80.0], [80.0, 20.0]]
world_size = 100.0
N = 1000

class robot:

    def __init__(self):
        self.x = random.random() * world_size
        self.y = random.random() * world_size
        self.orientation = random.random() * 2.0 * pi
        self.forward_noise = 0.0;
        self.turn_noise    = 0.0;
        self.sense_noise   = 0.0;
    
    def set(self, new_x, new_y, new_orientation):
        if new_x < 0 or new_x >= world_size:
            raise ValueError, 'X coordinate out of bound'
        if new_y < 0 or new_y >= world_size:
            raise ValueError, 'Y coordinate out of bound'
        if new_orientation < 0 or new_orientation >= 2 * pi:
            raise ValueError, 'Orientation must be in [0..2pi]'
        self.x = float(new_x)
        self.y = float(new_y)
        self.orientation = float(new_orientation)
    
    
    def set_noise(self, new_f_noise, new_t_noise, new_s_noise):
        # makes it possible to change the noise parameters
        # this is often useful in particle filters
        self.forward_noise = float(new_f_noise);
        self.turn_noise    = float(new_t_noise);
        self.sense_noise   = float(new_s_noise);
    
    
    def sense(self):
        Z = []
        for i in range(len(landmarks)):
            dist = sqrt((self.x - landmarks[i][0]) ** 2 + (self.y - landmarks[i][1]) ** 2)
            dist += random.gauss(0.0, self.sense_noise)
            Z.append(dist)
        return Z
    
    
    def move(self, turn, forward):
        if forward < 0:
            raise ValueError, 'Robot cant move backwards'         
        
        # turn, and add randomness to the turning command
        orientation = self.orientation + float(turn) + random.gauss(0.0, self.turn_noise)
        orientation %= 2 * pi
        
        # move, and add randomness to the motion command
        dist = float(forward) + random.gauss(0.0, self.forward_noise)
        x = self.x + (cos(orientation) * dist)
        y = self.y + (sin(orientation) * dist)
        x %= world_size    # cyclic truncate
        y %= world_size
        
        # set particle
        res = robot()
        res.set(x, y, orientation)
        res.set_noise(self.forward_noise, self.turn_noise, self.sense_noise)
        return res
    
    def Gaussian(self, mu, sigma, x):  
        # calculates the probability of x for 1-dim Gaussian with mean mu and var. sigma
        return exp(- ((mu - x) ** 2) / (sigma ** 2) / 2.0) / sqrt(2.0 * pi * (sigma ** 2))
    
    
    def measurement_prob(self, measurement):
        
        # calculates how likely a measurement should be    
        prob = 1.0;
        for i in range(len(landmarks)):
            dist = sqrt((self.x - landmarks[i][0]) ** 2 + (self.y - landmarks[i][1]) ** 2)
            prob *= self.Gaussian(dist, self.sense_noise, measurement[i])
        return prob
    
    def __repr__(self):
        return '[x=%.6s y=%.6s orient=%.6s]' % (str(self.x), str(self.y), str(self.orientation))

class plot:
    ind = 0
    myrobot = None
    Z = []
    particles = []
    prevLoc = []
    def __init__(self):
        self.ind = 0
        self.myrobot = robot()
        self.Z = self.myrobot.sense()
        self.particles = particles
        print "Initialized"
    def next(self, event):
        self.ind += 1
        
        self.prevLoc.append([self.myrobot.x, self.myrobot.y])
        
        self.myrobot = self.myrobot.move(0.1, 5)
        self.Z = self.myrobot.sense()

        p2 = []
        for i in range(N):
            p2.append(self.particles[i].move(0.1, 5))
        self.particles = p2

        w = []
        for i in range(N):
            w.append(self.particles[i].measurement_prob(self.Z))

        p3 = []
        index = int(random.random() * N)
        beta = 0.0
        mw = max(w)
        for i in range(N):
            beta += random.random() * 2.0 * mw
            while beta > w[index]:
                beta -= w[index]
                index = (index + 1) % N
            p3.append(self.particles[index])
        self.particles = p3
        
        xdata = []
        ydata = []
        for i in range(N):
            xdata.append(self.particles[i].x)
            ydata.append(self.particles[i].y)
        l.set_xdata(xdata)
        l.set_ydata(ydata)
        plt.draw()
        
        xdata1 = []
        ydata1 = []
        for i in range(len(self.prevLoc)):
            xdata1.append(self.prevLoc[i][0])
            ydata1.append(self.prevLoc[i][1])
        l1.set_xdata(xdata1)
        l1.set_ydata(ydata1)
        plt.draw()
        
        print self.ind
        print self.myrobot
        print self.Z
        #print self.prevLoc
        #print self.particles



particles = []
p = []
for i in range(N):
    r = robot()
    r.set_noise(0.1, 0.1, 5.0)
    particles.append(r)
    p.append([r.x, r.y])

#print particles
    
l, = plt.plot(*zip(*p), marker='o', color='r', ls='')
l1, = plt.plot(*zip(*p), marker='+', color='g', ls='')

ax = plt.subplot(111)
plt.subplots_adjust(bottom=0.2)
img = mpimg.imread('F:\PIS Project\line_ex.png')
plt.imshow(img)
callback = plot()
axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
bnext = Button(axnext, 'Next')
bnext.on_clicked(callback.next)
plt.show()

