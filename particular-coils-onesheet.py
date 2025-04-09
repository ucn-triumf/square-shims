#!/usr/bin/env python3

# Fri May 24 10:11:45 CDT 2019 Jeff added this line.

# Tue Feb 11 13:43:43 CST 2020 Jeff taking original patch.py and
# updating to solve the zero mode issue.  Will now update to use the
# patchlib submodule.

# Fri Feb 14 11:45:17 CST 2020 Jeff speeding up code by improving the
# sensorarray class to use numpy structures.

# Sat Aug 1 12:40:53 CDT 2020 Jeff added command line options and
# improved graphs

#Fri July 19 2024, Modeste added the opting to draw the msr and coroplast set up at line 1056 run with options -p.

# Wed Mar 26 16:46:14 PDT 2025 Jeff and Tahereh, fixing all the mutual mistakes of everyone involved in this mess of a project.
# We noticed also that about half the coils are wound in the wrong direction in Modeste's code.
# We are going to use a convention that a positive current in a coil generates a magnetic field
# at its centre in the direction pointing inward, i.e. into the room.

from scipy.constants import mu_0, pi
import numpy as np
from patchlib.patch import *
from Pis.Pislib import *
from dipole import *

#from pipesfitting import *

#from arduino_current_controller_routines import *


from optparse import OptionParser

parser = OptionParser()

parser.add_option("-s", "--nsensors", dest="nsensors", default=3,
                  help="ns where total sensor axes is s = 3*ns^3")

parser.add_option("-l", "--ell", dest="l", default=2,
                  help="l for spherical harmonic")

parser.add_option("-m", "--em", dest="m", default=0,
                  help="m for spherical harmonic")

parser.add_option("-M", "--matrices", dest="matrices", default=False,
                  action="store_true",
                  help="show matrices")

parser.add_option("-d", "--dipole", dest="dipole", default=False,
                  action="store_true",
                  help="use dipole field")

parser.add_option("-t", "--traces", dest="traces", default=False,
                  action="store_true",
                  help="show 3D view of coils and sensors")

parser.add_option("-r", "--residuals", dest="residuals", default=False,
                  action="store_true",
                  help="show residuals")

parser.add_option("-z", "--zoom", dest="zoom", default=False,
                  action="store_true",
                  help="zoom to ROI")

parser.add_option("-a", "--axes", dest="axes", default=False,
                  action="store_true",
                  help="make graphs along axes")

parser.add_option("-i", "--incells", dest="incells", default=False,
                  action="store_true",
                  help="ROI for statistics is in EDM cells")
parser.add_option("-p", "--makeplots", dest="makeplots", default=False,
                  action="store_true",
                  help="Make plots of walls")
parser.add_option("-w", "--wiggle", dest="wiggle",
                  action="store_true",
                  default=False, help="wiggle each point")

#d=dipole(1.2,0,0,0,0,100000)  # dipole1
#d=dipole(0,0,1.2,0,0,1)  # dipole2
d=dipole(0,0,1.2,1,0,0)  # dipole3

(options,args)=parser.parse_args()

l=int(options.l)
m=int(options.m)
sp=scalarpotential(l,m)
print("Sigma in spherical coordinates is %s"%sp.Sigma_spherical)
print("Sigma in cartesian coordinates is %s"%sp.Sigma)

print("Pix is %s"%sp.Pix)
print("Piy is %s"%sp.Piy)
print("Piz is %s"%sp.Piz)

if(options.dipole):
    bxtarget=d.bx
    bytarget=d.by
    bztarget=d.bz
else:
    bxtarget=sp.fPix
    bytarget=sp.fPiy
    bztarget=sp.fPiz

# Setup our coilset

myset=coilset()

# positions of faces (positive values -- faces will be at plus and
# minus of these)
a=2.2 # approximate position of layer 6 of MSR
xface=a/2 # m
yface=a/2 # m
zface=a/2 # m

#picture frame separated by 40mm

# set up the East wall
y1=300*.001
z1=300*.001
x1=xface
point1=(x1,y1,z1)
y2=y1
z2=600*0.001  # Jeff changed to fit on 4' sheet
x2=xface
point2=(x2,y2,z2)
y3=600*.001 # guess
z3=z2
x3=xface
point3=(x3,y3,z3)
y4=y3
z4=z1
x4=xface
point4=(x4,y4,z4)
points_ul=(point1,point2,point3,point4)
points_ul=np.array(points_ul)

# Now add mirror images of these
point1=(x1,-y1,z1)
point2=(x4,-y4,z4)
point3=(x3,-y3,z3)
point4=(x2,-y2,z2)
points_ur=(point1,point2,point3,point4)
points_ur=np.array(points_ur)

point1=(x1,-y1,-z1)
point2=(x2,-y2,-z2)
point3=(x3,-y3,-z3)
point4=(x4,-y4,-z4)
points_lr=(point1,point2,point3,point4)
points_lr=np.array(points_lr)

point1=(x1,y1,-z1)
point2=(x4,y4,-z4)
point3=(x3,y3,-z3)
point4=(x2,y2,-z2)
points_ll=(point1,point2,point3,point4)
points_ll=np.array(points_ll)

# now the central coil
x1=xface
y1=(235+5+20)*.001
z1=(157+5+98)*.001
point1=(x1,y1,z1)
point2=(x1,y1,-z1)
point3=(x1,-y1,-z1)
point4=(x1,-y1,z1)
points_c=(point1,point2,point3,point4)
points_c=np.array(points_c)

# now the middle left coil
x1=xface
y1=(235+5+20+40)*0.001
z1=(300-40)*.001
point1=(x1,y1,z1)
x2=xface
y2=600*.001 # guess
z2=z1
point2=(x2,y2,z2)
x3=xface
y3=y2
z3=-z2
point3=(x3,y3,z3)
x4=xface
y4=y1
z4=-z1
point4=(x4,y4,z4)
points_ml=(point1,point2,point3,point4)
points_ml=np.array(points_ml)

# now the middle right coil -- reflect and wind in same direction
point1=(x1,-y1,z1)
point2=(x4,-y4,z4)
point3=(x3,-y3,z3)
point4=(x2,-y2,z2)
points_mr=(point1,point2,point3,point4)
points_mr=np.array(points_mr)

# now the upper central coil
x1=xface
y1=(235+5+20)*0.001
z1=(300)*.001
point1=(x1,y1,z1)
x2=xface
y2=-y1
z2=z1
point2=(x2,y2,z2)
x3=xface
y3=y2
z3=600*.001
point3=(x3,y3,z3)
x4=xface
y4=y1
z4=z3
point4=(x4,y4,z4)
points_uc=(point1,point2,point3,point4)
print('points_uc',points_uc)
points_uc=np.array(points_uc)

# now the lower central coil -- reflect and wind in same direction
point1=(x1,y1,-z1)
point2=(x4,y4,-z4)
point3=(x3,y3,-z3)
point4=(x2,y2,-z2)
points_lc=(point1,point2,point3,point4)
points_lc=np.array(points_lc)

 
'''
myset.add_coil(points_ur)
myset.add_coil(points_ul)
myset.add_coil(points_ll)
myset.add_coil(points_lr)
myset.add_coil(points_c)
myset.add_coil(points_mr)
myset.add_coil(points_ml)
myset.add_coil(points_uc)
myset.add_coil(points_lc)
'''

# Modeste's suggestion on how to build it

# myset.add_coil(points_ur)
# myset.add_coil(points_uc)
# myset.add_coil(points_ul)
# myset.add_coil(points_mr)
# myset.add_coil(points_c)
# myset.add_coil(points_ml)
# myset.add_coil(points_lr)
# myset.add_coil(points_lc)
# myset.add_coil(points_ll)

# How it was actually built (Jeff and Tahereh convention)

myset.add_coil(points_ur) # coil 0
myset.add_coil(points_uc) # coil 1
myset.add_coil(points_ul) # coil 2
myset.add_coil(points_mr) # coil 3
myset.add_coil(points_c)  # coil 4
myset.add_coil(points_ml) # coil 5
myset.add_coil(points_lr) # coil 6
myset.add_coil(points_lc) # coil 7
myset.add_coil(points_ll) # coil 8


# now reflect them all to the other face: xface -> -xface
def reflect_x(points):
    newpoints=np.copy(points)
    newpoints[:,0]=-newpoints[:,0]
    newpoints=np.flip(newpoints,0) # wind them in the opposite direction
    return newpoints


# Jeff and Tahereh note that these names should likely be changed to reflect how someone
# standing inside the MSR would understand them.  But we will not fix that problem
# right now.
    
oside_ur=reflect_x(points_ur)
oside_ul=reflect_x(points_ul)
oside_ll=reflect_x(points_ll)
oside_lr=reflect_x(points_lr)
oside_c=reflect_x(points_c)
oside_ml=reflect_x(points_ml)
oside_mr=reflect_x(points_mr)
oside_uc=reflect_x(points_uc)
oside_lc=reflect_x(points_lc)


myset.add_coil(oside_ul)
myset.add_coil(oside_uc)
myset.add_coil(oside_ur)
myset.add_coil(oside_ml)
myset.add_coil(oside_c)
myset.add_coil(oside_mr)
myset.add_coil(oside_ll)
myset.add_coil(oside_lc)
myset.add_coil(oside_lr)



# Phew -- now onto the front/back  (south/north walls)

# coil 18

z1=(260-5-40)*.001
x1=(220-5)*.001
y1=yface
point1=(x1,y1,z1)
z2=600*.001 #guess
x2=x1
y2=yface
point2=(x2,y2,z2)
z3=z2
x3=600*.001 #guess
y3=yface
point3=(x3,y3,z3)
z4=z1
x4=x3
y4=yface
point4=(x4,y4,z4)
side_ur=(point4,point3,point2,point1) # Jeff and Tahereh changed the winding direction for these coils
side_ur=np.array(side_ur)

# now reflect around
point1=(-x1,y1,z1)
point2=(-x4,y4,z4)
point3=(-x3,y3,z3)
point4=(-x2,y2,z2)
side_ul=np.array((point4,point3,point2,point1))

point1=(-x1,y1,-z1)
point2=(-x2,y2,-z2)
point3=(-x3,y3,-z3)
point4=(-x4,y4,-z4)
side_ll=np.array((point4,point3,point2,point1))

point1=(x1,y1,-z1)
point2=(x4,y4,-z4)
point3=(x3,y3,-z3)
point4=(x2,y2,-z2)
side_lr=np.array((point4,point3,point2,point1))

# central coil
z1=(170+5)*.001
y1=yface
x1=(-221-5)*.001
point1=(x1,y1,z1)
point2=(x1,y1,-z1)
point3=(-x1,y1,-z1)
point4=(-x1,y1,z1)
side_c=np.array((point1,point2,point3,point4))

# middle right coil
x1=(221+5+40)*.001
y1=yface
z1=(170+5)*.001
point1=(x1,y1,z1)
x2=600*.001 #guess
y2=yface
z2=z1
point2=(x2,y2,z2)
point3=(x2,y2,-z2)
point4=(x1,y1,-z1)
side_mr=np.array((point4,point3,point2,point1))  # Jeff/Tahereh

# reflect it to middle left coil
point1=(-x1,y1,z1)
point2=(-x1,y1,-z1)
point3=(-x2,y2,-z2)
point4=(-x2,y2,z2)
side_ml=np.array((point4,point3,point2,point1))  # Jeff/Tahereh

# middle top
z1=600*.001
x1=(220-5-40)*.001
y1=yface
point1=(x1,y1,z1)
z2=z1
x2=-x1
y2=yface
point2=(x2,y2,z2)
z3=(260-5-40)*.001
x3=x2
y3=yface 
point3=(x3,y3,z3)
z4=z3
x4=x1
y4=yface
point4=(x4,y4,z4)
side_mt=np.array((point1,point2,point3,point4))

# mirror to middle bottom
point1=(x1,y1,-z1)
point2=(x4,y4,-z4)
point3=(x3,y3,-z3)
point4=(x2,y2,-z2)
side_mb=np.array((point1,point2,point3,point4))


myset.add_coil(side_ur)
myset.add_coil(side_mt)
myset.add_coil(side_ul)
myset.add_coil(side_mr)
myset.add_coil(side_c)
myset.add_coil(side_ml)
myset.add_coil(side_lr)
myset.add_coil(side_mb)
myset.add_coil(side_ll)


# now reflect them all to the other face: yface -> -yface
def reflect_y(points):
    newpoints=np.copy(points)
    newpoints[:,1]=-newpoints[:,1]
    newpoints=np.flip(newpoints,0) # wind them in the opposite direction
    return newpoints

oside_side_ur=reflect_y(side_ur)
oside_side_ul=reflect_y(side_ul)
oside_side_ll=reflect_y(side_ll)
oside_side_lr=reflect_y(side_lr)
oside_side_c=reflect_y(side_c)
oside_side_ml=reflect_y(side_ml)
oside_side_mr=reflect_y(side_mr)
oside_side_mt=reflect_y(side_mt)
oside_side_mb=reflect_y(side_mb)


myset.add_coil(oside_side_ul)
myset.add_coil(oside_side_mt)
myset.add_coil(oside_side_ur)
myset.add_coil(oside_side_ml)
myset.add_coil(oside_side_c)
myset.add_coil(oside_side_mr)
myset.add_coil(oside_side_ll)
myset.add_coil(oside_side_mb)
myset.add_coil(oside_side_lr)


# Double phew, now on to the top side  (Floor and ceiling)

x1=(300)*.001
y1=600*.001
z1=zface
point1=(x1,y1,z1)
x2=600*.001
y2=y1
z2=zface
point2=(x2,y2,z2)
x3=x2
y3=x1
z3=zface
point3=(x3,y3,z3)
x4=x1
y4=y3
z4=zface
point4=(x4,y4,z4)
top_ne=(point1,point2,point3,point4)
top_ne=np.array(top_ne)

# now reflect around
point1=(-x1,y1,z1)
point2=(-x4,y4,z4)
point3=(-x3,y3,z3)
point4=(-x2,y2,z2)
top_nw=np.array((point1,point2,point3,point4))

point1=(-x1,-y1,z1)
point2=(-x2,-y2,z2)
point3=(-x3,-y3,z3)
point4=(-x4,-y4,z4)
top_sw=np.array((point1,point2,point3,point4))

point1=(x1,-y1,z1)
point2=(x4,-y4,z4)
point3=(x3,-y3,z3)
point4=(x2,-y2,z2)
top_se=np.array((point1,point2,point3,point4))

# central coil
z1=zface
y1=(245+5+10)*.001
x1=(245+5+10)*.001
point1=(x1,y1,z1)
point2=(x1,-y1,z1)
point3=(-x1,-y1,z1)
point4=(-x1,y1,z1)
top_c=np.array((point1,point2,point3,point4))

# middle east coil
x1=300*.001
y1=(245+5+10)*.001
z1=zface
point1=(x1,y1,z1)
x2=600*.001
y2=y1
z2=zface
point2=(x2,y2,z2)
point3=(x2,-y2,z2)
point4=(x1,-y1,z1)
top_me=np.array((point1,point2,point3,point4))

# reflect it to middle left coil
point1=(-x1,y1,z1)
point2=(-x1,-y1,z1)
point3=(-x2,-y2,z2)
point4=(-x2,y2,z2)
top_mw=np.array((point1,point2,point3,point4))

# middle north
x1=(300-40)*.001
y1=600*.001
z1=zface
point1=(x1,y1,z1)
x2=-x1
y2=y1
z2=zface
point2=(x2,y2,z2)
x3=x2
y3=(245+5+10+40)*.001
z3=zface
point3=(x3,y3,z3)
x4=x1
y4=y3
z4=zface
point4=(x4,y4,z4)
top_mn=np.array((point4,point3,point2,point1))  # Jeff/Tahereh

# mirror to middle bottom
point1=(x1,-y1,z1)
point2=(x4,-y4,z4)
point3=(x3,-y3,z3)
point4=(x2,-y2,z2)
top_ms=np.array((point4,point3,point2,point1)) # Jeff/Tahereh


# Jeff and Tahereh completely changed the notation to indicate as-built geometry

myset.add_coil(top_nw)
myset.add_coil(top_mn)
myset.add_coil(top_ne)
myset.add_coil(top_mw)
myset.add_coil(top_c)
myset.add_coil(top_me)
myset.add_coil(top_sw)
myset.add_coil(top_ms)
myset.add_coil(top_se)



# now reflect them all to the other face: zface -> -zface
def reflect_z(points):
    newpoints=np.copy(points)
    newpoints[:,2]=-newpoints[:,2]
    newpoints=np.flip(newpoints,0) # wind them in the opposite direction
    return newpoints

bott_nw=reflect_z(top_nw)
bott_mn=reflect_z(top_mn)
bott_ne=reflect_z(top_ne)
bott_mw=reflect_z(top_mw)
bott_c=reflect_z(top_c)
bott_me=reflect_z(top_me)
bott_sw=reflect_z(top_sw)
bott_ms=reflect_z(top_ms)
bott_se=reflect_z(top_se)


myset.add_coil(bott_se)
myset.add_coil(bott_ms)
myset.add_coil(bott_sw)
myset.add_coil(bott_me)
myset.add_coil(bott_c)
myset.add_coil(bott_mw)
myset.add_coil(bott_ne)
myset.add_coil(bott_mn)
myset.add_coil(bott_nw)


# output geometry to a file
myset.output_coils('coil_points/coil_')


class sensor:
    def __init__(self,pos):
        self.pos = pos

class sensorarray:
    def __init__(self,xdim,ydim,zdim,corners):
        x = corners[1]-corners[0]
        y = corners[2]-corners[0]
        z = corners[3]-corners[0]
        #self.sensorgrid=np.mgrid[-a:a:xdim*1j,-a:a:ydim*1j,-a:a:zdim*1j]
        #print(self.sensorgrid)
        self.sensors = []
        if(xdim==1 and ydim==1 and zdim==1):
            pos = corners[0]+x/2+y/2
            self.sensors.append(sensor(pos))
            pos = corners[0]+x/2+y/2+z
            self.sensors.append(sensor(pos))
            pos = corners[0]+y/2+z/2
            self.sensors.append(sensor(pos))
            pos = corners[0]+y/2+z/2+x
            self.sensors.append(sensor(pos))
            pos = corners[0]+x/2+z/2
            self.sensors.append(sensor(pos))
            pos = corners[0]+x/2+z/2+y
            self.sensors.append(sensor(pos))
        else:
            for i in range(xdim):
                for j in range(ydim):
                    for k in range(zdim):
                        pos = corners[0]+x*i/(xdim-1)+y*j/(ydim-1)+z*k/(zdim-1)
                        self.sensors.append(sensor(pos))
        self.numsensors = len(self.sensors)
    def draw_sensor(self,number,ax):
        x = self.sensors[number].pos[0]
        y = self.sensors[number].pos[1]
        z = self.sensors[number].pos[2]
        c = 'r'
        m = 'o'
        ax.scatter(x,y,z,c=c,marker=m)
    def draw_sensors(self,ax):
        for number in range(self.numsensors):
            self.draw_sensor(number,ax)
    def vec_b(self):
        # makes a vector of magnetic fields in the same ordering as
        # the_matrix class below
        the_vector=np.zeros((self.numsensors*3))
        for j in range(myarray.numsensors):
            r = myarray.sensors[j].pos
            b=np.array([bxtarget(r[0],r[1],r[2]),
                        bytarget(r[0],r[1],r[2]),
                        bztarget(r[0],r[1],r[2])])
            for k in range(3):
                the_vector[j*3+k]=b[k]
        return the_vector


# set up array of sensors
a_sensors=0.5
p0=np.array([-a_sensors/2,-a_sensors/2,-a_sensors/2])
p1=p0+np.array([a_sensors,0,0])
p2=p0+np.array([0,a_sensors,0])
p3=p0+np.array([0,0,a_sensors])
points=(p0,p1,p2,p3)

nsensors=int(options.nsensors)
myarray=sensorarray(nsensors,nsensors,nsensors,points)
print(myarray.sensors[0].pos)
print(myarray.numsensors)
print(myarray.sensors[myarray.numsensors-1].pos)
print(myarray.sensors[myarray.numsensors-2].pos)

print('the vector test')
print(len(myarray.vec_b()),myarray.vec_b())

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['legend.fontsize'] = 10

if(options.traces):
    fig = plt.figure()
    ax=fig.add_subplot(111,projection='3d')
    myset.draw_coils(ax)
    myarray.draw_sensors(ax)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')
    plt.show()
    
#####################################################################

#drawing the coroplast and msr walls

if(options.makeplots):   #drawing the coroplast and msr walls
    # Draw wires
    fig1 = plt.figure(figsize=(9.5,9.5))
    ax4 = myset.draw_layout(fig1,title_add = " - OneWall",arrow=False,poslabel=False,drawL5=True,drawL6=True)
    
    #draw corroplast
    DrawCoruplast(ax4[0],ax4[1],ax4[2])
    
    # initialize pipes
    backpipesDraw = sparsepipelist()
    wallpipesDraw = pipelist()

    QuickFeedThroughs(wallpipesDraw,backpipesDraw,color="grey") #plotting unrerouted holes in new color.
    
    # draw pipes
    text=False
    for piplist in [backpipesDraw,wallpipesDraw]:
        piplist.draw_yz(ax4[1],text=text,div_rad=51)
        piplist.draw_xy(ax4[2],text=text,div_rad=51)
        piplist.draw_xz(ax4[0],text=text,div_rad=51)
    plt.savefig("/Users/modestekatotoka/Desktop/tucan_2024/tucan/squares_with_pipes/msr_walls/msr_walls_figures.png",dpi=300,bbox_inches='tight')
    plt.show()



from matplotlib import cm

class the_matrix:
    def __init__(self,myset,myarray):
        self.m=np.zeros((myset.numcoils,myarray.numsensors*3))
        #self.fill(myset,myarray)
        self.fillspeed(myset,myarray)
        self.condition = np.linalg.cond(self.m)

        # for some reason I chose to create the transpose of the usual
        # convention, when I first wrote the fill method
        self.capital_M=self.m.T # M=s*c=sensors*coils Matrix

        # Do the svd
        self.U,self.s,self.VT=np.linalg.svd(self.capital_M)

        print('s is',self.s)
        # s is just a list of the diagonal elements, rather than a matrix
        # You can make the matrix this way:
        self.S=np.zeros(self.capital_M.shape)
        self.S[:self.capital_M.shape[1],:self.capital_M.shape[1]]=np.diag(self.s)
        # Or, I've seen people use "full_matrices=True" in the svd command

        # Start to calculate the inverse explicitly
        # list of reciprocals
        d=1./self.s
        self.D=np.zeros(self.capital_M.shape)
        # matrix of reciprocals
        self.D[:self.capital_M.shape[1],:self.capital_M.shape[1]]=np.diag(d)

        # inverse of capital_M
        self.Minv=self.VT.T.dot(self.D.T).dot(self.U.T)
        #self.Minv=np.linalg.pinv(self.capital_M)
        
        # now gets to fixin'
        # remove just the last mode
        n_elements=myset.numcoils-1

        self.Dp=self.D[:,:n_elements]
        self.VTp=self.VT[:n_elements,:]
        self.Minvp=self.VTp.T.dot(self.Dp.T).dot(self.U.T)
        
    def fill(self,myset,myarray):
        for i in range(myset.numcoils):
            myset.set_independent_current(i,1.0)
            for j in range(myarray.numsensors):
                r = myarray.sensors[j].pos
                b = myset.b(r)
                for k in range(3):
                    self.m[i,j*3+k]=b[k]
            myset.set_independent_current(i,0.0)

    def fillspeed(self,myset,myarray):
        myset.set_common_current(1.0)
        for i in range(myset.numcoils):
            print("Doing coil %d"%i)
            for j in range(myarray.numsensors):
                r = myarray.sensors[j].pos
                bx,by,bz=myset.coil[i].b_prime(r[0],r[1],r[2])
                b=[bx,by,bz]
                for k in range(3):
                    self.m[i,j*3+k]=b[k]
        myset.zero_currents()
            
    def check_field_graphically(self,myset,myarray):
        # test each coil by graphing field at each sensor
        for i in range(myset.numcoils):
            fig = plt.figure()
            ax=fig.add_subplot(111,projection='3d')
            myset.draw_coil(i,ax)
            myset.coil[i].set_current(1.0)
            for j in range(myarray.numsensors):
                r = myarray.sensors[j].pos
                b=myset.b(r)
                bhat=b*5.e4
                points = []
                points.append(r)
                points.append(r+bhat)
                xs = ([p[0] for p in points])
                ys = ([p[1] for p in points])
                zs = ([p[2] for p in points])
                ax.plot(xs,ys,zs)
            myset.coil[i].set_current(0.0)
            ax.legend()
            plt.show()

    def show_matrices(self):
        fig1,ax1=plt.subplots()
        fig2,ax2=plt.subplots()
        fig3,ax3=plt.subplots()
        fig4,ax4=plt.subplots()
        fig5,ax5=plt.subplots()
        fig6,ax6=plt.subplots()
        fig7,ax7=plt.subplots()
        fig8,ax8=plt.subplots()
        fig9,ax9=plt.subplots()

        ax1.imshow(self.capital_M,cmap=cm.bwr)
        ax2.imshow(self.U,cmap=cm.bwr)
        ax3.imshow(self.S,cmap=cm.bwr)
        ax4.imshow(self.VT,cmap=cm.bwr)
        ax5.imshow(self.D,cmap=cm.bwr)
        ax6.imshow(self.Minv,cmap=cm.bwr)
        ax7.imshow(self.Dp,cmap=cm.bwr)
        ax8.imshow(self.VTp,cmap=cm.bwr)
        ax9.imshow(self.Minvp,cmap=cm.bwr)

        ax1.set_title('M')
        ax2.set_title('U')
        ax3.set_title('S')
        ax4.set_title('VT')
        ax5.set_title('D')
        ax6.set_title('Minv')
        ax7.set_title('Dp')
        ax8.set_title('VTp')
        ax9.set_title('Minvp')

        plt.show()

        
mymatrix=the_matrix(myset,myarray)

print('The condition number is %f'%mymatrix.condition)
if(options.matrices):
    mymatrix.show_matrices()

# Set up vector of desired fields

vec_i=mymatrix.Minv.dot(myarray.vec_b())


# Assign currents to coilcube

myset.set_currents(vec_i)

# Check the field at the center of the coilcube
r=np.array([0,0,0])
print(myset.b(r))
print(myset.b_prime(0,0,0))

from scipy.optimize import curve_fit

def fiteven(x,p0,p2,p4,p6):
    return p0+p2*x**2+p4*x**4+p6*x**6

def fitodd(x,p1,p3,p5,p7):
    return p1*x+p3*x**3+p5*x**5+p7*x**7

def fitgraph(xdata,ydata,ax):
    popt,pcov=curve_fit(fiteven,xdata[abs(xdata)<.5],ydata[abs(xdata)<.5])
    print(popt)
    ax.plot(points1d,fiteven(xdata,*popt),'r--',label='$p_0$=%2.1e,$p_2$=%2.1e,$p_4$=%2.1e,$p_6$=%2.1e'%tuple(popt))

# scans along each axis
points1d=np.mgrid[-1:1:101j]
bx1d_xscan,by1d_xscan,bz1d_xscan=myset.b_prime(points1d,0.,0.)
bx1d_yscan,by1d_yscan,bz1d_yscan=myset.b_prime(0.,points1d,0.)
bx1d_zscan,by1d_zscan,bz1d_zscan=myset.b_prime(0.,0.,points1d)

# target field
bx1d_target_xscan=bxtarget(points1d,0.,0.)*np.ones(np.shape(points1d))
bx1d_target_yscan=bxtarget(0.,points1d,0.)*np.ones(np.shape(points1d))
bx1d_target_zscan=bxtarget(0.,0.,points1d)*np.ones(np.shape(points1d))

by1d_target_xscan=bytarget(points1d,0.,0.)*np.ones(np.shape(points1d))
by1d_target_yscan=bytarget(0.,points1d,0.)*np.ones(np.shape(points1d))
by1d_target_zscan=bytarget(0.,0.,points1d)*np.ones(np.shape(points1d))

bz1d_target_xscan=bztarget(points1d,0.,0.)*np.ones(np.shape(points1d))
bz1d_target_yscan=bztarget(0.,points1d,0.)*np.ones(np.shape(points1d))
bz1d_target_zscan=bztarget(0.,0.,points1d)*np.ones(np.shape(points1d))

if(options.zoom):
    mask=(points1d>=-a_sensors/2)&(points1d<=a_sensors/2)
else:
    mask=np.full(np.shape(points1d),True)

if(options.axes):
    fig7,(ax71)=plt.subplots(nrows=1)
    fig8,(ax81)=plt.subplots(nrows=1)
    fig9,(ax91)=plt.subplots(nrows=1)
    
    ax71.plot(points1d[mask],bz1d_xscan[mask],label='$B_z(x,0,0)$')
    ax71.plot(points1d[mask],bz1d_target_xscan[mask],label='target $B_z(x,0,0)$')
    ax71.plot(points1d[mask],bz1d_yscan[mask],label='$B_z(0,y,0)$')
    ax71.plot(points1d[mask],bz1d_target_yscan[mask],label='target $B_z(0,y,0)$')
    ax71.plot(points1d[mask],bz1d_zscan[mask],label='$B_z(0,0,z)$')
    ax71.plot(points1d[mask],bz1d_target_zscan[mask],label='target $B_z(0,0,z)$')
    ax71.set_xlabel('x, y, or z (m)')
    from sympy import latex
    if(options.dipole):
        ax71.set_ylabel('$B_z=dipole$')
    else:
        ax71.set_ylabel('$B_z=\Pi_{z,%d,%d}=%s$'%(l,m,latex(sp.Piz)))
    if(not options.zoom):
        ax71.axvline(x=a/2,color='black',linestyle='--')
        ax71.axvline(x=-a/2,color='black',linestyle='--')
        ax71.axvline(x=a_sensors/2,color='red',linestyle='--')
        ax71.axvline(x=-a_sensors/2,color='red',linestyle='--')

    ax81.plot(points1d[mask],by1d_xscan[mask],label='$B_y(x,0,0)$')
    ax81.plot(points1d[mask],by1d_target_xscan[mask],label='target $B_y(x,0,0)$')
    ax81.plot(points1d[mask],by1d_yscan[mask],label='$B_y(0,y,0)$')
    ax81.plot(points1d[mask],by1d_target_yscan[mask],label='target $B_y(0,y,0)$')
    ax81.plot(points1d[mask],by1d_zscan[mask],label='$B_y(0,0,z)$')
    ax81.plot(points1d[mask],by1d_target_zscan[mask],label='target $B_y(0,0,z)$')
    ax81.set_xlabel('x, y, or z (m)')
    if(options.dipole):
        ax81.set_ylabel('$B_y=dipole$')
    else:
        ax81.set_ylabel('$B_y=\Pi_{y,%d,%d}=%s$'%(l,m,latex(sp.Piy)))
    if(not options.zoom):
        ax81.axvline(x=a/2,color='black',linestyle='--')
        ax81.axvline(x=-a/2,color='black',linestyle='--')
        ax81.axvline(x=a_sensors/2,color='red',linestyle='--')
        ax81.axvline(x=-a_sensors/2,color='red',linestyle='--')

    ax91.plot(points1d[mask],bx1d_xscan[mask],label='$B_x(x,0,0)$')
    ax91.plot(points1d[mask],bx1d_target_xscan[mask],label='target $B_x(x,0,0)$')
    ax91.plot(points1d[mask],bx1d_yscan[mask],label='$B_x(0,y,0)$')
    ax91.plot(points1d[mask],bx1d_target_yscan[mask],label='target $B_x(0,y,0)$')
    ax91.plot(points1d[mask],bx1d_zscan[mask],label='$B_x(0,0,z)$')
    ax91.plot(points1d[mask],bx1d_target_zscan[mask],label='target $B_x(0,0,z)$')
    ax91.set_xlabel('x, y, or z (m)')
    if(options.dipole):
        ax91.set_ylabel('$B_x=dipole$')
    else:
        ax91.set_ylabel('$B_x=\Pi_{x,%d,%d}=%s$'%(l,m,latex(sp.Pix)))
    if(not options.zoom):
        ax91.axvline(x=a/2,color='black',linestyle='--')
        ax91.axvline(x=-a/2,color='black',linestyle='--')
        ax91.axvline(x=a_sensors/2,color='red',linestyle='--')
        ax91.axvline(x=-a_sensors/2,color='red',linestyle='--')

    ax71.axhline(y=0,color='black')
    ax81.axhline(y=0,color='black')
    ax91.axhline(y=0,color='black')
    
    ax71.legend()
    ax81.legend()
    ax91.legend()


if(options.residuals):

    ax101=plt.figure(101)
    plt.plot(points1d[mask],bz1d_xscan[mask]-bz1d_target_xscan[mask],label='residual $B_z(x,0,0)$')
    plt.plot(points1d[mask],bz1d_yscan[mask]-bz1d_target_yscan[mask],label='residual $B_z(0,y,0)$')
    plt.plot(points1d[mask],bz1d_zscan[mask]-bz1d_target_zscan[mask],label='residual $B_z(0,0,z)$')
    plt.xlabel('x, y, or z (m)')
    plt.ylabel('residual $B_z$ (true-target)')
    plt.legend()
    if(not options.zoom):
        plt.axvline(x=a/2,color='black',linestyle='--')
        plt.axvline(x=-a/2,color='black',linestyle='--')
        plt.axvline(x=a_sensors/2,color='red',linestyle='--')
        plt.axvline(x=-a_sensors/2,color='red',linestyle='--')

    ax102=plt.figure(102)
    plt.plot(points1d[mask],by1d_xscan[mask]-by1d_target_xscan[mask],label='residual $B_y(x,0,0)$')
    plt.plot(points1d[mask],by1d_yscan[mask]-by1d_target_yscan[mask],label='residual $B_y(0,y,0)$')
    plt.plot(points1d[mask],by1d_zscan[mask]-by1d_target_zscan[mask],label='residual $B_y(0,0,z)$')
    plt.xlabel('x, y, or z (m)')
    plt.ylabel('residual $B_y$ (true-target)')
    plt.legend()
    if(not options.zoom):
        plt.axvline(x=a/2,color='black',linestyle='--')
        plt.axvline(x=-a/2,color='black',linestyle='--')
        plt.axvline(x=a_sensors/2,color='red',linestyle='--')
        plt.axvline(x=-a_sensors/2,color='red',linestyle='--')

    ax103=plt.figure(103)
    plt.plot(points1d[mask],bx1d_xscan[mask]-bx1d_target_xscan[mask],label='residual $B_x(x,0,0)$')
    plt.plot(points1d[mask],bx1d_yscan[mask]-bx1d_target_yscan[mask],label='residual $B_x(0,y,0)$')
    plt.plot(points1d[mask],bx1d_zscan[mask]-bx1d_target_zscan[mask],label='residual $B_x(0,0,z)$')
    plt.xlabel('x, y, or z (m)')
    plt.ylabel('residual $B_x$ (true-target)')
    plt.legend()
    if(not options.zoom):
        plt.axvline(x=a/2,color='black',linestyle='--')
        plt.axvline(x=-a/2,color='black',linestyle='--')
        plt.axvline(x=a_sensors/2,color='red',linestyle='--')
        plt.axvline(x=-a_sensors/2,color='red',linestyle='--')

plt.show()

print('vec_i is:',vec_i)

#generating a current.csv file

max_unnormalized_current=np.max(np.abs(vec_i)) # arb. units
max_normalized_current=0.04 # Amperes
calibration_factor=max_normalized_current/max_unnormalized_current
calibrated_vec_i=vec_i*calibration_factor # Amperes

channel_number=np.arange(54)
#channel_number =np.arange(50)
my_calibrated_array_i=calibrated_vec_i.reshape(-1,1) # Amperes
print('my calibrated currents array',my_calibrated_array_i)
print('my calibrated currents array',size(my_calibrated_array_i))

import csv

with open('current.csv','w', newline='') as csvfile:
        fields=['channel_number','calibrated_vec_i']
        writer=csv.DictWriter(csvfile, delimiter=',', lineterminator='\n',fieldnames=fields)
        writer.writeheader()
        data=(channel_number,my_calibrated_array_i)
        for cn,ci in zip(channel_number, my_calibrated_array_i):
            writer.writerow({'channel_number':cn, 'calibrated_vec_i':ci[0]})

# Now let's check what the field should be after setting these currents

myset.set_currents(calibrated_vec_i)

# the field at the centre of the coilcube
r=np.array([0,0,0])
print('Check field at centre of the coilcube')
print(myset.b(r))
print(myset.b_prime(0,0,0))


# scans along each axis
points1d=np.mgrid[-a:a:101j]
bx1d_xscan,by1d_xscan,bz1d_xscan=myset.b_prime(points1d,0.,0.)
bx1d_yscan,by1d_yscan,bz1d_yscan=myset.b_prime(0.,points1d,0.)
bx1d_zscan,by1d_zscan,bz1d_zscan=myset.b_prime(0.,0.,points1d)

# target field
bx1d_target_xscan=bxtarget(points1d,0.,0.)*np.ones(np.shape(points1d))*calibration_factor
bx1d_target_yscan=bxtarget(0.,points1d,0.)*np.ones(np.shape(points1d))*calibration_factor
bx1d_target_zscan=bxtarget(0.,0.,points1d)*np.ones(np.shape(points1d))*calibration_factor

by1d_target_xscan=bytarget(points1d,0.,0.)*np.ones(np.shape(points1d))*calibration_factor
by1d_target_yscan=bytarget(0.,points1d,0.)*np.ones(np.shape(points1d))*calibration_factor
by1d_target_zscan=bytarget(0.,0.,points1d)*np.ones(np.shape(points1d))*calibration_factor

bz1d_target_xscan=bztarget(points1d,0.,0.)*np.ones(np.shape(points1d))*calibration_factor
bz1d_target_yscan=bztarget(0.,points1d,0.)*np.ones(np.shape(points1d))*calibration_factor
bz1d_target_zscan=bztarget(0.,0.,points1d)*np.ones(np.shape(points1d))*calibration_factor


if(options.zoom):
    mask=(points1d>=-a_sensors/2)&(points1d<=a_sensors/2)
else:
    mask=np.full(np.shape(points1d),True)

if(options.axes and not options.wiggle):
    fig7,(ax71)=plt.subplots(nrows=1)
    fig8,(ax81)=plt.subplots(nrows=1)
    fig9,(ax91)=plt.subplots(nrows=1)
    #xscan plots
    ax71.plot(points1d[mask],bz1d_xscan[mask],label='$B_z(x,0,0)$')
    ax71.plot(points1d[mask],bz1d_target_xscan[mask],label='target $B_z(x,0,0)$')
    ax71.plot(points1d[mask],bz1d_yscan[mask],label='$B_z(0,y,0)$')
    ax71.plot(points1d[mask],bz1d_target_yscan[mask],label='target $B_z(0,y,0)$')
    ax71.plot(points1d[mask],bz1d_zscan[mask],label='$B_z(0,0,z)$')
    ax71.plot(points1d[mask],bz1d_target_zscan[mask],label='target $B_z(0,0,z)$')
    ax71.set_xlabel('x, y, or z (m)')
    from sympy import latex
    if(options.dipole):
        ax71.set_ylabel('$B_z=dipole$')
    else:
        ax71.set_ylabel('$B_z=\Pi_{z,%d,%d}=%s$'%(l,m,latex(sp.Piz)))
    if(not options.zoom):
        ax71.axvline(x=a/2,color='black',linestyle='--')
        ax71.axvline(x=-a/2,color='black',linestyle='--')
        ax71.axvline(x=a_sensors/2,color='red',linestyle='--')
        ax71.axvline(x=-a_sensors/2,color='red',linestyle='--')
    #yscan plots
    ax81.plot(points1d[mask],by1d_xscan[mask],label='$B_y(x,0,0)$')
    ax81.plot(points1d[mask],by1d_target_xscan[mask],label='target $B_y(x,0,0)$')
    ax81.plot(points1d[mask],by1d_yscan[mask],label='$B_y(0,y,0)$')
    ax81.plot(points1d[mask],by1d_target_yscan[mask],label='target $B_y(0,y,0)$')
    ax81.plot(points1d[mask],by1d_zscan[mask],label='$B_y(0,0,z)$')
    ax81.plot(points1d[mask],by1d_target_zscan[mask],label='target $B_y(0,0,z)$')
    ax81.set_xlabel('x, y, or z (m)')
    if(options.dipole):
        ax81.set_ylabel('$B_y=dipole$')
    else:
        ax81.set_ylabel('$B_y=\Pi_{y,%d,%d}=%s$'%(l,m,latex(sp.Piy)))
    if(not options.zoom):
        ax81.axvline(x=a/2,color='black',linestyle='--')
        ax81.axvline(x=-a/2,color='black',linestyle='--')
        ax81.axvline(x=a_sensors/2,color='red',linestyle='--')
        ax81.axvline(x=-a_sensors/2,color='red',linestyle='--')
    #zscan plots
    ax91.plot(points1d[mask],bx1d_xscan[mask],label='$B_x(x,0,0)$')
    ax91.plot(points1d[mask],bx1d_target_xscan[mask],label='target $B_x(x,0,0)$')
    ax91.plot(points1d[mask],bx1d_yscan[mask],label='$B_x(0,y,0)$')
    ax91.plot(points1d[mask],bx1d_target_yscan[mask],label='target $B_x(0,y,0)$')
    ax91.plot(points1d[mask],bx1d_zscan[mask],label='$B_x(0,0,z)$')
    ax91.plot(points1d[mask],bx1d_target_zscan[mask],label='target $B_x(0,0,z)$')
    ax91.set_xlabel('x, y, or z (m)')
    if(options.dipole):
        ax91.set_ylabel('$B_x=dipole$')
    else:
        ax91.set_ylabel('$B_x=\Pi_{x,%d,%d}=%s$'%(l,m,latex(sp.Pix)))
    if(not options.zoom):
        ax91.axvline(x=a/2,color='black',linestyle='--')
        ax91.axvline(x=-a/2,color='black',linestyle='--')
        ax91.axvline(x=a_sensors/2,color='red',linestyle='--')
        ax91.axvline(x=-a_sensors/2,color='red',linestyle='--')

    ax71.axhline(y=0,color='black')
    ax81.axhline(y=0,color='black')
    ax91.axhline(y=0,color='black')
    
    ax71.legend()
    ax81.legend()
    ax91.legend()

    np.savetxt('xscan_onesheet.out',np.transpose((points1d[mask],bx1d_xscan[mask],by1d_xscan[mask],bz1d_xscan[mask])))
    np.savetxt('xscan_onesheet_target.out',np.transpose((points1d[mask],bx1d_target_xscan[mask],by1d_target_xscan[mask],bz1d_target_xscan[mask])))

    np.savetxt('yscan_onesheet.out',np.transpose((points1d[mask],bx1d_yscan[mask],by1d_yscan[mask],bz1d_yscan[mask])))
    np.savetxt('yscan_onesheet_target.out',np.transpose((points1d[mask],bx1d_target_yscan[mask],by1d_target_yscan[mask],bz1d_target_yscan[mask])))

    np.savetxt('zscan_onesheet.out',np.transpose((points1d[mask],bx1d_zscan[mask],by1d_zscan[mask],bz1d_zscan[mask])))
    np.savetxt('zscan_onesheet_target.out',np.transpose((points1d[mask],bx1d_target_zscan[mask],by1d_target_zscan[mask],bz1d_target_zscan[mask])))

plt.show()


# Jeff says:
# Below is some stuff from Modeste, I think.  It doesn't seem to do
# anything.

# ##############################################################
# #Coils Deformation studies, run with -w

# # Now let's move some coils
# myset.set_currents(calibrated_vec_i)
# #myset.coil[0].move(-0.1,0,0)
# myset.wiggle(0.1)

# if(options.traces and options.wiggle):
#     fig = plt.figure()
#     ax=fig.add_subplot(111,projection='3d')
#     myset.draw_coils(ax)
#     #myarray.draw_sensors(ax)
#     ax.set_xlabel('x (m)')
#     ax.set_ylabel('y (m)')
#     ax.set_zlabel('z (m)')
#     plt.show()

# # output metadata, so that we know what's in these data files

# from sympy import latex

# data={
#     "l":l,
#     "m":m,
#     "Pix":latex(sp.Pix),
#     "Piy":latex(sp.Piy),
#     "Piz":latex(sp.Piz)
# }

# import json
# with open('data.json', 'w') as f:
#     json.dump(data, f)

# #########################################################################################

# #graphing the data from x-,y-,z-scan.out containing the the theoretical target fields and compare with simulations.

# # load theoretical fields

# data=np.transpose(np.loadtxt('xscan_onesheet.out'))
# x_sim,bx_sim,by_sim,bz_sim=data
# bx_sim=bx_sim*1e9 # convert to nT
# by_sim=by_sim*1e9
# bz_sim=bz_sim*1e9

# data=np.transpose(np.loadtxt('yscan_onesheet.out'))
# x_sim,bx_sim,by_sim,bz_sim=data
# bx_sim=bx_sim*1e9 # convert to nT
# by_sim=by_sim*1e9
# bz_sim=bz_sim*1e9

# data=np.transpose(np.loadtxt('zscan_onesheet.out'))
# x_sim,bx_sim,by_sim,bz_sim=data
# bx_sim=bx_sim*1e9 # convert to nT
# by_sim=by_sim*1e9
# bz_sim=bz_sim*1e9

# #meta data from theoretical fields.
# with open('data.json') as json_file:
#     graphdata=json.load(json_file)
