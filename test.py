# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 15:37:33 2018

@author: Ethan P. Kamphaus
"""
#calculate the parameters for angle term
def calcangleparam(ind,atomtypes):
    global zlist
    ans=[]
    #bond angle list
    tlist=[109.47,83.5,104.51]
    i=0
    while i < len(ind):
        atom1=ind[i]
        atom2=ind[i+1]-1
        atom3=ind[i+2]-1
        type1=atomtypes[atom1]-1
        type2=atomtypes[atom2]-1
        type3=atomtypes[atom3]-1
        thet=tlist[type2]
        #need to fix this later
        ord1=1
        ord2=ord1
        #pulling terms
        rone =calcbondst([type1,type2],ord1)
        rij=rone[0]
        rtwo=calcbondst([type2,type3],ord2)
        rjk=rtwo[0]
        rik=rij**2+rjk**2+-2*rij*rjk*math.cos(thet)
        
        Kijk=(664.12)*((zlist[type1]*zlist[type3])/(rik**5))*(rij*rjk*(1-math.cos(thet)**2)-rij**2*math.cos(thet))
        C2=1/(4*math.sin(thet)**2)
        C1=-4*C2*math.cos(thet)
        C0=C2*(2*math.cos(thet)**2+1)
        ans.append([Kijk,C0,C1,C2])
        i=i+3
        
    return ans
#calculating the "equilibrium" bond distance
def calcbondst(atomtype,order):
    global zlist
    #this is the list of the bond specification 
    #0=carbon tetrahedral
    #1=Hydrogen
    rlist=[0.757,0.460,0.658]
    #this is a list for electronegativity matches bond specification
    xlist=[5.343,4.5280,8.741]
    #this is a list for effective charges matches bond specification
    zlist=[1.912,0.712,2.3]
   
    #pulling values from list
    ri=rlist[atomtype[0]-1]
    rj=rlist[atomtype[1]-1]
    xi=xlist[atomtype[0]-1]
    xj=xlist[atomtype[1]-1]
    zi=zlist[atomtype[0]-1]
    zj=zlist[atomtype[1]-1]
   
    #calculating correction for bond order
    #not sure if this changes or not
    lam=0.1332
    rbo=-lam*(ri+rj)*np.log(1)

    #calculating correction for electronegativity
    ren=((ri*rj)*((xi)**(1/2)+(xj)**(1/2))**2)/(xi*ri+xj*rj)
    #putting terms together
    rij=ri+rj+rbo+ren
    
    #calculating force constant
    kij=664.12*((zi*zj)/(rij)**3)
    val=[rij,kij]
    return val

#to determine LJ parameters
def calcnb(atomtype):
    #this is a list for LJ 12-6 parameters following same order
    #hydrogens not included here... probably hydrogen bonding duplicating carbon?
    siglist=[3.851,3.851,3.5]
    #this is the list for eps
    elist=[0.105,0.105,0.06]
    #this is list for parital charges
    clist=[0.2,0.2,-.45]
    
    val=[]
    i=0
    for item in atomtype:
        j=0
        for item2 in atomtype:
            if j>i:
                #pulling values
                sigi=siglist[atomtype[i]-1]
                sigj=siglist[atomtype[j]-1]
                ei=elist[atomtype[i]-1]
                ej=elist[atomtype[j]-1]
                ci=clist[atomtype[i]-1]
                cj=clist[atomtype[j]-1]
                
                #arthimetic mean for sig
                sigij=(1/2)*(sigi+sigj)
                #geometric mean for eps
                eij=(ei*ej)**(1/2)
                
                val.append(sigij)
                val.append(eij)
                val.append(ci)
                val.append(cj)
            j=j+1
        i=i+1
    return val
    
#loop to determine which pairs are actually bonded
#then determine the parameters to only do once a simulation 
def calcparamst(connect,atomtype):
    rep=[]
    othrep=[]
    z=0
    for flist in connect:
        i=0
        j=0
        if flist[0]!=0:
            oth=[]
            othb=[]
            for val in np.linspace(0,1,len(flist)/2): 
                bondatom=flist[i]-1
                bondatomtype=atomtype[bondatom]
                oth.append(bondatomtype)
                firstatomtype=atomtype[z]
                bondorder=flist[i+1]
                othb.append(bondorder)
                var=calcbondst([firstatomtype,bondatomtype],bondorder)
                var.append(z)
                var.append(bondatom)
                rep.append(var)
                i=i+2
            k=0
            for item in np.linspace(0,1,len(oth)/2):
                if (abs(oth[k]-oth[k+1])==1):
                    othrep.append([z,othb[k],k+1,othb[k+1],k+2])
                k=k+2
        z=z+1    
    return rep

#function for making a circle out of a point
#use for visualization of 2-d particles
def circlepoints(cent,rad):
    t1=[]
    s1=[]
    s2=[]
    t=[]
    s=[]
    ans=[]
    t1=np.arange(-rad,rad,.1)+cent[1]
    s1=(rad**2-(t1-cent[1])**2)**(1/2)+cent[0]
    s2=-(rad**2-(t1-cent[1])**2)**(1/2)+cent[0]
    t.append(t1)
    t.append(t1)
    s.append(s1)
    s.append(s2)
    ans.append(t)
    ans.append(s)
    return ans

#to form PBC need to move atom from one side of box to other
def checkpos(pos):
    global L
    posmod=[]
    for atom in pos:
        temp=[]
        #checking the x position
        if atom[0]>L:
            temp.append(atom[0]-2*L)
        elif atom[0]<-L:
            temp.append(atom[0]+2*L)
        else:
            temp.append(atom[0])
        #checking the y position
        if atom[1]>L:
            temp.append(atom[1]-2*L)
        elif atom[1]<-L:
            temp.append(atom[1]+2*L)
        else:
            temp.append(atom[1])
        posmod.append(temp)
    return posmod

#function to determine which atom is closer in the case of PBC
def closerpos(pos1,pos2):
    global L
    rcheck=[]
    xtemp=pos2[0]
    ytemp=pos2[1]
    xtempinc=xtemp+2*L
    ytempinc=ytemp+2*L
    xtempdec=xtemp-2*L
    ytempdec=ytemp-2*L
    rcheck.append([xtemp,ytemp])
    rcheck.append([xtempinc,ytemp])
    rcheck.append([xtemp,ytempinc])
    rcheck.append([xtempinc,ytempinc])
    rcheck.append([xtempdec,ytemp])
    rcheck.append([xtemp,ytempdec])
    rcheck.append([xtempdec,ytempdec])
    rtest=[]
    rmin=L
    remx=xtemp
    remy=ytemp
    for pair in rcheck:
        rtest=dist(pos1,pair)
        if rtest < rmin:
            rmin=rtest
            remx=pair[0]
            remy=pair[1]
    closeval=[remx,remy]
    return closeval

#in order to make the Evdw and Ees =0 for bonded atoms
def connect2rpair(connect):
    i=0
    convrpair=[]
    for thing in connect:
        j=0
        if thing[0]!=0:
            for val in np.linspace(0,1,len(thing)/2): 
                atom2=thing[j]-1
                z=1
                possible=0
                while z<i:
                    possible=possible+(len(connect)-z)
                    z=z+1
                convrpair.append(possible+atom2-1)
                j=j+2
        i=i+1
    return convrpair

#to determine which atoms require angles
def connect2three(connect,pos):
    i=1
    ans=[]
    realans=[]
    
    for brack in connect:
        j=0
        atom1=i
        numcon=0   
        if brack[0]!=0:
            while j<(len(brack)-2):
                atom2=brack[j]
                j=j+2
                f=j
                while f<=(len(brack)-2):
                    atom3=brack[f]
                    ans.append(atom1)
                    mem1=atom1
                    mem2=atom2
                    mem3=atom3
                    ans.append(atom2)
                    ans.append(atom3)
                    f=f+2
                    numcon=numcon+1      
            i=i+1
            rmax=0
            rmaxpos=0
            z=0
            g=[]
            if numcon==3:
                 while z <len(ans):
                     r23=dist(pos[ans[z+1]-1],pos[ans[z+2]-1])
                     if r23 > rmax:
                         rmax=r23
                         rmaxpos=[ans[z],ans[z+1],ans[z+2]]
                     z=z+3
                 z=0
                 while z<len(ans):
                    atom1pos=ans[z]
                    atom2pos=ans[z+1]
                    atom3pos=ans[z+2]
                    if atom1pos==rmaxpos[0] and atom2pos==rmaxpos[1] and atom3pos==rmaxpos[2]:
                        f=1
                    else:
                        realans.append(atom1pos)
                        realans.append(atom2pos)
                        realans.append(atom3pos)
                    z=z+3
            else:
                j=0
                flip=5
                while j < len(realans):
                    if mem1!=realans[j] or mem2!=realans[j+1] or mem3!=realans[j+2]:
                        flip=0
                    else:
                        flip=1
                    j=j+3
                if flip==0:
                    realans.append(mem1)
                    realans.append(mem2)
                    realans.append(mem3)   
                    
            
    #for "linear connections"
    i=0
    f=[]
    while i < len(connect):
        numcon2=0
        list1=connect[i]
        if list1[0]!=0:
            x=0
            while x < len(list1):
                val1=list1[x]-1
                list2=connect[val1]
                if list2[0]!=0:
                    numcon2=numcon2+1
                    y=0
                    while y <len(list2):
                        f.append(i+1)
                        mem1=i+1
                        mem2=val1+1
                        mem3=list2[y]
                        f.append(val1+1)
                        f.append(list2[y])
                        y=y+2
                        numcon2=numcon2+1
                x=x+2
        i=i+1
        z=0
        rmax=0
        rmaxpos=0
        if numcon2==3:
             while z <len(f):
                 r23=dist(pos[f[z+1]-1],pos[f[z+2]-1])
                 if r23 > rmax:
                     rmax=r23
                     rmaxpos=[ans[z],ans[z+1],ans[z+2]]
                 z=z+3
             z=0
             while z<len(f):
                atom1pos=f[z]
                atom2pos=f[z+1]
                atom3pos=f[z+2]
                if atom1pos==rmaxpos[0] and atom2pos==rmaxpos[1] and atom3pos==rmaxpos[2]:
                    ftemp=1
                else:
                    realans.append(atom1pos)
                    realans.append(atom2pos)
                    realans.append(atom3pos)
                z=z+3
        else:
            j=0
            flip=5
            while j < len(realans):
                if mem1!=realans[j] or mem2!=realans[j+1] or mem3!=realans[j+2]:
                    flip=0
                else:
                    flip=1
                j=j+3
            if flip==0:
                  realans.append(mem1)
                  realans.append(mem2)
                  realans.append(mem3)
           
             
    return realans
#2-d Distance calculation
def dist(pos1,pos2):
    x1=pos1[0]
    x2=pos2[0]
    y1=pos1[1]
    y2=pos2[1]
    dist=((x2-x1)**2+(y2-y1)**2)**(1/2)
    return dist

#general energy expression that calls the other functions
def energy(pos,dat,ljparam,convrpair,angleparam,threes):
    E=0
    Epair=paircalc(pos,dat,ljparam,convrpair)
    Etrip=threecalc(pos,angleparam,threes)
    Equat=fourcalc(pos)
    E=Epair+Etrip+Equat
    return E

#for 4 atom interactions (dihedral angles) only for molecular systems
def fourcalc(pos):
    Equat=0
    return Equat

#this will return the interatomic distances
#for a PBC with this atom in the center
def interatomdistgen(pos):
    minr=[]
    #need to generate all other possible atoms
    i=0
    for atom in pos:
       j=0 
       for oatom in pos:
           if j>i:
              newval=closerpos(atom,oatom)
              minr.append(dist(atom,newval))
           j=j+1
       i=i+1
    return minr

#update overall temporary positiion matrix
def move(delta,pos):
    postomove=randomint(len(pos)-1)
    newpos=updatepos(pos[postomove],delta)
    mod=pos.copy()
    del mod[postomove]
    mod.insert(postomove,newpos)
    return mod

#this is where the metropolis algorithim is implemented for a NVT ensemble
def nextconfig(delta,pos,dat,ljparam,convrpair,angleparam,threes):
    global count
    new=[]
    boltz=8.617e-5
    #assumming Temperature here 
    T=300
    new=move(delta,pos)
    Eold=energy(pos,dat,ljparam,convrpair,angleparam,threes)
    Enew=energy(new,dat,ljparam,convrpair,angleparam,threes)
    if Enew<Eold:
        nextpos=new
        count=count+1
    else:
        term=math.exp(-(1/(boltz*T))*(Enew-Eold))
        test=random.random()
        if test <term:
            nextpos=new
            count=count+1
        else:
            nextpos=pos
     #making sure that PBC are kept
    fnextpos=checkpos(nextpos)
    return fnextpos

#the part of the energy due to pairs of molecules
#starting with LJ only
def paircalc(pos,dat,nbparam,convrpair):
    #LJ code
    #calculating the interatomic distances
    #uses a atom center in a box approach
    rpair=interatomdistgen(pos)
    #sig is in A
    #eps is in kcal/mol
    #actual code start
    Epair=0
    #calculating VdW forces based on LJ 12-6 potential
    #calculating electrostatic interaction
    vlj=0
    Epot=0
    i=0
    z=0
    for r in rpair:
        for valch in convrpair:
            if valch==z:
                Epot=Epot
                vlj=vlj
            else:
                eps=nbparam[i+1]
                sig=nbparam[i]
                zi=nbparam[i+2]
                zj=nbparam[i+3]
                lj=eps*((sig/r)**12-(sig/r)**6)
                vlj=vlj+lj
                cint=332.0637*((zi*zj)/(r))
                Epot=Epot+cint
        i=i+4
        z=z+1
    Epair=Epair+vlj+Epot
   
    #other pair interactions
    #bond vibration for molecular system
    i=0
    Est=0
    for data in dat:
        rij=data[0]
        kij=data[1]
        atom1=data[2]
        pos1=pos[atom1]
        atom2=data[3]
        pos2=pos[atom2]
        closeval=closerpos(pos1,pos2)
        Est=(1/2)*kij*(dist(pos1,closeval)-rij)**2+Est
    
    Epair=Epair+Epot+Est
    return Epair

#code to actually visualize the circles
def plotpoints(pos,radius):
    global L
    #create plot variables
    fig, ax = plt.subplots()
    #loop over positions
    i=0
    for rad in radius:
        #call function to get points for circle
        ans=circlepoints(pos[i],rad)
        tt=ans[0]
        st=ans[1]
        ax.plot(tt, st,'b')
        i=i+1
    #label axis
    ax.set(xlabel='x', ylabel='y',
    title='pairs')
    #ax.grid()
    #make plot square
    ax.axis('square')
    plt.xlim(-L,L)
    plt.ylim(-L,L)

    fig.savefig("test.png")
    #actually show
    plt.show()
    
#generate random int upto the max
#will create from 0 ... max
def randomint(max):
    f=random.random()
    frac=f*10%1
    numfrac=max*frac
    num=round(numfrac,0)
    num=int(num)
    return num

#for 3 atom interactions (bond angles) only for molecular systems
def threecalc(pos,angleparam,connect):
    Etrip=0
    i=0
    for val in angleparam:
        atom1=connect[i+0]-1
        atom2=connect[i+1]-1
        atom3=connect[i+2]-1
        pos1=pos[atom1]
        pos2=pos[atom2]
        pos3=pos[atom3]
        p12=dist(pos1,pos2)
        p13=dist(pos1,pos3)
        p23=dist(pos2,pos3)
        #using law of cosines to find angle
        term1=(p23**2-p12**2-p13**2)/(-2*p12*p13)
        term1=round(term1,2)
        thet=((math.acos(term1))/(math.acos(0)*2))*180
        if thet > 180:
            thet=360-thet
        E=val[0]*(val[1]+val[2]*math.cos(thet)+val[3]*math.cos(2*thet))
        Etrip=Etrip+E
        i=i=3
    return Etrip

#ck new position for atom and make move
#note that move is dependent on the energy-statistics stuff
def updatepos(pos,delta):
    newpos=[]
    #the value to move can be positive or negative
    #this should generate either a negative or positive number
    posneg1=random.random()
    check=round(posneg1,0)
    if check ==0:
        delta1=delta*-1
    else:
        delta1=delta
    posneg2=random.random()
    check=round(posneg2,0)
    if check ==0:
        delta2=delta*-1
    else:
        delta2=delta
    #actually determining the move
    xmove=random.random()*delta1
    ymove=random.random()*delta
    newpos.append(pos[0]+xmove)
    newpos.append(pos[1]+ymove)
    return newpos

#importing libraries for graphing
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import random

#size of L/2 in Angstroms
L=15
count=0
pos=[[0,0],[3,0],[-2,-2],[3,5],[7,5],[-10 ,-10],[-12,-13],[9,9],[10,0],[13,10],[10,13],[-5,-3],[-3,-5]]
#which atom is bonded with which
#specify in order (no duplicate bond description)
#also have to specify the bond order as well
connect=[[2,1,3,1,],[4,1],[12,1,13,1],[5,1],[8,1],[7,1],[0],[10,1,11,1],[0],[0],[0],[0],[0]]
#order=[[1,1],[0],[0],[0],[0],[1],[0],[0],[0]]
#need atom type for atom parameter
#assuming all tetrahedral carbon atoms with index of 1
atomtype=[3,1,1,1,3,1,1,1,1,2,2,2,2]
#calculating parameters
dat=calcparamst(connect,atomtype)
threes=connect2three(connect,pos)
angleparam=calcangleparam(threes,atomtype)
nbparam=calcnb(atomtype)
convrpair=connect2rpair(connect)

radius=[.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3]
plotpoints(pos,radius)

Elist=[]
n=10000
number=np.linspace(1,n,n)
for i in number:
    pos=nextconfig(.015,pos,dat,nbparam,convrpair,angleparam,threes)
    #probably ultra inefficient in kcal/mol now
    Elist.append(energy(pos,dat,nbparam,convrpair,angleparam,threes))
plotpoints(pos,radius)

fig, ax = plt.subplots()
ax.plot(number, Elist)
fig.savefig("test.png")
#actually show
plt.show()

print(count/n)