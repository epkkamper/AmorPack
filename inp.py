# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 12:24:51 2018

@author: Ethan Kamphaus
"""
#this code is to generate the input that would go into the monte carlo code
#this is the first step of the amorophous packing calculation
#the Monte carlo should do the rest if its not too slow

#function to add solvent
def addsolvent(rlist,newcm,atomtype):
    move=rotate(rlist)
    newpos=[]
    newatomtype=[]
    i=0
    for val in move:
        cm=newcm
        valx=val[0]
        valy=val[1]
        valz=val[2]
        cmx=cm[0]
        cmy=cm[1]
        cmz=cm[2]
        newposx=cmx+valx
        newposy=cmy+valy
        newposz=cmz+valz
        newpos.append([newposx,newposy,newposz])
        newatomtype.append(atomtype[i])
        i=i+1
        
    return newpos,newatomtype
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
        #checking z position
        if atom[2]> L:
            temp.append(atom[2]-2*L)
        elif atom[2]<-L:
            temp.append(atom[2]+2*L)
        else:
            temp.append(atom[2])
        posmod.append(temp)
    return posmod

#calculate the center of mass of a particular molecule
def com(pos):
    i=0
    count=0
    bigsumx=0
    bigsumy=0
    bigsumz=0
    while i < len(pos):
        count=count+1
        r=pos[i]
        rx=r[0]
        ry=r[1]
        rz=r[2]
        bigsumx=bigsumx+rx
        bigsumy=bigsumy+ry
        bigsumz=bigsumz+rz
        i=i+1
    bigsum=[bigsumx*(1/count),bigsumy*(1/count),bigsumz*(1/count)]
    return bigsum
    
#to determine the sizes of atoms based on VdW radii
def size(atom):
    r=0
    #carbon, hydrogen,oxygen
    if atom==1:
        r=1.7
    elif atom==2:
        r=1.09
    elif atom==3:
        r=1.52
    return r

#to make connect matrix from position data
def makeconnect(pos,length,connect):
    i=length
    newconnect=connect
    newconnect1=[]
    j=0
    s=1
    while i < len(pos):
        if j >=length:
            j=0
            s+=1
        add=connect[j]
        if add[0]!=0:   
            newconnect1=[]
            
            z=0
            while z < len(add):
                addlist=add[z]+length*(s)
                newconnect1.append(addlist)
                newconnect1.append(add[z+1])
                z=z+2
            newconnect.append(newconnect1)
        else:
            newconnect.append([0])
        i=i+1
        j=j+1
    return newconnect
#determining the size of molecule based on largest R
def molecsize(rlist,atomtype):
    #first is finding the largest magnitude
    mag=rmag(rlist)

    maxmagx=0
    maxmagy=0
    maxmagz=0
    remitx=-1
    remity=-1
    remitz=-1
    i=0
    while i < len(rlist):
        val=rlist[i]
        magx=abs(val[0])
        magy=abs(val[1])
        magz=abs(val[2])
        
        if magx>maxmagx:
            maxmagx=magx
            remitx=i
        if magy>maxmagy:
            maxmagy=magy
            remity=i
        if magz>maxmagz:
            maxmagz=magz
            remitz=i   
        i=i+1
    atomtx=atomtype[remitx]
    atomty=atomtype[remity]
    atomtz=atomtype[remitz]
    sizex=size(atomtx)
    sizey=size(atomty)
    sizez=size(atomtz)
    #I took half of the VdW to get close to Materials studio isosurfaces
    #This will need to be checked 
    #NOT ROBUST
    valx=(maxmagx+sizex)*2
    valy=(maxmagy+sizey)*2
    valz=(maxmagz+sizez)*2
    maxmag=[valx,valy,valz]
    sa=2*valx*valy+2*valx*valz+2*valy*valz
    
    prod=1
    for val in maxmag:
        prod=prod*val
        
    return prod,maxmag

#random number for positions
def randomnum(max):
    f=random.random()
    num=max*f
    return num

#determing the position vectors based on the center of mass
#this will give transformation info
def rgen(com,pos):
    rlist=[]
    for val in pos:
        valx=val[0]
        valy=val[1]
        valz=val[2]
        rx=valx-com[0]
        ry=valy-com[1]
        rz=valz-com[2]
        rlist.append([rx,ry,rz])
    return rlist

#calculate magnitudes of r transformation matirx
def rmag(rlist):
    mag=[]
    for val in rlist:
        calc=(val[0]**2+val[1]**2+val[2]**2)**(1/2)
        mag.append(calc)
    return mag

#function to determine number of solvent molecules and place them accordingly
def solplace(pos,atomtype):
    global L
    centerofm=com(pos)
    rlist=rgen(centerofm,pos)
    sizes,maxmag=molecsize(rlist,atomtype)
    cellvol=(2*L)**3
    #density won't be exactly right
    solventnum=round(cellvol/sizes,0)
    
    newpos,newatomtype=pickpos(pos,[],atomtype,solventnum,maxmag)
    newpos=checkpos(newpos)

    return newpos,newatomtype

#to rotate cartesian coordiantes randomly
def rotate(rlist):
    theta=randomnum(math.acos(0)*2)
    phi=randomnum(math.acos(0)*4)
    psi=randomnum(math.acos(0)*4)
    newrlist=[]
    for val in rlist:
       rx=val[0]
       ry=val[1]
       rz=val[2]
       R1=rx*(math.cos(phi)*math.cos(psi)-math.sin(phi)*math.cos(theta)*math.sin(psi))+ry*(-math.cos(phi)*math.sin(psi)-math.sin(phi)*math.cos(theta)*math.cos(psi))+rz*(math.sin(phi)*math.sin(theta))
       R2=rx*(math.sin(phi)*math.cos(psi)+math.cos(phi)*math.cos(theta)*math.sin(psi))+ry*(-math.sin(phi)*math.sin(psi)+math.cos(phi)*math.cos(theta)*math.cos(psi))+rz*(-math.sin(theta)*math.cos(phi))
       R3=rx*(math.sin(theta)*math.sin(phi))+ry*(math.sin(theta)*math.cos(phi))+rz*(math.cos(theta))
       newrlist.append([R1,R2,R3])
    return newrlist

def randomsph(maxmag):
    maxv=0
    maxmagscale=0.35
    for val in maxmag:
        if val > maxv:
            maxv=val
    maxv=maxv+maxmagscale*maxv
    radius=maxv/2
    theta=randomnum(math.acos(0)*4)
    phi=randomnum(math.acos(0)*2)
    cmx=radius*math.sin(theta)*math.cos(phi)
    cmy=radius*math.sin(theta)*math.sin(phi)
    cmz=radius*math.cos(theta)  
    cm=[cmx,cmy,cmz]
    return cm

#to randomly pack inital config
def pickpos(pos,limit,atomtype,solventnum,maxmag):
    #starting position
    checkerx=randomnum(1)
    if checkerx >=0.5:
        val=1
    else:
        val=-1
    startx=val*randomnum(L)
    checkery=randomnum(1)
    if checkery >=0.5:
        val=1
    else:
        val=-1
    starty=val*randomnum(L)
    checkerz=randomnum(1)
    if checkerz >=0.5:
        val=1
    else:
        val=-1
    startz=val*randomnum(L)
    
    if len(limit)>0:
        #for available space
        #will need a [] for x y and z that contains available space 
        nothing=0
    else:
        cmlist=[]
        posfinal=[]
        atomtypefinal=[]
        #adding the first molecule center of mass in
        cmlist.append([startx,starty,startz])
        
        #time to implement random rotation of molecule
        rlist=rgen(cmlist[0],pos)
        newpos,newatomtype=addsolvent(rlist,cmlist[0],atomtype)
        newpos=checkpos(newpos)
        q=0
        for val in newpos:
            posfinal.append(val)
            atomtypefinal.append(newatomtype[q])
            q+=1
        #start placing next solvent atoms
        #make sure that next cm is within radius of 1st
        solvent=solventnum-1
        x=1
        badval=0
        bigbreak=0
        warning='NEXT'
        while x<solvent and bigbreak==0:
            breakcounter=0
            breakval=0
            totval=25000
            while breakcounter<totval and breakval==0 and bigbreak==0:
                tol=0.023
                print(breakcounter)
                print(x)
                print(badval)
                if badval!=0 and (x-1-badval)>=0:
                    clist=cmlist[x-1-badval]
                    cm=randomsph(maxmag)
                    cmlist.append([cm[0]+clist[0],cm[1]+clist[1],cm[2]+clist[2]])
                    newpos,newatomtype=addsolvent(rlist,cmlist[x],atomtype)
                    newpos=checkpos(newpos)
                    #checking to see if its too close to other molecules
                    cval=0
                    for val in newpos:
                        j=0
                        n=len(posfinal)
                        while j <n:
                            listy=posfinal[j]
                            if abs(val[0]-listy[0]) <tol or abs(val[1]-listy[1]) <tol or abs(val[2]-listy[2]) <tol:
                                cval=1
                            j=j+1
                else:
                    clist=cmlist[x-1]
                    cm=randomsph(maxmag)
                    cmlist.append([cm[0]+clist[0],cm[1]+clist[1],cm[2]+clist[2]])
                    newpos,newatomtype=addsolvent(rlist,cmlist[x],atomtype)
                    newpos=checkpos(newpos)
                    #checking to see if its too close to other molecules
                    cval=0
                    for val in newpos:
                        j=0
                        n=len(posfinal)
                        while j <n:
                            listy=posfinal[j]
                            if abs(val[0]-listy[0]) <tol or abs(val[1]-listy[1]) <tol or abs(val[2]-listy[2]) <tol:
                                cval=1
                            j=j+1
                if cval==0:
                    q=0
                    for val in newpos:
                        posfinal.append(val)
                        atomtypefinal.append(newatomtype[q])
                        q+=1
                    breakval=1
                    badval=0
                else:
                    del cmlist[x]
                    breakcounter+=1
                    if breakcounter==totval and (x-1-badval)>=0:
                        badval+=1
                        breakcounter=0
                    elif (x-1-badval<0):
                        warning='TERMINATE'
                        bigbreak=1
                        
            x=x+1
        
    return posfinal,atomtypefinal

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

#write output in format of VASP POSCAR/CONTCAR
def writetofile(pos,atomtype):
    global L
    o=open('POSCAR','w')
    title='Single DME'+'\n'
    o.write(title)
    scale=1.0
    string=str(scale)+'\n'
    o.write(string)
    dim=[[2*L,0.0,0.0],[0.0,2*L,0.0],[0.0,0.0,2*L]]
    for val in dim:
        string2=''
        for item in val:
            string2=string2+str(item)+' '
        string2=string2+'\n'
        o.write(string2)
    #time for a crappy sort mechanism

    temp=atomtype
    postemp=pos
    swapped=1
    n=len(atomtype)
    while swapped==1:
        swapped=0
        i=1
        while i < n:
            val=temp[i-1]
            posval=postemp[i-1]
            if temp[i-1] >= temp[i]:
                del temp[i-1]
                del postemp[i-1]
                postemp.insert(i,posval)
                temp.insert(i,val)
                swapped=1 
            i=i+1
        n=n-1

    listch=[]
    for item in temp:
        if len(listch)==0:
            listch.append(item)
        else:
            i=0 
            chval=0
            while i < len(listch):
                if item==listch[i]:
                    chval=1
                i=i+1
            if chval==0:
                listch.append(item)
    nextline=[]
    string=''
    for item in listch:
        if item==1:
            pval="C"
        elif item==2:
            pval="H"
        elif item==3:
            pval="O"
        string=string+pval+' '
    string=string+'\n'
    o.write(string)
    
    cotrack=[]
    for item in listch:
        i=0
        count=0
        while i < len(atomtype):
            if atomtype[i]==item:
                count=count+1
            i=i+1
        cotrack.append(count)
    string=''
    for val in cotrack:
        string=string+str(val)+' '
    string=string+'\n'
    o.write(string)
    o.write('Direct \n')
    #now time to put the coordinates in
    for val in postemp:
        string=''
        for item in val:
            string=string+str((item+L)/(2*L))+' '
        string=string+'\n'
        o.write(string)
    o.close()
        
#start of actual code
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import random

L=7.5
pos=[[8.906807481342875, 6.028197086559096, 4.977726457970697], [5.430271215178088, 4.507575739351316, 4.978192487449214], [7.84977985753009, 5.100017808462681, 4.966892364549554], [6.387078239787913, 5.582142364829742, 4.988406681076546], [9.998588819138213, 5.172396239453822, 5.0083706839329825], [4.139714917326096, 5.199166861219596, 5.073900467971685], [7.820146137055629, 4.00484340907629, 3.8670227577654206], [7.862633686569854, 4.086322540300075, 6.167705356893444], [6.513454752451159, 6.647829130514011, 3.832222023984171], [6.512555866559883, 6.606235273418772, 6.2027896046438], [11.2974192711846, 6.038327599737135, 5.021865806249617], [10.35025526900572, 4.299659588629879, 3.823862868081893], [10.291127439570296, 4.326597146486957, 6.227472212337933], [4.1130983615888, 6.406991596996446, 6.169119545550114], [3.2954513788963204, 3.47841757117105, 4.985273777003181], [4.109939648795231, 6.338618887695795, 3.9887222050164066]]
atomtype=[3,3,1,1,1,1,2,2,2,2,2,2,2,2,2,2]
connect=[[3,1,5,1],[4,1,6,1],[4,1,7,1,8,1],[9,1,10,1],[11,1,12,1,13,1],[14,1,15,1,16,1],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0]]
radius=[]
cm=com(pos)
rlist=rgen(cm,pos)

newpos,ats=solplace(pos,atomtype)
print(ats)
print(newpos)
radius=[]
for val in newpos:
    radius.append(0.3)
plotpoints(newpos,radius)
writetofile(newpos,ats)
newconect=makeconnect(newpos,len(pos),connect)