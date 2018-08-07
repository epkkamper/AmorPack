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
        rik=rij*rij+rjk*rjk+-2*rij*rjk*math.cos(thet)
        
        Kijk=(664.12)*((zlist[type1]*zlist[type3])/(rik*rik*rik*rik*rik))*(rij*rjk*(1-math.cos(thet)**2)-rij*rij*math.cos(thet))
        C2=1/(4*math.sin(thet)*math.sin(thet))
        C1=-4*C2*math.cos(thet)
        C0=C2*(2*math.cos(thet)*math.cos(thet)+1)
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
    ren=((ri*rj)*(((xi)**(1/2)+(xj)**(1/2))*((xi)**(1/2)+(xj)**(1/2))))/(xi*ri+xj*rj)
    #putting terms together
    rij=ri+rj+rbo+ren
    
    #calculating force constant
    kij=664.12*((zi*zj)/(rij*rij*rij))
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
    s1=(rad*rad-(t1-cent[1])*(t1-cent[1]))**(1/2)+cent[0]
    s2=-(rad**2-(t1-cent[1])*(t1-cent[1]))**(1/2)+cent[0]
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
        #checking z position
        if atom[2]> L:
            temp.append(atom[2]-2*L)
        elif atom[2]<-L:
            temp.append(atom[2]+2*L)
        else:
            temp.append(atom[2])
        posmod.append(temp)
    return posmod

#function to determine which atom is closer in the case of PBC
def closerpos(pos1,pos2):
    global L
    rcheck=[]
    xtemp=pos2[0]
    ytemp=pos2[1]
    ztemp=pos2[2]
    xtempinc=xtemp+2*L
    ytempinc=ytemp+2*L
    xtempdec=xtemp-2*L
    ytempdec=ytemp-2*L
    ztempinc=ztemp+2*L
    ztempdec=ztemp-2*L
    rcheck.append([xtemp,ytemp,ztemp])
    rcheck.append([xtempinc,ytemp,ztemp])
    rcheck.append([xtemp,ytempinc,ztemp])
    rcheck.append([xtemp,ytemp,ztempinc])
    rcheck.append([xtempinc,ytempinc,ztemp])
    rcheck.append([xtemp,ytempinc,ztempinc])
    rcheck.append([xtempinc,ytemp,ztempinc])
    rcheck.append([xtempinc,ytempinc,ztempinc])
    rcheck.append([xtempdec,ytemp,ztemp])
    rcheck.append([xtemp,ytempdec,ztemp])
    rcheck.append([xtemp,ytemp,ztempdec])
    rcheck.append([xtemp,ytempdec,ztempdec])
    rcheck.append([xtempdec,ytemp,ztempdec])
    rcheck.append([xtempdec,ytempdec,ztemp])
    rcheck.append([xtempdec,ytempdec,ztempdec])
    rcheck.append([xtempdec,ytempinc,ztemp])
    rcheck.append([xtempdec,ytemp,ztempinc])
    rcheck.append([xtempinc,ytempdec,ztemp])
    rcheck.append([xtemp,ytempdec,ztempinc])
    rcheck.append([xtempinc,ytemp,ztempdec])
    rcheck.append([xtemp,ytempinc,ztempdec])
    rcheck.append([xtempinc,ytempinc,ztempdec])
    rcheck.append([xtempinc,ytempdec,ztempinc])
    rcheck.append([xtempdec,ytempinc,ztempinc])
    rcheck.append([xtempdec,ytempdec,ztempinc])
    rcheck.append([xtempdec,ytempinc,ztempdec])
    rcheck.append([xtempinc,ytempdec,ztempdec])
    rtest=[]
    rmin=1000
    remx=xtemp
    remy=ytemp
    remz=ztemp
    for pair in rcheck:
        rtest=effdist(pos1,pair)
        if rtest < rmin:
            rmin=rtest
            remx=pair[0]
            remy=pair[1]
            remz=pair[2]
    closeval=[remx,remy,remz]
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
                     r23=effdist(pos[ans[z+1]-1],pos[ans[z+2]-1])
                     if r23 > rmax:
                         rmax=r23
                         rmaxpos=[ans[z],ans[z+1],ans[z+2]]
                     z=z+3
                 z=0
                 rmindin=-10
                 while z<len(ans):
                     rag=effdist(pos[ans[z+1]-1],pos[ans[z+2]-1])
                     if  (rag-rmax) > rmindin:
                         rmindin=rag-rmax
                         ca2=[ans[z],ans[z+1],ans[z+2]]
                     z=z+3
                 z=0
                 while z<len(ans):
                    atom1pos=ans[z]
                    atom2pos=ans[z+1]
                    atom3pos=ans[z+2]
                    if atom1pos==rmaxpos[0] and atom2pos==rmaxpos[1] and atom3pos==rmaxpos[2] and atom1pos==ca2[0] and atom2pos==ca2[1] and atom2pos==ca2[2]:
                        nothing=1
                    else:
                        k=0
                        if len(realans)==0:
                            realans.append(atom1pos)
                            realans.append(atom2pos)
                            realans.append(atom3pos)
                        flip2=0
                        while k<len(realans):
                            if atom1pos!=realans[k] or atom2pos!=realans[k+1] or atom3pos!=realans[k+2]:
                                flip=0
                            else:
                                flip2=1
                            k=k+3
                        if flip2==0:
                            realans.append(atom1pos)
                            realans.append(atom2pos)
                            realans.append(atom3pos)
                    z=z+3
            else:
                j=0
                flip2=0
                if len(realans)==0:
                    realans.append(mem1)
                    realans.append(mem2)
                    realans.append(mem3)
                while j < len(realans):
                    if mem1!=realans[j] or mem2!=realans[j+1] or mem3!=realans[j+2]:
                        flip=0
                    else:
                        flip2=1
                    j=j+3
                if flip2==0:
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
                        mem2=i+1
                        mem1=val1+1
                        mem3=list2[y]
                        f.append(val1+1)
                        f.append(i+1)
                        f.append(list2[y])
                        y=y+2
                        numcon2=numcon2+1
                x=x+2
        i=i+1
        z=0
        rmax=0
        rmindin=-10
        if numcon2==4:
             while z <len(f):
                 r23=effdist(pos[f[z+1]-1],pos[f[z+2]-1])
                 if r23 > rmax:
                     rmax=r23
                     rmaxpos=[ans[z],ans[z+1],ans[z+2]]
                 z=z+3
             while z<len(ans):
                 rag=effdist(pos[ans[z+1]-1],pos[ans[z+2]-1])
                 if  (rag-rmax) > rmindin:
                     rmindin=rag-rmax
                     ca2=[ans[z],ans[z+1],ans[z+2]]
                 z=z+3
             z=0
             while z<len(f):
                atom1pos=f[z]
                atom2pos=f[z+1]
                atom3pos=f[z+2]
                if atom1pos==rmaxpos[0] and atom2pos==rmaxpos[1] and atom3pos==rmaxpos[2] and atom1pos==ca2[0] and atom2pos==ca2[1] and atom2pos==ca2[2]:
                    ftemp=1
                else:
                    k=0
                    if len(realans)==0:
                        realans.append(atom1pos)
                        realans.append(atom2pos)
                        realans.append(atom3pos)
                    flip2=0
                    while k<len(realans):
                        if atom1pos!=realans[k] or atom2pos!=realans[k+1] or atom3pos!=realans[k+2]:
                            flip=0
                        else:
                            flip2=1
                        k=k+3
                    if flip2==0:
                        realans.append(atom1pos)
                        realans.append(atom2pos)
                        realans.append(atom3pos)
                z=z+3
        else:
            z=0
            while z<len(f):
                atom1pos=f[z]
                atom2pos=f[z+1]
                atom3pos=f[z+2]
                mem1=atom1pos
                mem2=atom2pos
                mem3=atom3pos
                j=0
                flip2=0
                if len(realans)==0:
                    realans.append(mem1)
                    realans.append(mem2)
                    realans.append(mem3)
                while j < len(realans):
                    if mem1!=realans[j] or mem2!=realans[j+1] or mem3!=realans[j+2]:
                        flip=0
                    else:
                        flip2=1
                    j=j+3
                if flip2==0:
                      realans.append(mem1)
                      realans.append(mem2)
                      realans.append(mem3)
                z=z+3
           
    return realans

#determine if you need torsional paramters
def connect2fourt(connect):
    ans=[]
    realans=[]
    i=0
    for clist in connect:
        if clist[0]!=0:
            j=0
            while j < len(clist):
                nextatom=clist[j]-1
                nextc=connect[nextatom]
                if nextc[0]!=0:
                    z=0
                    while z < len(nextc):
                        nextatom2=nextc[z]-1
                        nextc2=connect[nextatom2]
                        if nextc2[0]!=0:
                            h=0
                            while h < len(nextc2):
                                nextatom3=nextc2[h]-1
                                p=0
                                if len(realans)==0:
                                    realans.append(i+1)
                                    realans.append(nextatom+1)
                                    realans.append(nextatom2+1)
                                    realans.append(nextatom3+1)
                                else:
                                    val=0
                                    while p < len(realans):
                                        if (i+1)!=realans[p] or (nextatom+1)!=realans[p+1] or (nextatom2+1)!=realans[p+2] or (nextatom3+1)!=realans[p+3]:
                                            val2=0
                                        else:
                                            val=1
                                        p=p+4
                                    if val==0:
                                        realans.append(i+1)
                                        realans.append(nextatom+1)
                                        realans.append(nextatom2+1)
                                        realans.append(nextatom3+1)
                                h=h+2
                        else:
                            h=j+2
                            while h < len(clist):
                                nextatom3=clist[h]-1
                                
                                p=0
                                if len(realans)==0:
                                    realans.append(nextatom3+1)
                                    realans.append(i+1)
                                    realans.append(nextatom+1)
                                    realans.append(nextatom2+1)
                                else:
                                    val=0
                                    while p < len(realans):
                                        if (nextatom3+1)!=realans[p] or (i+1)!=realans[p+1] or (nextatom+1)!=realans[p+2] or (nextatom2+1)!=realans[p+3]:
                                            val2=0
                                        else:
                                            val=1
                                        p=p+4
                                    if val==0:
                                        realans.append(nextatom3+1)
                                        realans.append(i+1)
                                        realans.append(nextatom+1)
                                        realans.append(nextatom2+1)
                                h=h+2
                            h=j-2
                            while h >=0 :
                                nextatom3=clist[h]-1
                                
                                p=0
                                if len(realans)==0:
                                    realans.append(nextatom3+1)
                                    realans.append(i+1)
                                    realans.append(nextatom+1)
                                    realans.append(nextatom2+1)
                                else:
                                    val=0
                                    while p < len(realans):
                                        if (nextatom3+1)!=realans[p] or (i+1)!=realans[p+1] or (nextatom+1)!=realans[p+2] or (nextatom2+1)!=realans[p+3]:
                                            val2=0
                                        else:
                                            val=1
                                        p=p+4
                                    if val==0:
                                        realans.append(nextatom3+1)
                                        realans.append(i+1)
                                        realans.append(nextatom+1)
                                        realans.append(nextatom2+1)
                                h=h-2
                        z=z+2
                else:
                    z=j+2
                    while z < len(clist):
                        nextatom2=clist[z]-1
                        nextc2=connect[nextatom2]
                        
                        if nextc2[0]!=0:
                            
                            h=0
                            while h< len(nextc2):
                                nextatom3=nextc2[h]-1
                                
                                p=0
                                if len(realans)==0:
                                    realans.append(nextatom+1)
                                    realans.append(i+1)
                                    realans.append(nextatom2+1)
                                    realans.append(nextatom3+1)
                                    
                                else:
                                    val=0
                                    while p < len(realans):
                                        if (nextatom+1)!=realans[p] or (i+1)!=realans[p+1] or (nextatom2+1)!=realans[p+2] or (nextatom3+1)!=realans[p+3]:
                                            val2=0
                                        else:
                                            val=1
                                        p=p+4
                                    if val==0:
                                        realans.append(nextatom+1)
                                        realans.append(i+1)
                                        realans.append(nextatom2+1)
                                        realans.append(nextatom3+1)
                                h=h+2
                        z=z+2
                    z=j-2
                    while z >=0:
                        nextatom2=clist[z]-1
                        nextc2=connect[nextatom2]
                        if nextc2[0]!=0:
                            h=0
                            while h< len(nextc2):
                                nextatom3=nextc2[h]-1    
                                p=0
                                
                                if len(realans)==0:
                                    realans.append(nextatom+1)
                                    realans.append(i+1)
                                    realans.append(nextatom2+1)
                                    realans.append(nextatom3+1)
                                else:
                                    val=0
                                    while p < len(realans):
                                        if (nextatom+1)!=realans[p] or (i+1)!=realans[p+1] or (nextatom2+1)!=realans[p+2] or (nextatom3+1)!=realans[p+3]:
                                            val2=0
                                        else:
                                            val=1
                                        p=p+4
                                    if val==0:
                                        realans.append(nextatom+1)
                                        realans.append(i+1)
                                        realans.append(nextatom2+1)
                                        realans.append(nextatom3+1)
                                h=h+2
                        z=z-2
                j=j+2
        i=i+1
        
    i=0
    while i <len(connect):
        clist=connect[i]
        if clist[0]!=0:
            g=0
            while g < len(clist):                
                f=0
                nextatom=clist[g]-1
                nextc=connect[nextatom]
                while f < len(nextc):
                    if nextc[0]!=0:
                        valse=nextc[f]
                        x=0
                        while x < len(connect):
                            ch=connect[x]
                            if ch[0]!=0:
                                y=0
                                while y < len(ch):
                                    if ch[y]==valse and x!=i and x!=nextatom:
                                        nextatom2=x
                                        nextc2=connect[nextatom2]
                                        if nextc2[0]!=0:
                                            z=0
                                            while z <len(nextc2):
                                                nextatom3=valse-1
                                                if nextatom3==nextatom2 or nextatom3==nextatom or nextatom3==i:
                                                    nothing=0
                                                else:
                                                    p=0
                                                    if len(realans)==0:
                                                        realans.append(i+1)
                                                        realans.append(nextatom+1)
                                                        realans.append(nextatom3+1)
                                                        realans.append(nextatom2+1)
                                                    val=0
                                                    while p < len(realans):
                                                        if (i+1)!=realans[p] or (nextatom+1)!=realans[p+1] or (nextatom3+1)!=realans[p+2] or (nextatom2+1)!=realans[p+3]:
                                                            val2=0
                                                        else:
                                                            val=1
                                                        p=p+4
                                                    if val==0:
                                                        realans.append(i+1)
                                                        realans.append(nextatom+1)
                                                        realans.append(nextatom3+1)
                                                        realans.append(nextatom2+1)
                                                z=z+2
                                    y=y+2
                            x=x+1
                    f=f+2
                g=g+2
        i=i+1
    i=0
    while i <len(connect):
        clist=connect[i]
        if clist[0]!=0:
            g=0
            while g < len(clist):                
                f=0
                valse=clist[g]
                nextatom=clist[g]-1
                while f < len(connect):
                    ch=connect[f]
                    if ch[0]!=0:
                        x=0
                        while x < len(ch):
                            if ch[x]==valse and f!=i and f!=nextatom:
                                nextatom2=f
                                nextc2=connect[nextatom2]
                                if nextc2[0]!=0:
                                    z=0
                                    while z <len(nextc2):
                                        nextatom3=nextc2[z]-1
                                        if nextatom3==nextatom2 or nextatom3==nextatom or nextatom3==i:
                                            nothing=0
                                        else:
                                            p=0
                                            if len(realans)==0:
                                                realans.append(i+1)
                                                realans.append(nextatom+1)
                                                realans.append(nextatom2+1)
                                                realans.append(nextatom3+1)
                                            val=0
                                            while p < len(realans):
                                                if (i+1)!=realans[p] or (nextatom+1)!=realans[p+1] or (nextatom3+1)!=realans[p+2] or (nextatom2+1)!=realans[p+3]:
                                                    val2=0
                                                else:
                                                    val=1
                                                p=p+4
                                            if val==0:
                                                realans.append(i+1)
                                                realans.append(nextatom+1)
                                                realans.append(nextatom2+1)
                                                realans.append(nextatom3+1)
                                        z=z+2
                            x=x+2
                    f=f+1
                g=g+2
        i=i+1
        
    return realans

#for atoms involved in inversion not the same as torsional
def connect2fouri(connect):
    i=0
    j=0
    ans=[]
    for clist in connect:
        if len(clist) >= 6:
            if j < len(clist):
                atom2=clist[j]
                atom3=clist[j+2]
                atom4=clist[j+4]
            ans.append(i+1)
            ans.append(atom2)
            ans.append(atom3)
            ans.append(atom4)
        i=i+1
    
    u=0    
    for clist in connect:
        numcon=0
        if clist[0]!=0:
            j=0
            while j < len(clist):
                t=[]
                numcon=0
                nextval=clist[j]-1
                nextlist=connect[nextval]
                numcon=numcon+1
                if nextlist[0]!=0:
                    z=0
                    t.append(nextval+1)
                    t.append(u+1)
                    while z < len(nextlist):
                        numcon=numcon+1
                        t.append(nextlist[z])
                        z=z+2    
                j=j+2
                if numcon>=4:
                    atomo=[]
                    k=1
                    atom1=t[0]
                    cap=len(t)-2
                    while k < len(t):
                        atomo.append(t[k])
                        k=k+1
                    l=0
                    while l <= len(atomo)-3:
                        if l==len(atomo)-3:
                            tt=[]
                            tt.append(atom1)
                            tt.append(atomo[l])
                            z=l+1
                            while z <= len(atomo)-1:
                                tt.append(atomo[z])
                                z=z+1
                            cval2=0
                            h=0
                            while h <len(ans):
                                if tt[0]!=ans[h] or tt[1]!=ans[h+1] or tt[2]!=ans[h+2] or tt[3]!=ans[h+3]:
                                    cval=1
                                else:
                                    cval2=1
                                h=h+4
                            if cval2==0:
                                ans.append(tt[0])
                                ans.append(tt[1])
                                ans.append(tt[2])
                                ans.append(tt[3])
                        else:
                            d=l+1
                            while d <= len(atomo)-1:
                                x=l+1
                                tt=[]
                                tt.append(atom1)
                                tt.append(atomo[l])
                                while x <= len(atomo)-1:
                                    if x !=d:
                                        tt.append(atomo[x])
                                    x=x+1
                                d=d+1
                                cval2=0
                                h=0
                                while h <len(ans):
                                    if tt[0]!=ans[h] or tt[1]!=ans[h+1] or tt[2]!=ans[h+2] or tt[3]!=ans[h+3]:
                                        cval=1
                                    else:
                                        cval2=1
                                    h=h+4
                                if cval2==0:
                                    ans.append(tt[0])
                                    ans.append(tt[1])
                                    ans.append(tt[2])
                                    ans.append(tt[3])
                        l=l+1
                else:
                    k=0
                    con=0
                    track=[]
                    while k<len(connect):
                        q=0
                        test=connect[k]
                        if test[0]!=0:
                            while q < len(test):
                                scanval=test[q]
                                p=k+1
                                while p <len(connect):
                                    nextlist=connect[p]
                                    f=0
                                    if nextlist[0]!=0:
                                        while f<len(nextlist):
                                            nextval=nextlist[f]
                                            if scanval==nextval:
                                                con=con+1
                                                track.append(k+1)
                                                track.append(p+1)
                                                track.append(scanval)
                                            f=f+2
                                    p=p+1
                                q=q+2
                        k=k+1   
                    if con!=0:
                        t=[]
                        
                        t.append(track[2])
                        t.append(track[0])
                        t.append(track[1])
                        k=0
                        listy=connect[track[2]-1]
                        while k<len(listy)-1:
                            t.append(listy[k])
                            k=k+2   

                        atomo=[]
                        k=1
                        atom1=t[0]
                        cap=len(t)-2
                        while k < len(t):
                            atomo.append(t[k])
                            k=k+1
                        l=0
                        while l <= len(atomo)-3:
                            if l==len(atomo)-3:
                                tt=[]
                                tt.append(atom1)
                                tt.append(atomo[l])
                                z=l+1
                                while z <= len(atomo)-1:
                                    tt.append(atomo[z])
                                    z=z+1
                                cval2=0
                                h=0
                                while h <len(ans):
                                    if tt[0]!=ans[h] or tt[1]!=ans[h+1] or tt[2]!=ans[h+2] or tt[3]!=ans[h+3]:
                                        cval=1
                                    else:
                                        cval2=1
                                    h=h+4
                                if cval2==0:
                                    ans.append(tt[0])
                                    ans.append(tt[1])
                                    ans.append(tt[2])
                                    ans.append(tt[3])
                            else:
                                d=l+1
                                while d <= len(atomo)-1:
                                    x=l+1
                                    tt=[]
                                    tt.append(atom1)
                                    tt.append(atomo[l])
                                    while x <= len(atomo)-1:
                                        if x !=d:
                                            tt.append(atomo[x])
                                        x=x+1
                                    d=d+1
                                    cval2=0
                                    h=0
                                    while h <len(ans):
                                        if tt[0]!=ans[h] or tt[1]!=ans[h+1] or tt[2]!=ans[h+2] or tt[3]!=ans[h+3]:
                                            cval=1
                                        else:
                                            cval2=1
                                        h=h+4
                                    if cval2==0:
                                        ans.append(tt[0])
                                        ans.append(tt[1])
                                        ans.append(tt[2])
                                        ans.append(tt[3])
                            l=l+1
        u=u+1
    return ans

#function to determine dihedral
def dihedral(pos1,pos2,pos3,pos4):
    #determine formula for plane 1
    vec12=[pos1[0]-pos2[0],pos1[1]-pos2[1],pos1[2]-pos2[2]]
    vec13=[pos1[0]-pos3[0],pos1[1]-pos3[1],pos1[2]-pos3[2]]
    vecn1=[(vec12[1]*vec13[2]-vec12[2]*vec13[1]),(vec12[0]*vec13[2]-vec12[2]*vec13[0]),(vec12[0]*vec13[1]-vec12[1]*vec13[0])]
    d1=-1*(pos1[0]*vecn1[0]+pos1[1]*vecn1[1]+pos1[2]*vecn1[2])
    #determine formula for plane 2
    vec14=[pos1[0]-pos4[0],pos1[1]-pos4[1],pos1[2]-pos4[2]]
    vecn2=[(vec12[1]*vec14[2]-vec12[2]*vec14[1]),(vec12[0]*vec14[2]-vec12[2]*vec14[0]),(vec12[0]*vec14[1]-vec12[1]*vec14[0])]
    d2=-1*(pos1[0]*vecn1[0]+pos1[1]*vecn1[1]+pos1[2]*vecn1[2])
    
    #now to calculate the angle between plane 1 and 2 
    thetad=math.acos((vecn1[0]*vecn2[0]+vecn1[1]*vecn2[1]+vecn1[2]*vecn2[2])/((vecn1[0]*vecn1[0]+vecn1[1]*vecn1[1]+vecn1[2]*vecn1[2])**(1/2)*(vecn2[0]*vecn2[0]+vecn2[1]*vecn2[1]+vecn2[2]*vecn2[2])**(1/2)))
    return thetad
#2-d Distance calculation
def dist(pos1,pos2):
    x1=pos1[0]
    x2=pos2[0]
    y1=pos1[1]
    y2=pos2[1]
    z1=pos1[2]
    z2=pos2[2]
    dist=((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))**(1/2)
    return dist

#2-d Distance calculation without square root for some situations
def effdist(pos1,pos2):
    x1=pos1[0]
    x2=pos2[0]
    y1=pos1[1]
    y2=pos2[1]
    z1=pos1[2]
    z2=pos2[2]
    dist=((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))
    return dist

#general energy expression that calls the other functions
def energy(pos,dat,ljparam,convrpair,angleparam,threes,fourst,foursi):
    E=0
    Epair=paircalc(pos,dat,ljparam,convrpair)
    Etrip=threecalc(pos,angleparam,threes)
    Equat=fourcalc(pos,foursi,fourst)
    E=Epair+Etrip+Equat
    return E

#for 4 atom interactions (dihedral angles) only for molecular systems inversion
def fourcalc(pos,foursi,fourst):
    Equat=0
    #torsional calculation 1st
    i=0
    while i < len(fourst):
        posnum1=fourst[i]
        posnum2=fourst[i+1]
        posnum3=fourst[i+2]
        posnum4=fourst[i+3]
        #this switch up is to make the dihedral function work correclty 
        pos3=pos[posnum1-1]
        pos1=pos[posnum2-1]
        pos2=pos[posnum3-1]
        pos4=pos[posnum4-1]
        thetac=dihedral(pos1,pos2,pos3,pos4)
        thetac=(thetac/(2*math.acos(0)))*180
        #the Vs are very dependent on atom type and bonding
        #parameters listed here are for two sp3 bonded 
        thetas=180
        Vjk=2
        njk=3
        Equat=Equat+(1/2)*Vjk*(1-math.cos(njk*(thetac-thetas)))
        i=i+4
    #now time for inversion calc
    j=0
    while j < len(foursi):
        posnum1=foursi[j]
        posnum2=foursi[j+1]
        posnum3=foursi[j+2]
        posnum4=foursi[j+3]
        pos1=pos[posnum1-1]
        pos2=pos[posnum2-1]
        pos3=pos[posnum3-1]
        pos4=pos[posnum4-1]
        thetai=inverang(pos1,pos2,pos3,pos4)
        thetai=(thetai/(2*math.acos(0)))*180
        #parameters are basically only for carbon-carbon
        #the force constant will be significanlty different for C-O bonding
        Kijkl=6
        C1=-1
        C0=1
        C2=0
        Equat=Equat+Kijkl*(C0+C1*math.sin(thetai)+C2*math.cos(2*thetai))
        j=j+4
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
              minr.append(effdist(atom,newval))
           j=j+1
       i=i+1
    return minr

#calculate inversion angle
def inverang(pos1,pos2,pos3,pos4):
    vec1=[pos1[0]-pos2[0],pos1[1]-pos2[1],pos1[2]-pos2[2]]
    vec2=[pos2[0]-pos3[0],pos2[1]-pos3[1],pos2[2]-pos3[2]]
    vec3=[pos2[0]-pos4[0],pos2[1]-pos4[1],pos2[2]-pos4[2]]
    vecn=[(vec2[1]*vec3[2]-vec2[2]*vec3[1]),(vec2[0]*vec3[2]-vec2[2]*vec3[0]),(vec2[0]*vec3[1]-vec2[1]*vec3[0])]
    
    thetar=math.asin((vec1[0]*vecn[0]+vec1[1]*vecn[1]+vec1[2]*vecn[2])/((vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2])**(1/2)*(vecn[0]*vecn[0]+vecn[1]*vecn[1]+vecn[2]*vecn[2])**(1/2)))
    return thetar
#update overall temporary positiion matrix
def move(delta,pos):
    postomove=randomint(len(pos)-1)
    newpos=updatepos(pos[postomove],delta)
    mod=pos.copy()
    del mod[postomove]
    mod.insert(postomove,newpos)
    return mod

#this is where the metropolis algorithim is implemented for a NVT ensemble
def nextconfig(delta,pos,dat,ljparam,convrpair,angleparam,threes,fourst,foursi):
    global count
    new=[]
    repE=0
    boltz=8.617e-5
    #assumming Temperature here 
    T=300
    new=move(delta,pos)
    Eold=energy(pos,dat,ljparam,convrpair,angleparam,threes,fourst,foursi)
    Enew=energy(new,dat,ljparam,convrpair,angleparam,threes,fourst,foursi)
    if Enew<Eold:
        nextpos=new
        count=count+1
        repE=Enew
    else:
        term=math.exp(-(1/(boltz*T))*(Enew-Eold))
        test=random.random()
        if test <term:
            nextpos=new
            count=count+1
            repE=Enew
        else:
            nextpos=pos
            repE=Eold
     #making sure that PBC are kept
    fnextpos=checkpos(nextpos)
    return fnextpos,repE

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
                lj=eps*((((sig*sig*sig*sig*sig*sig)*(sig*sig*sig*sig*sig*sig))/(r*r*r*r*r*r))-((sig*sig*sig*sig*sig*sig)/(r*r*r)))
                vlj=vlj+lj
                cutoff=1.2
                if r < ((cutoff*sig)*(cutoff*sig)): 
                    cint=332.0637*((zi*zj)/(r**(1/2)))
                else:
                    cint=0
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
        p12=effdist(pos1,pos2)
        p13=effdist(pos1,pos3)
        p23=effdist(pos2,pos3)
        #using law of cosines to find angle
        term1=(p23-p12-p13)/(-2*p12**(1/2)*p13**(1/2))
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
    posneg3=random.random()
    check=round(posneg3,0)
    if check==0:
        delta3=delta*-1
    else:
        delta3=delta
    #actually determining the move
    xmove=random.random()*delta1
    ymove=random.random()*delta2
    zmove=random.random()*delta3
    newpos.append(pos[0]+xmove)
    newpos.append(pos[1]+ymove)
    newpos.append(pos[2]+zmove)
    return newpos

#write positions in a format that can be copy and pasted into gaussian
def writepos(pos,atomtype):
    o=open('outpos','w')
    i=0
    for item in atomtype:
        if item==1:
            pval="C"
        elif item==2:
            pval="H"
        elif item==3:
            pval="O"
        posv=pos[i]
        posx=posv[0]
        posy=posv[1]
        posz=posv[2]
        string=pval+' '+str(posx)+' '+str(posy)+' '+str(posz)+ '\n'
        o.write(string)
        i=i+1
    o.close()

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
        
#importing libraries for graphing
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import random

#size of L/2 in Angstroms
L=7.5
count=0
#which atom is bonded with which
connect=[[3, 1, 5, 1], [4, 1, 6, 1], [4, 1, 7, 1, 8, 1], [9, 1, 10, 1], [11, 1, 12, 1, 13, 1], [14, 1, 15, 1, 16, 1], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [19, 1, 21, 1], [20, 1, 22, 1], [20, 1, 23, 1, 24, 1], [25, 1, 26, 1], [27, 1, 28, 1, 29, 1], [30, 1, 31, 1, 32, 1], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [35, 1, 37, 1], [36, 1, 38, 1], [36, 1, 39, 1, 40, 1], [41, 1, 42, 1], [43, 1, 44, 1, 45, 1], [46, 1, 47, 1, 48, 1], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [51, 1, 53, 1], [52, 1, 54, 1], [52, 1, 55, 1, 56, 1], [57, 1, 58, 1], [59, 1, 60, 1, 61, 1], [62, 1, 63, 1, 64, 1], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [67, 1, 69, 1], [68, 1, 70, 1], [68, 1, 71, 1, 72, 1], [73, 1, 74, 1], [75, 1, 76, 1, 77, 1], [78, 1, 79, 1, 80, 1], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [83, 1, 85, 1], [84, 1, 86, 1], [84, 1, 87, 1, 88, 1], [89, 1, 90, 1], [91, 1, 92, 1, 93, 1], [94, 1, 95, 1, 96, 1], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [99, 1, 101, 1], [100, 1, 102, 1], [100, 1, 103, 1, 104, 1], [105, 1, 106, 1], [107, 1, 108, 1, 109, 1], [110, 1, 111, 1, 112, 1], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [115, 1, 117, 1], [116, 1, 118, 1], [116, 1, 119, 1, 120, 1], [121, 1, 122, 1], [123, 1, 124, 1, 125, 1], [126, 1, 127, 1, 128, 1], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [131, 1, 133, 1], [132, 1, 134, 1], [132, 1, 135, 1, 136, 1], [137, 1, 138, 1], [139, 1, 140, 1, 141, 1], [142, 1, 143, 1, 144, 1], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [147, 1, 149, 1], [148, 1, 150, 1], [148, 1, 151, 1, 152, 1], [153, 1, 154, 1], [155, 1, 156, 1, 157, 1], [158, 1, 159, 1, 160, 1], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0]]
#specify in order (no duplicate bond description)
#also have to specify the bond order as well
atomtype=[3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
#order=[[1,1],[0],[0],[0],[0],[1],[0],[0],[0]]
#need atom type for atom parameter
#assuming all tetrahedral carbon atoms with index of 1
pos=[[-0.6721305228640979, -1.1220702419685615, 4.90561573968102], [-3.286679829208907, 1.6279057318347554, 4.908579963398111], [-1.9050578361740662, -0.44482811381054466, 4.9006186375654615], [-1.9502964161135594, 1.0947602226441544, 4.910956188735048], [-1.1044073468625832, -2.4397466664071454, 4.9500039737986175], [-3.076020144568357, 3.0778595287296113, 4.991520278747659], [-2.949191139150223, -0.8021911917802222, 3.812638507651922], [-2.848796495760988, -0.7892919282631556, 6.112449672906257], [-0.9100417722679772, 1.3265423853104423, 3.743747492515883], [-0.9398545051063216, 1.3390464745591346, 6.114604522118206], [0.15243157757268655, -3.3654861299321897, 4.959319139086327], [-1.809798916869406, -3.0806754584127845, 3.776510789518996], [-1.7948989126301704, -2.9897163729848732, 6.179424784167859], [-1.9451757789129083, 3.526442441275318, 6.073394833963086], [-4.981839807704482, 3.2841328296159444, 4.918234481685998], [-2.019356512647464, 3.4823487207012587, 3.8938771731380655], [-7.027062120826987, 5.96527011887337, -3.4929384519323854], [4.598011602337139, -7.366397455773628, -3.801325623509147], [6.968127694288299, 6.872563255832203, -3.820159475449053], [5.497535602477027, 6.613713158070558, -3.3833560813094508], [-5.914823011770766, 6.553948051413824, -4.116683915208924], [3.278494765673993, 7.139903116881923, -3.3418490384492614], [7.138544846001903, -6.624116102188072, -3.460390547745831], [6.856659552158186, 7.103610030068035, -5.361306255912492], [5.738324441535445, 6.3001310100983385, -1.8720975195167533], [5.409010700013838, 5.100386157906437, -3.8899417718115217], [-4.668513270721869, 5.671432854306646, -3.857814851267861], [-5.364833168990206, -7.128755659505446, -3.622861375895205], [-5.759269488181349, 6.602662495775277, -5.62471564971213], [3.048287994108132, 5.5498200587091535, -3.63117650733048], [2.5286940823565516, -6.264720642443848, -4.046119311862977], [3.3517716478610415, 6.74401704735012, -1.8302049348489078], [-4.07242466277997, 0.3167534410279158, 0.7297629959371577], [-5.824776210347398, 3.0262577300234943, 2.338345949366726], [-4.3270206383494285, 1.4954985545907302, 1.4407339758881577], [-5.732951484718868, 1.778353332438555, 1.5874881402876966], [-2.7422141288249957, 0.5101295909140564, 0.8068164531022433], [-7.242768056047861, 3.1051427839257606, 2.281143564219043], [-3.4128499173084474, 1.8529222177700584, 2.883145982815094], [-4.0342197262636725, 2.8209916806715736, 0.8791917078954974], [-5.962061567061805, 0.3716026758096742, 2.050723039485682], [-6.562174263720795, 1.4858027018219904, 0.05494378052900473], [-2.267151777132343, -0.726023473948656, 0.05833642265577765], [-1.6464768800801224, 0.48877791978702945, 2.1193049665150596], [-2.3401896456513445, 1.5908018499432792, 0.07949759032628378], [6.734817533988155, 2.683381962414222, 0.8060902339371623], [-6.846243172155214, 4.766880924506808, 3.363873738349156], [7.346767399357773, 1.7427128473055857, 2.69146848615798], [0.6145335730361676, -2.482438511667679, 2.6291302664734957], [3.4323507069248347, -1.6944995266222804, 3.70197882588608], [1.8029637108514116, -2.5698923494657833, 2.753871854898673], [2.1922959850961163, -1.4484919476697176, 3.6941942780476893], [0.6475313994441496, -3.69388205683443, 1.6697973583021763], [3.5507096407909007, -0.5498867474130638, 4.594878644927785], [3.14055972336576, -2.839642008508489, 3.086126787714928], [1.899235995682762, -3.8393068474024146, 1.4071940849318096], [1.966843847348633, -0.23706295620811346, 4.96102368119557], [0.798244222708977, -1.3325057256747614, 3.197201786407721], [-0.6206219294786184, -3.782169072164861, 1.3928809292543818], [1.6798609611426514, -3.972341173892941, 1.9717053533896056], [0.4759507115431667, -5.0024600593913195, 0.24380981767827592], [2.164968327824763, -0.16748501041831698, 4.345149486462809], [5.2260276613526315, -1.281691001388177, 4.312864136518776], [3.316483277311118, 0.766240479399884, 5.921527192832102], [-3.443896776919104, -5.987863341013394, 5.88996576074161], [-4.611923408330096, 6.105651733250541, 5.19399883265589], [-4.261358159795131, -6.888504036658287, 5.417763935503338], [-3.639073557409782, 6.932380041336977, 5.756713706777045], [-4.403669315572186, -5.159565874624423, 5.3120154254035015], [-3.802408466289126, 5.034220126217381, 5.584635036717405], [-5.288434555298313, -6.331001913983758, 5.680862276686092], [-5.333575088474108, 7.376515206439022, 3.8532333020258385], [-2.5312958717588447, -7.242516001640879, 7.307625134091957], [-2.6981570205717684, 6.383094315259704, 5.3542076148669615], [-3.6720653509254753, -4.066727614225267, 5.733222361834301], [-5.242847812366476, -4.2315626532445645, 5.752494815794069], [-5.337613859148485, -5.669047722976552, 3.814921628878146], [-2.6575935677338194, 4.439247009702111, 5.393480467397502], [-5.4243252502739985, 4.31218046980667, 4.705903494838523], [-2.6096838838942116, 5.694806703195621, 7.128708690621021], [-4.672993057202433, 5.201833176659436, -0.39586769238803754], [-6.631283427387572, -6.7292329500584085, -1.7977957745240758], [-5.6053011731444045, 6.254938608234756, -0.5545926760081397], [-5.621230899190415, 7.254132958239646, -1.7911410366307337], [-4.981266772222931, 4.617961518527055, 0.8394630149055224], [-6.4063609994339, -5.874969492745773, -3.0084588706118396], [-7.113064259778199, 5.930951916651113, -0.45186544131040307], [-5.527618693799283, 7.223202286185174, 0.613770638586967], [-5.597351119479978, 6.191516762953057, -2.87708313504028], [-4.063933809659503, -7.408773306299402, -1.7350402222683705], [-4.029874431270073, 3.3985796234186645, 1.1953207715885465], [-6.290634350523266, 3.8927560492552815, 1.0234654974064146], [-4.701713941063019, 5.337693534928425, 2.1025826608945497], [-4.836224129011271, -5.541306264203222, -3.2164153548681043], [7.103169330336698, -4.814837957737748, -2.5996195738096084], [-6.323098818122752, -6.796744482030347, -4.207048606585417], [-2.999620435470373, -4.666299570025781, -2.5178027320383265], [-1.1762160517227773, -1.4042707327700026, -4.613158686651246], [-2.067555374863428, -3.7226379503784477, -3.2092627040436525], [-2.2127120426159053, -2.233538473581601, -3.9556608181693775], [-2.4575030280158288, -5.836996165224316, -2.066948522032959], [-1.481852352148188, -0.07250914943100173, -5.299018022166564], [-1.7580363768505194, -3.6603014501501874, -2.463216545853946], [-0.5828222122506723, -4.00692794592085, -4.335390199958049], [-3.7388910307246555, -2.086920318921689, -2.783825390507779], [-2.417251708799828, -2.4134510039841195, -4.7542333367504455], [-3.3765686872853866, -7.025811770060353, -1.2548735798105923], [-2.4397389941876817, -6.114977121841196, -1.013415289123973], [-1.1470322341115562, -6.381355504429413, -3.0334832379670047], [-1.8774286235888589, -0.06680527454630547, -6.0564068197531995], [0.025724099091865682, 0.5877359083662066, -5.918471903661896], [-3.004578553272232, 0.22609329639041498, -4.260759943285649], [-1.477709975364668, 1.9897411348054526, -7.324406124231281], [-1.46336413865302, -1.7989664633730946, -9.707565340506406], [-1.797426181372085, 0.6348482256286512, -8.058897022174584], [-1.014318089750283, -0.4702212195233706, -9.033933363801694], [-2.3682192389173657, 2.5962554131344087, -6.650532506646236], [-0.5401287387343281, -2.6607328642381765, -10.614521231083685], [-3.2683473476075733, 0.1727052798053319, -7.3042204487666], [-1.6878668743765104, 0.15730837965319333, -8.981490426765141], [-1.1006640596352852, 0.15343207546184434, -8.04940448289517], [0.45592751525169106, 0.04259875713355887, -9.81235243890351], [-2.1619373088584553, 4.141208881854517, -5.758797692609178], [-3.859411554690311, 2.577228476699621, -5.569491501550019], [-2.2165839344208074, 2.442696117761212, -7.393248418515448], [1.0212347023501955, -2.1997883823831557, -11.39897850358561], [-1.514876865279195, -4.167391764477415, -11.179724405734504], [-0.48178496868532505, -2.1479809232613727, -9.78425300177377], [5.28671107442501, -5.454587224831318, -5.766210315085679], [7.0691124085362524, -6.279682887642622, -6.066255928709651], [6.229654260638579, -5.934978890987182, -5.446613288758936], [6.0113883383298425, -5.733858655977926, -6.562082241042557], [5.904086789408708, -5.797615199578466, -4.498279196303226], [6.680605162621432, -5.908993973260756, -7.287672514708348], [6.669689935224763, -7.4218388984996935, -4.457278897591813], [-7.3217434721913355, -5.354984488069682, -4.589304689864939], [4.520229202824227, -6.245675413307644, -7.401084321475437], [5.674766620705409, -4.1766772236442335, -7.447576167525742], [4.98494007297557, -5.336755906879162, -4.649069953079135], [6.060202305068211, -7.24223985966308, -3.5221246309246084], [7.178211431469265, -5.113592117361794, -3.6576164427828246], [6.1475334038084855, -4.374810915083085, 6.597369373718479], [-6.749990267166316, -6.825381383034021, -6.166687840208599], [5.179113243827639, -6.328278930650679, 6.732200977217605], [1.2659802501385293, -6.018361224083765, 7.094634512301266], [0.2446516336159128, 6.663494343182961, 3.318782005415148], [0.7300205521829497, -7.104774539753613, 5.729488597974509], [0.8387119458354135, -7.158706885068966, 4.664066947433511], [0.9963313017582465, -6.411650394247424, -7.345043092108297], [0.531483263798119, 6.801028325484113, 2.498253769716431], [-0.7495186200498374, 7.491347997853277, 5.262556972165589], [1.320297529678971, 6.496105379695201, 5.188991658035178], [0.3245715642174396, -5.695019486274559, 5.348722641264585], [2.3950329259008285, -6.838283534143148, 5.174831663605897], [1.527887630473188, -5.304739657530048, -5.796011642976046], [-0.41998113743844456, -6.475519873265583, -7.366337865266396], [1.7060315231441763, 7.40497831394485, 7.4383418157714685], [2.05395981844255, 7.283519287882171, 2.968573133514923], [-0.4207035195384776, 5.157629858002299, 0.955114801853501], [0.09984445229502015, -6.7546609877613975, 3.0751414674099617]]
#calculating parameters
dat=calcparamst(connect,atomtype)
threes=connect2three(connect,pos)
fourst=connect2fourt(connect)
foursi=connect2fouri(connect)
angleparam=calcangleparam(threes,atomtype)
nbparam=calcnb(atomtype)
convrpair=connect2rpair(connect)
radius=[]
for val in pos:
    radius.append(0.3)
plotpoints(pos,radius)

Elist=[]
n=1000
number=np.linspace(1,n,n)
for i in number:
    pos,Elistv=nextconfig(.0015,pos,dat,nbparam,convrpair,angleparam,threes,fourst,foursi)
    Elist.append(Elistv)
    #probably ultra inefficient in kcal/mol now
    #Elist.append(energy(pos,dat,nbparam,convrpair,angleparam,threes,fourst,foursi))
scale=Elist[0]
Elist2=Elist/scale
i=1
Elist3=[]
num=[]
f=0
while i < len(Elist):
    Elist3.append((Elist2[i]-Elist2[i-1]))
    num.append(i)
    i=i+1

plotpoints(pos,radius)
#print(pos)
fig, ax = plt.subplots()
ax.plot(number, Elist2)
fig, ax = plt.subplots()
ax.plot(num, Elist3)
fig.savefig("test.png")
#actually show
plt.show()

writepos(pos,atomtype)
writetofile(pos,atomtype)
print(count/n)