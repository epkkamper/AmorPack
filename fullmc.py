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
    rmin=L
    remx=xtemp
    remy=ytemp
    remz=ztemp
    for pair in rcheck:
        rtest=dist(pos1,pair)
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
                     r23=dist(pos[ans[z+1]-1],pos[ans[z+2]-1])
                     if r23 > rmax:
                         rmax=r23
                         rmaxpos=[ans[z],ans[z+1],ans[z+2]]
                     z=z+3
                 z=0
                 rmindin=-10
                 while z<len(ans):
                     rag=dist(pos[ans[z+1]-1],pos[ans[z+2]-1])
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
                 r23=dist(pos[f[z+1]-1],pos[f[z+2]-1])
                 if r23 > rmax:
                     rmax=r23
                     rmaxpos=[ans[z],ans[z+1],ans[z+2]]
                 z=z+3
             while z<len(ans):
                 rag=dist(pos[ans[z+1]-1],pos[ans[z+2]-1])
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
    thetad=math.acos((vecn1[0]*vecn2[0]+vecn1[1]*vecn2[1]+vecn1[2]*vecn2[2])/((vecn1[0]**2+vecn1[1]**2+vecn1[2]**2)**(1/2)*(vecn2[0]**2+vecn2[1]**2+vecn2[2]**2)**(1/2)))
    return thetad
#2-d Distance calculation
def dist(pos1,pos2):
    x1=pos1[0]
    x2=pos2[0]
    y1=pos1[1]
    y2=pos2[1]
    z1=pos1[2]
    z2=pos2[2]
    dist=((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)**(1/2)
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
              minr.append(dist(atom,newval))
           j=j+1
       i=i+1
    return minr

#calculate inversion angle
def inverang(pos1,pos2,pos3,pos4):
    vec1=[pos1[0]-pos2[0],pos1[1]-pos2[1],pos1[2]-pos2[2]]
    vec2=[pos2[0]-pos3[0],pos2[1]-pos3[1],pos2[2]-pos3[2]]
    vec3=[pos2[0]-pos4[0],pos2[1]-pos4[1],pos2[2]-pos4[2]]
    vecn=[(vec2[1]*vec3[2]-vec2[2]*vec3[1]),(vec2[0]*vec3[2]-vec2[2]*vec3[0]),(vec2[0]*vec3[1]-vec2[1]*vec3[0])]
    
    thetar=math.asin((vec1[0]*vecn[0]+vec1[1]*vecn[1]+vec1[2]*vecn[2])/((vec1[0]**2+vec1[1]**2+vec1[2]**2)**(1/2)*(vecn[0]**2+vecn[1]**2+vecn[2]**2)**(1/2)))
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
L=15
count=0
pos=[[8.9348,5.8488,4.9916],[5.4583,4.8609,4.9918],[7.8384,4.9466,4.9995],[6.5548,5.7630,4.9995],[10.1756,5.1579,5.0112],[4.2179,5.5521,5.014],[7.87,4.2894,4.1078],[7.8761,4.3002,5.8986],[6.5233,6.4205,4.1079],[6.5173,6.4098,5.8984],[10.9739,5.9112,4.9984],[10.2815,4.4955,4.1293],[10.2733,4.5338,5.9214],[4.1199,6.1762,5.9216],[3.4193,3.7987,4.9984],[4.1118,6.2147,4.1297]]
#which atom is bonded with which
#specify in order (no duplicate bond description)
#also have to specify the bond order as well
connect=[[3,1,5,1],[4,1,6,1],[4,1,7,1,8,1],[9,1,10,1],[11,1,12,1,13,1],[14,1,15,1,16,1],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0]]
#order=[[1,1],[0],[0],[0],[0],[1],[0],[0],[0]]
#need atom type for atom parameter
#assuming all tetrahedral carbon atoms with index of 1
atomtype=[3,3,1,1,1,1,2,2,2,2,2,2,2,2,2,2]
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
n=15000
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
print(pos)
writetofile(pos,atomtype)
print(count/n)