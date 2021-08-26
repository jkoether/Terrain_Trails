import matplotlib.pyplot as plt
# import scipy.special as spec
import triangle as tr

import numpy as np
import itertools as itertools

class loop():
    def __init__(self,pts):
        self.pts = pts
        self.edg=np.arange(pts.shape[0]-1)
        self.edg=np.transpose(np.vstack((self.edg,self.edg+1)))
        self.edg=np.vstack((self.edg,[pts.shape[0]-1,0]))
        self.inverted=np.ones(self.edg.shape[0],bool)


def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0



def RemovePoints(pts,edg,rem):
    # rem is the point numbers to remove.

    if rem.dtype=='bool':
        rem=np.flatnonzero(rem)
    edg=np.flip(np.sort(rem))
    pts=np.delete(pts,rem)
    edg=np.delete(edg,np.isin(edg[:,0],rem)|np.isin(edg[:,1],rem),axis=0)
    for r in rem:
        edg[edg>r]=edg[edg>r]-1


def split_crossings(pts,edg,inv=np.array([])):
    x_edges,x_pt=LineIntersection(pts,edg)
    #print(" *"+str(x_edges.shape[0])+" intersections")

    # EdgePlot1(pts,edg,h_edg=x_edges.flatten(),show_labels=True)
    # EdgePlot1(pts,edg,h_edg=inv,show_labels=True)
    #split edges
    if x_pt.shape[0]>0:
        if inv.shape[0]==edg.shape[0]:
            inv=np.hstack((inv,inv[x_edges[:,0]])) #inherit inverted flag.
        edg=np.vstack((edg,np.transpose(np.vstack((np.arange(x_pt.shape[0])+pts.shape[0],edg[x_edges[:,0],1])))))
        edg[x_edges[:,0],1]=np.arange(x_pt.shape[0])+pts.shape[0]

        pts=np.vstack((pts,x_pt))
    # EdgePlot1(pts,edg,h_edg=inv,show_labels=True)
    return pts,edg,inv


def ScrubNetwork(pts,edg,d_tol):
    d_tol=d_tol/10

    # Add intersections for any lines that cross
    pts,edg,_=split_crossings(pts,edg)





    #merge close points < 10% of offset.
    comb = np.vstack(list(itertools.combinations(range(pts.shape[0]), 2)))
    d=np.sqrt((pts[comb[:,0],0]-pts[comb[:,1],0])**2+(pts[comb[:,0],1]-pts[comb[:,1],1])**2)

    comb=comb[d<d_tol]
    #print(" *"+str(comb.shape[0])+" merged points")
    comb=comb[np.flip(np.argsort(comb[:,1])),:]
    for c in comb:
        pts[c[0]]=np.mean(pts[c,:],axis=0)


        edg=np.delete(edg,np.all(np.equal(edg,c),axis=1),axis=0) #remove any lines between merged points, could skip and look for self interesections later..
        edg[edg==c[1]]=c[0] #merge endpoints
        edg[edg>c[1]]=edg[edg>c[1]]-1  #renumber

    pts=np.delete(pts,comb[:,1],axis=0)


    # merge close points and lines
    d= ptDist(pts,edg)
    x_loc=np.vstack(np.where(d<d_tol))#locations where points are very close to another, merge these.
    #print(" *"+str(x_loc.shape[0])+" merged points/lines")

    for i in range(x_loc.shape[1]):
        edg=np.vstack((edg,[edg[x_loc[1,i],1],x_loc[0,i]]))
        edg[x_loc[1,i],1]=x_loc[0,i]


    # delete duplicate edges
    l1=edg.shape[0]
    edg=np.unique(np.sort(edg,axis=1), axis=0) #remove any duplicate edges.
    #print(" *"+str(l1-edg.shape[0]) + " duplicatge edges")

    # delete zero length edges
    idx=edg[:,0]==edg[:,1]
    np.delete(edg,idx,axis=0)
    #rint(" *"+str(sum(idx)) + " zero-length edges")
    # EdgePlot1(pts,edg,show_labels=True)

    #turn deadends into loops
    n=0
    for i in range(pts.shape[0]):
        if np.sum(edg==i)==1:

            n=n+1
            idx=np.argmax(np.any(edg==i,1))
            if edg[idx,1]==i:

                a=np.arctan2(np.diff(pts[edg[idx],1]),np.diff(pts[edg[idx],0])).item()
                d=np.array([np.cos(a+np.pi/2)*d_tol/2,np.sin(a+np.pi/2)*d_tol/2])
                pts=np.vstack((pts,pts[edg[idx,1]]+d))
                pts[edg[idx,1]]=pts[edg[idx,1]]-d
                edg=np.vstack((edg,[edg[idx,1].item(),pts.shape[0]-1]))
                edg=np.vstack((edg,[pts.shape[0]-1,edg[idx,0].item()]))
            else:
                a=np.arctan2(-np.diff(pts[edg[idx],1]),-np.diff(pts[edg[idx],0])).item()
                d=np.array([np.cos(a+np.pi/2)*d_tol/2,np.sin(a+np.pi/2)*d_tol/2])
                pts=np.vstack((pts,pts[edg[idx,0]]+d))
                pts[edg[idx,0]]=pts[edg[idx,0]]-d
                edg=np.vstack((edg,[edg[idx,0].item(),pts.shape[0]-1]))
                edg=np.vstack((edg,[pts.shape[0]-1,edg[idx,1].item()]))

    #print(" *"+str(n) + " deadends")

    return pts,edg

def LineIntersection(pts,edg,pts2=np.array([]),edg2=np.array([])):
    if len(edg2.shape)==1:
        edg2=np.expand_dims(edg2, 0)

    if pts2.any(): #intersections between two sets of lines, no self intersections, no existing intersections
        A=pts[edg[:,0],1]-pts[edg[:,1],1]
        B=pts[edg[:,1],0]-pts[edg[:,0],0]
        C=-A*pts[edg[:,0],0]-B*pts[edg[:,0],1]

        #if  Ax+By+C values are different signs for each point, then one splits the two points.  If same is true for transpose value then the line segments intersect.
        #this will check all combinations
        split1=np.outer(pts2[edg2[:,0],0],A)+np.outer(pts2[edg2[:,0],1],B)+C
        split2=np.outer(pts2[edg2[:,1],0],A)+np.outer(pts2[edg2[:,1],1],B)+C
        split=((split1>0) & (split2<0)) | ((split1<0) & (split2>0))

        #swap which point are used for line equations and repeat
        A2=pts2[edg2[:,0],1]-pts2[edg2[:,1],1]
        B2=pts2[edg2[:,1],0]-pts2[edg2[:,0],0]
        C2=-A2*pts2[edg2[:,0],0]-B2*pts2[edg2[:,0],1]

        #if  Ax+By+C values are different signs for each point, then one splits the two points.  If same is true for transpose value then the line segments intersect.
        #this will check all combinations
        split1=np.outer(pts[edg[:,0],0],A2)+np.outer(pts[edg[:,0],1],B2)+C2
        split2=np.outer(pts[edg[:,1],0],A2)+np.outer(pts[edg[:,1],1],B2)+C2
        split_T=((split1>0) & (split2<0)) | ((split1<0) & (split2>0))
        #check transpose values
        #lines intersect when one sets of points are  split by the ocrresponding line segment and vise versa.
        split=split&np.transpose(split_T) #If transpose is not true then only one set of points is split by other line, lines do not cross.

    	#rows and columns of true elements, correspond to intersecting edge numbers.  will include a duplicate for each.
    	#may also include some already intersecting lines as well as self intersections.


        if np.any(split):

            #remove duplicates
            x_edges=np.unique(np.sort(np.vstack(np.nonzero(split)),axis=0),axis=1)


            x_edges=np.transpose(x_edges)


            x_pt=np.zeros(x_edges.shape)
            for i in range(x_edges.shape[0]):
                x_pt[i,:]=[(B2[x_edges[i,0]]*C[x_edges[i,1]]-B[x_edges[i,1]]*C2[x_edges[i,0]])/((A2[x_edges[i,0]]*B[x_edges[i,1]]-A[x_edges[i,1]]*B2[x_edges[i,0]])),
                           (A[x_edges[i,1]]*C2[x_edges[i,0]]-A2[x_edges[i,0]]*C[x_edges[i,1]])/((A2[x_edges[i,0]]*B[x_edges[i,1]]-A[x_edges[i,1]]*B2[x_edges[i,0]]))]
            return x_edges,x_pt
        else:
            return np.array([]),np.array([])


    else: #self intersections
        A=pts[edg[:,0],1]-pts[edg[:,1],1]
        B=pts[edg[:,1],0]-pts[edg[:,0],0]
        C=-A*pts[edg[:,0],0]-B*pts[edg[:,0],1]

        # split=(np.outer(pts[edg[:,0],0],A)+np.outer(pts[edg[:,0],1],B)+C<0)!=(np.outer(pts[edg[:,1],0],A)+np.outer(pts[edg[:,1],1],B)+C<0)
        #if  Ax+By+C values are different signs for each point, then one splits the two points.  If same is true for transpose value then the line segments intersect.
        #this will check all combinations
        split1=np.outer(pts[edg[:,0],0],A)+np.outer(pts[edg[:,0],1],B)+C
        split2=np.outer(pts[edg[:,1],0],A)+np.outer(pts[edg[:,1],1],B)+C
        split=((split1>0) & (split2<0)) | ((split1<0) & (split2>0))
        #check transpose values
        #lines intersect when one sets of points are  split by the ocrresponding line segment and vise versa.
        split=split&np.transpose(split) #If transpose is not true then only one set of points is split by other line, lines do not cross.

    	#rows and columns of true elements, correspond to intersecting edge numbers.  will include a duplicate for each.
    	#may also include some already intersecting lines as well as self intersections.



        if np.any(split):

            #remove duplicates
            x_edges=np.unique(np.sort(np.vstack(np.nonzero(split)),axis=0),axis=1)

        	#remove predefined intersections, look for any duplicates between intersecting edges (may not be needed)
            #this is probably not needed with the updated logic.
            x_edges=x_edges[:,np.all(np.diff(np.sort(np.concatenate((edg[x_edges]),axis=1), axis=1), axis=1) !=0,axis=1)]

            x_edges=np.transpose(x_edges)
            x_pt=np.zeros(x_edges.shape)
            for i in range(x_edges.shape[0]):
                x_pt[i,:]=[(B[x_edges[i,0]]*C[x_edges[i,1]]-B[x_edges[i,1]]*C[x_edges[i,0]])/((A[x_edges[i,0]]*B[x_edges[i,1]]-A[x_edges[i,1]]*B[x_edges[i,0]])),
                   (A[x_edges[i,1]]*C[x_edges[i,0]]-A[x_edges[i,0]]*C[x_edges[i,1]])/((A[x_edges[i,0]]*B[x_edges[i,1]]-A[x_edges[i,1]]*B[x_edges[i,0]]))]
            return x_edges,x_pt
        else:
            return np.array([]),np.array([])

def aCorr(a):
    #correction needed for inputs for arccos, numerical errors for vert lines
    tol=10**-6
    a[(a>1) & (a<1+tol)]=1.0
    a[(a<-1) & (a>(-1-tol))]=-1.0
    return a


def ptDist(pts,edg,p2=np.array([]),infLength=False):
    #get distance between a set of points, and a set edges/points
    #if no second set of points input, compare lines to points and return nan when point is referenced by edge.

    if edg.ndim==1:
        edg=np.array([edg])
    if p2.size == 0: #intersections between two
        p2=pts
        p_num=np.transpose(np.tile(np.arange(p2.shape[0]),(edg.shape[0],1)))
        on_line=(np.tile(edg[:,0],(p2.shape[0],1))==p_num)|(np.tile(edg[:,1],(p2.shape[0],1))==p_num)
    else:
        on_line=np.zeros((p2.shape[0], edg.shape[0]),dtype=bool)

    #Law of cosines for triangle with 2 points of original line (A,B) and offset point (C). (line a is opposite angle A)
    c=np.tile(np.sqrt(np.sum((pts[edg[:,0]]-pts[edg[:,1]])**2,axis=1)),(p2.shape[0],1)) #between two line segment points
    a=np.sqrt(np.sum((np.tile(pts[edg[:,1]],(p2.shape[0],1,1))-np.swapaxes(np.tile(p2,(edg.shape[0],1,1)),0,1))**2,axis=2))
    b=np.sqrt(np.sum((np.tile(pts[edg[:,0]],(p2.shape[0],1,1))-np.swapaxes(np.tile(p2,(edg.shape[0],1,1)),0,1))**2,axis=2))


    A=np.arccos(aCorr((b**2+c**2-a**2)/(2*b*c)))
    B=np.arccos(aCorr((a**2+c**2-b**2)/(2*a*c)))

    d=np.full(A.shape, -1,float)

    if not infLength:
        #if either angle to offset point from original line is > 90 deg, the point
        #is past one end of the line, in that case just use the min distance to one of the end points.
        d[A>np.pi/2]=b[A>np.pi/2]
        d[B>np.pi/2]=a[B>np.pi/2]

    #point distance to line segment.
    d[d<0]=a[d<0]*np.sin(B[d<0])

    d[on_line]=np.inf #!!!should this be zero?

    return d

def MaxWidth(pts,offset):
    edg=seqEdges(pts.shape[0])

    comb = np.vstack(list(itertools.combinations(range(pts.shape[0]), 2)))
    d=np.sqrt((pts[comb[:,0],0]-pts[comb[:,1],0])**2+(pts[comb[:,0],1]-pts[comb[:,1],1])**2)

    d2=ptDist(pts,edg)
    w=np.where(d2<offset)
    #x_edg,x_pt=LineIntersection(pts,edg,pts,comb)
    for i in range(d2.shape[0]):
        idx=np.flatnonzero(d2[i,:]<offset)
        if np.any(idx):
            x=1

def PolygonArea(pts):
    n = len(pts) # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += pts[i][0] * pts[j][1]
        area -= pts[j][0] * pts[i][1]

    return area / 2.0


def EdgePlot1(pts,edg=np.array([],int),h_edg=np.array([],int),h_pts=np.array([],int),title="",show_labels=False,label_nums=[]):
    if h_edg.dtype==bool:
        h_edg=np.flatnonzero(h_edg)

    if h_pts.dtype==bool:
        h_pts=np.flatnonzero(h_pts)

    #highlight_edges=np.flatnonzero(inverted)
    fig4, ax = plt.subplots()

    if not edg.any():
        edg=seqEdges(pts.shape[0])
    ax.plot(pts[:,0],pts[:,1], 'o', color='black',markersize=3);
    ax.plot(pts[h_pts,0],pts[h_pts,1], 'o',markersize=8,color='red',alpha=0.4);
    if show_labels:
        for i in range(pts.shape[0]):
            plt.annotate(str(i), # this is the text
                     (pts[i,0],pts[i,1]))
        for i in range(edg.shape[0]):
            plt.annotate(str(i),
                     (np.mean(pts[edg[i],0]),np.mean(pts[edg[i],1])),
                     color='blue')
    if np.any(label_nums):
        label_nums=(label_nums*180/np.pi).astype(int)
        for i in range(pts.shape[0]):
            plt.annotate(str(label_nums[i]), # this is the text
                     (pts[i,0],pts[i,1]))
    for i in range(edg.shape[0]):
        #if np.sqrt(np.sum((points3[e[0],:]-points3[e[1],:])**2))>=offset*2:
        ax.plot([pts[edg[i,0],0],pts[edg[i,1],0]],[pts[edg[i,0],1],pts[edg[i,1],1]],'g-')
        if np.isin(i,h_edg):
            ax.plot([pts[edg[i,0],0],pts[edg[i,1],0]],[pts[edg[i,0],1],pts[edg[i,1],1]],'b-',linewidth=4,alpha=0.6)

    plt.title(title)
    plt.axis('equal')
    plt.show()
    plt.pause(0.1)


def EdgePlot2(pts,edg,pts2,edg2,h_edg=np.array([],int),h_pts=np.array([],int),title="",show_labels=False,holes=np.array([])):

    if h_edg.dtype==bool:
        h_edg=np.flatnonzero(h_edg)

    if h_pts.dtype==bool:
        h_pts=np.flatnonzero(h_pts)

    fig4, ax = plt.subplots()

    if np.any(holes):
        ax.plot(holes[:,0],holes[:,1], '.', color='red',markersize=1);
    ax.plot(pts[:,0],pts[:,1], 'o', color='blue',markersize=3);
    ax.plot(pts2[:,0],pts2[:,1], 'o', color='black',markersize=3);

    ax.plot(pts2[h_pts,0],pts2[h_pts,1], 'o',markersize=8,color='blue',alpha=0.4);

    if show_labels:
        for i in range(pts.shape[0]):
            plt.annotate(str(i),
                     (pts[i,0],pts[i,1]))

        for i in range(pts2.shape[0]):
            plt.annotate(str(i),
                     (pts2[i,0],pts2[i,1]))

        # for i in range(edg.shape[0]):
        #     plt.annotate(str(i),
        #              (np.mean(pts[edg[i],0]),np.mean(pts[edg[i],1])),
        #              color='blue')

        for i in range(edg2.shape[0]):
            plt.annotate(str(i),
                     (np.mean(pts2[edg2[i],0]),np.mean(pts2[edg2[i],1])),
                     color='blue')


    for i in range(edg.shape[0]):
        #if np.sqrt(np.sum((points3[e[0],:]-points3[e[1],:])**2))>=offset*2:
        ax.plot([pts[edg[i,0],0],pts[edg[i,1],0]],[pts[edg[i,0],1],pts[edg[i,1],1]],'r-')


    for i in range(edg2.shape[0]):
        #if np.sqrt(np.sum((points3[e[0],:]-points3[e[1],:])**2))>=offset*2:
        # ax.plot([pts2[edg2[i,0],0],pts2[edg2[i,1],0]],[pts2[edg2[i,0],1],pts2[edg2[i,1],1]],'k-')
        ax.plot([pts2[edg2[i,0],0],pts2[edg2[i,1],0]],[pts2[edg2[i,0],1],pts2[edg2[i,1],1]],'k-')
        if np.isin(i,h_edg):
            #ax.plot([pts2[edg2[i,0],0],pts2[edg2[i,1],0]],[pts2[edg2[i,0],1],pts2[edg2[i,1],1]],'g-')
            ax.plot([pts2[edg2[i,0],0],pts2[edg2[i,1],0]],[pts2[edg2[i,0],1],pts2[edg2[i,1],1]],'b-',linewidth=4,alpha=0.6)

    plt.title(title)
    plt.axis('equal')
    plt.show()
    plt.pause(0.1)

def splitLoops(pts,inverted):
    loops=[pts]
    inverted=[inverted]
    split=True
    # n=0#debugging
    while split:
        # if n==10:
        #     x=1
        split=False
        for i in range(len(loops)):
            # print(i)
            edg=seqEdges(loops[i].shape[0])
            # EdgePlot1(loops[i],edg,h_pts=inverted[i],show_labels=True)
            for j in range(loops[i].shape[0]-1):
                # print("-"+str(j))
                x_edg,x_pt=LineIntersection(loops[i],edg,loops[i][j:j+2,:],np.array([0,1]))
                if x_edg.any():
                   split=True
                   loops.append(np.vstack((loops[i][0:j+1,:],x_pt[0,:],loops[i][x_edg[0,1]+1:,:]))) #pinch off the closed loop
                   loops[i]=np.vstack((x_pt[0,:],loops[i][j+1:x_edg[0,1]+1,:]))
                   inverted.append(np.hstack((inverted[i][0:j+1],[False],inverted[i][x_edg[0,1]+1:]))) #pinch off the closed loop
                   inverted[i]=np.hstack(([False],inverted[i][j+1:x_edg[0,1]+1]))
                   break
            if split:
                # n=n+1
                break
        # for i in range(len(loops)):
        #     print(str(loops[i].shape[0]) + " - " + str(inverted[i].shape[0]))

    return loops,inverted

def seqEdges(n):
    edg=np.arange(n-1)
    edg=np.transpose(np.vstack((edg,edg+1)))
    edg=np.vstack((edg,[n-1,0]))
    return edg







def getOffsetPoints(loop,offset):
    pts=np.empty([loop.shape[0],2])
    for i in range(loop.shape[0]):

        d=np.diff(loop[np.mod([i-1,i,i+1],loop.shape[0]),:],axis=0)
        d[0,:]=-d[0,:]
        a=np.arctan2(d[:,1],d[:,0])
        if a[0]>a[1]: #a[1] must be greater than a[0]
            a[1]=a[1]+2*np.pi
        a2=np.mean(a)
        if np.sin(a2-a[0])==0:
            x=1
        pts[i,:]=loop[i,:]+[np.cos(a2),np.sin(a2)]*np.asarray(offset/np.sin(a2-a[0]))
    return pts

def FourPtInt(pts):
#check if first and last (third) lines made from 4 pts are crossed.

    A1=pts[0,1]-pts[1,1]
    B1=pts[1,0]-pts[0,0]
    C1=-A1*pts[0,0]-B1*pts[0,1]

    A2=pts[2,1]-pts[3,1]
    B2=pts[3,0]-pts[2,0]
    C2=-A2*pts[2,0]-B2*pts[2,1]
    if A1*B2-A2*B1 == 0:
        x=1
    return np.asarray([(B1*C2-B2*C1)/((A1*B2-A2*B1)),(A2*C1-A1*C2)/((A1*B2-A2*B1))])


def RefPtDist(pts,refPt):
    return np.sqrt(np.sum((pts-refPt)**2,axis=1))

def CrossCheck(pts,extend1=False):
#check if first and last (third) lines made from 4 pts are crossed.

#Extend1 checks if either one the infinite length lines interesects the other finite length lines.
    A1=pts[0,1]-pts[1,1]
    B1=pts[1,0]-pts[0,0]
    C1=-A1*pts[0,0]-B1*pts[0,1]

    A2=pts[2,1]-pts[3,1]
    B2=pts[3,0]-pts[2,0]
    C2=-A2*pts[2,0]-B2*pts[2,1]

    split1=A1*pts[2:,0]+B1*pts[2:,1]+C1
    split2=A2*pts[0:2,0]+B2*pts[0:2,1]+C2

    split1=((split1[0]>0)&(split1[1]<0))|((split1[0]<0)&(split1[1]>0))
    split2=((split2[0]>0)&(split2[1]<0))|((split2[0]<0)&(split2[1]>0))

    if extend1:
        if split1 | split2:
            return np.asarray([(B1*C2-B2*C1)/((A1*B2-A2*B1)),(A2*C1-A1*C2)/((A1*B2-A2*B1))])
        else:
            return np.asarray([])
    else:
        if split1 & split2:
            return np.asarray([(B1*C2-B2*C1)/((A1*B2-A2*B1)),(A2*C1-A1*C2)/((A1*B2-A2*B1))])
        else:
            return np.asarray([])

def InvertCheck(org_pts,loop,pts):
    inverted=np.zeros(pts.shape[0],bool)
    loop=np.squeeze(org_pts[loop,:])

    #look for offset edges that are inverted

    for i in range(loop.shape[0]):
        idx= np.squeeze([i-1,i])%loop.shape[0]
        #get two angles to from one base point to two offset points these
        #should always be in CCW order unless they are inverted
        # inverted[i]=np.diff(np.arctan2(org_pts[l_edg[i],1]-pts[idx,1],org_pts[l_edg[i],0]-pts[idx,0]))>0
        p=np.vstack((loop[idx],pts[idx,:]))

        inverted[i]=CrossCheck(p[[0,2,1,3],:]).shape[0]>0
    return inverted

def isCCW(pts):
    d=np.diff(np.vstack((pts,pts[0,:])),axis=0)
    return np.sum(d[:,0]*d[:,1])>0

def getAngles(pts):
    d=np.vstack((pts[-1,:],pts))-np.vstack((pts,pts[0,:]))
    a=np.arctan2(d[:,1],d[:,0])
    a=(np.diff(a)+np.pi)%(2*np.pi)
    a[a<0]=a[a<0]+2*np.pi
    return a

def absAngles(pts):
    #tuple input
    pts=np.vstack(pts)
    d=pts[0:-1,:]-pts[1:,:]
    d=d[[0,-1]]
    d[-1]=-d[-1]
    a=np.arctan2(d[:,1],d[:,0])
    # a=(np.diff(a)+np.pi)%(2*np.pi)
    a=a[[0,-1]]
    if a[0]>a[1]:
        a[1]=a[1]+2*np.pi
    return a

def getNewPts(pt1,pt2,offset,lim_angle):

    # if lim_angle > pi, we are truncating exterior corners.
    intPt=FourPtInt(np.vstack((pt1,pt2)))
    d1=RefPtDist(pt1,intPt)
    d2=RefPtDist(pt2,intPt)

    x_pt=CrossCheck(np.vstack((pt1,pt2)),extend1=True)


    #if either projected lines intersects the other edge,
    #or both projected lines do and intersection is closer to inner points
    if np.any(x_pt) or (np.diff(d1)<0 and np.diff(d2)>0):
        a=absAngles((pt1,pt2))
        if (lim_angle<np.pi and (np.diff(a) > lim_angle)) or (lim_angle>np.pi and (np.diff(a) < lim_angle)):
            return intPt
        else:
            l=(offset/2)/np.sin(np.diff(a)/2)
            return intPt+l*(np.vstack([np.cos(a),np.sin(a)])).T
    else:
        newPt=np.array([])
        if np.diff(d1)<0:
            newPt=np.hstack((newPt,pt1[1,:]))
        if np.diff(d2)>0:
            newPt=np.hstack((newPt,pt2[0,:]))
        if not np.any(newPt):
            x=1
        return newPt

def getTrucPts(pts,idx,offset, lim_angle):
    #if angle > pi (180), then it checks external corners.
    c_pt=np.array([idx-1,idx+1])%pts.shape[0]
    corner=pts[idx,:]
    if c_pt[0]>600 and c_pt[0]<605:
        x=1
    intPt=corner


    a=np.diff(absAngles((pts[c_pt[0]],intPt,pts[c_pt[1]]))).item()

    d1=np.sqrt(np.sum((pts[c_pt[0],:]-intPt)**2))
    d2=np.sqrt(np.sum((pts[c_pt[1],:]-intPt)**2))

    min_dist=2*np.sin(a/2)*min(d1,d2) #max truncated distance using the closest endpoint

    while min_dist<offset and ((lim_angle<np.pi and a<lim_angle) or (lim_angle>np.pi and a>lim_angle)):
        if c_pt[0]==604:
            x=1
        # EdgePlot1(pts,h_pts=c_pt)
        # find shortest edge
        if c_pt[0]==(c_pt[1]+1)%pts.shape[0]:
            return -1,[]
        if d1<d2:
            c_pt[0]=(c_pt[0]-1)%pts.shape[0]
            #distance from two new edge points to the other edge line
            # if the first point is closer, than the new segment getting closer and we need to skip ahead to avoid errors.
            idx1=np.array([c_pt[1]-1,c_pt[1]])%pts.shape[0]
            idx2=np.array([c_pt[0]+1,c_pt[0]])%pts.shape[0]
            if np.diff(ptDist(pts[idx1,:],np.array([0,1]),p2=pts[idx2,:],infLength=True),axis=0)<0:
                intPt=corner
            else:
                intPt=FourPtInt(pts[np.hstack((idx1,idx2)),:])

        else:
            c_pt[1]=(c_pt[1]+1)%pts.shape[0]
            idx1=np.array([c_pt[0]+1,c_pt[0]])%pts.shape[0]
            idx2=np.array([c_pt[1]-1,c_pt[1]])%pts.shape[0]
            if np.diff(ptDist(pts[idx1,:],np.array([0,1]),p2=pts[idx2,:],infLength=True),axis=0)<0:
                intPt=corner
            else:
                intPt=FourPtInt(pts[np.hstack((idx1,idx2)),:])

        if c_pt[0]==602:
            x=1
        if c_pt[0]==(c_pt[1]+1)%pts.shape[0]:
            return -1,[]
        a=np.diff(absAngles((pts[c_pt[0]],intPt,pts[c_pt[1]]))).item()

        d1=np.sqrt(np.sum((pts[c_pt[0],:]-intPt)**2))
        d2=np.sqrt(np.sum((pts[c_pt[1],:]-intPt)**2))
        min_dist=2*np.sin(a/2)*min(d1,d2) #max truncated distance using the closest endpoint

    idx1=np.array([c_pt[0],c_pt[0]+1])%pts.shape[0]
    idx2=np.array([c_pt[1]-1,c_pt[1]])%pts.shape[0]
    return c_pt,getNewPts(pts[idx1,:],pts[idx2,:],offset,lim_angle)





def CleanLoop(org_pts,org_edg,loops,offset,inverted,lim_angle,debug=False):
    tol=.001*offset
    del_Loop=False
    #delete inverted points and eliminate 2 point loops
    for i in range(len(loops)-1,-1,-1):
          loops[i]=np.delete(loops[i],inverted[i],axis=0)
          if loops[i].shape[0]<3:
              loops.pop(i)

    # truncate sharp internal corners
    for i in range(len(loops)-1,-1,-1):
        # print(' -'+str(i))
        # EdgePlot1(loops[i])
        #get all angles
        a=getAngles(loops[i])
        while any(a<lim_angle) or  any(a>(2*np.pi-lim_angle)):
            # EdgePlot1(loops[i],label_nums=a)
            if any(a<lim_angle): #check for internal angle first (arbitrary)
                idx=np.argmax(a<lim_angle)
                cutPts,newPts=getTrucPts(loops[i],idx,offset,lim_angle)
            else:
                idx=np.argmax(a>(2*np.pi-lim_angle))
                cutPts,newPts=getTrucPts(loops[i],idx,offset,2*np.pi-lim_angle)
                # if np.all(cutPts!=-1):
                #     if np.diff(cutPts)<0:
                #         idx=np.arange(cutPts[0]-1,cutPts[1]+loops[i].shape[0]+2)%loops[i].shape[0]
                #     else:
                #         idx=np.arange(cutPts[0]-1,cutPts[1]+2)%loops[i].shape[0]
                #     fig4, ax = plt.subplots()

                #     ax.plot(loops[i][idx,0],loops[i][idx,1],'-k')
                #     if np.any(newPts):
                #         p=np.vstack((loops[i][cutPts[0],:],newPts,loops[i][cutPts[1],:]))
                #     else:
                #         p=np.vstack((loops[i][cutPts[0],:],loops[i][cutPts[1],:]))
                #     ax.plot(p[:,0],p[:,1], '-r')
                #     # ax.plot(loops[i][h_pts,0],pts[h_pts,1], 'o',markersize=8,color='red',alpha=0.4);
                #     plt.axis('equal')
                #     plt.title(str(cutPts[0])+" / "+str(cutPts[1]))

                #     plt.show()
                #     plt.pause(0.1)
                #     x=1
            if np.any(cutPts==-1): #loop too small
                del_Loop=True
                loops.pop(i)
                a=0 #to skip next while loop
                break
            if cutPts[0]>cutPts[1]: #wrapped around
                if np.any(newPts):
                    loops[i]=np.vstack((loops[i][cutPts[1]:cutPts[0]+1,:],newPts))
                else:
                    loops[i]=np.vstack((loops[i][cutPts[1]:cutPts[0]+1,:]))
            else:
                if np.any(newPts):
                    loops[i]=np.vstack((loops[i][0:cutPts[0]+1,:],newPts,loops[i][cutPts[1]:,:]))
                else:
                    loops[i]=np.vstack((loops[i][0:cutPts[0]+1,:],loops[i][cutPts[1]:,:]))
            a=getAngles(loops[i])




    for i in range(len(loops)):
        if loops[i].shape[0]>2:

            edg2=seqEdges(loops[i].shape[0])
            d=np.min(ptDist(org_pts,org_edg,loops[i]),axis=1)
            # d2=np.min(ptDist(loops[i],edg2,org_pts),axis=0)
            x_edg,x_pt=LineIntersection(org_pts,org_edg,loops[i],edg2)
            # d3=ptDist(loops[i],edg2)

            hp=(d<offset*.999)

            # if x_edg.shape[0]>0:
            #     xb=np.zeros(edg2.shape[0],bool)
            #     xb[x_edg[:,0]]=True
            #     he=xb | (d2<offset*.999)
            # else:
            #     he=(d2<offset*.999)
            # EdgePlot2(org_pts,org_edg,loops[i],edg2,h_edg=he,h_pts=hp)
            loops[i]=np.delete(loops[i],hp,axis=0)



    for i in range(len(loops)-1,-1,-1):
        if loops[i].shape[0]<3:
            loops.pop(i)
        elif loops[i].shape[0]==3:
            d=ptDist(loops[i],seqEdges(loops[i].shape[0]))
            d[d==np.inf]=0
            if np.min(d) < offset*3:
                loops.pop(i)
        elif loops[i].shape[0]<7:
            d=ptDist(loops[i],seqEdges(loops[i].shape[0]))
            d[d==np.inf]=0
            if np.mean(d) < offset*3:
                loops.pop(i)
        # else:
        #     EdgePlot2(org_pts,org_edg,loops[i],seqEdges(loops[i].shape[0]))
    holes=[]


    for i in range(len(loops)):
        #maybe this doesn't matter?  holes should always be to the "right"
         direction=PolygonArea(loops[i])>0 #CCW = true
         #find longest edge
         idx=np.argmax(np.sqrt(np.sum(np.diff(np.vstack((loops[i],loops[i][0,:])),axis=0)**2,axis=1)))
         idx=np.arange(idx,idx+2)%loops[i].shape[0]
         a=np.arctan2(np.diff(loops[i][idx,1]),np.diff(loops[i][idx,0])).item()
         #get mid point
         #project point out perpindicular from midpoint (.1*offset)
         holes.append(np.mean(loops[i][idx,:],axis=0)+0.02*offset*np.asarray([np.cos(a-np.pi/2),np.sin(a-np.pi/2)]))
         #add to hole list.





    return loops,holes

def getLoops(pts,edg,offset):
    completed=np.zeros(edg.shape,dtype=bool)
    line_len=np.sqrt((pts[edg[:,1],0]-pts[edg[:,0],0])**2+(pts[edg[:,1],1]-pts[edg[:,0],1])**2)
    n=0
    loops_all=[]
    holes_all=[]
    while not np.all(completed):
        closed=False
        # print("n="+str(n))

        #find the longest uncompleted edge
        idx = np.where(completed==False)
        # idx2=np.argmax(line_len[idx[0]])
        # base_edge=idx[0][idx2]
        base_edge=idx[0][0]
        loop_edges=[]
        loop=[]

        if idx[1][0]==0:
            base_pts=edg[base_edge,:]
            rev=False
        else:
            rev=True
            base_pts=np.flip(edg[base_edge,:])
        while not closed:


            #print(str(base_edge)+" (" +str(base_pts[0])+", " +str(base_pts[1])+")")
            conn_edges=np.flatnonzero(np.any(np.isin(edg,base_pts[1]),1))
            if rev:
                loop_edges.append(np.flip(edg[base_edge,:]))
                loop.append(edg[base_edge,0])
            else:
                loop_edges.append(edg[base_edge,:])
                loop.append(edg[base_edge,1])

            if conn_edges.shape[0]>2: #junction point, need to find next connection CCW from base edge
                order=np.hstack((np.argmax(conn_edges==base_edge),np.flatnonzero(conn_edges!=base_edge)))
                conn_edges=conn_edges[order]
                conn_points=edg[conn_edges]
                conn_points=conn_points[conn_points!=base_pts[1]]
                conn_angles=np.arctan2(pts[conn_points,1]-pts[base_pts[1],1],pts[conn_points,0]-pts[base_pts[1],0])
                conn_angles[conn_angles<conn_angles[0]]=conn_angles[conn_angles<conn_angles[0]]+2*np.pi
                idx=np.argmin(conn_angles[1:])+1
                conn_edges=conn_edges[idx]
                conn_points=conn_points[idx]

            else:
                conn_edges=conn_edges[conn_edges!=base_edge].item()
                conn_points=edg[conn_edges]
                conn_points=conn_points[conn_points!=base_pts[1]]
            base_pts=np.hstack((base_pts,conn_points)) #includes base edge plus next point.

            if rev:
                completed[base_edge,1]=True
            else:
                completed[base_edge,0]=True

            #get next edge ready
            base_edge=conn_edges
            if base_pts[1]==edg[conn_edges,0]: #next points are in correct order
                base_pts=edg[base_edge,:]
                rev=False

            else: ##next edge is in reverse direction
                base_pts=np.flip(edg[base_edge,:])
                rev=True

            if (rev and completed[base_edge,1]) or (not rev and completed[base_edge,0]):
                   closed=True
        loop_edges=np.vstack(loop_edges)
        loop=np.vstack(loop)
        if n==5:
            x=1
        if np.squeeze(pts[loop,:]).shape[0]>=3:
            pts2=getOffsetPoints(np.squeeze(pts[loop,:]),offset)
            inverted=InvertCheck(pts,loop,pts2)
            loops,inverted=splitLoops(pts2,inverted)
            # if n==5:
            #     loops,holes=CleanLoop(pts,edg,loops,offset,inverted,debug=True)
            # else:
            loops,holes=CleanLoop(pts,edg,loops,offset,inverted,np.pi/3)
            loops_all.extend(loops)
            holes_all.extend(holes)
        n=n+1



    n=0
    edges2=[]
    for i in range(len(loops_all)):
        edges2.append(n+seqEdges(loops_all[i].shape[0]))
        n=n+loops_all[i].shape[0]
    edges2=np.vstack(edges2)
    points2=np.vstack(loops_all)
    holes_all=np.vstack(holes_all)

    #EdgePlot2(pts,edg,points2,edges2,holes=holes_all)

    return points2,edges2,holes_all

def InflateNetwork(points,edges,offset):

    #print("Cleaning up input network")

    points,edges=ScrubNetwork(points,edges,offset)

    # EdgePlot1(points,edges,show_labels=True)



    pts2,edg2,holes=getLoops(points,edges,offset)
# x_edges,x_pt=LineIntersection(pts2,edg2)
# EdgePlot1(pts2,edg2,h_edg=x_edges)
    face = {
      "vertices": pts2,
      "segments": edg2,
      "holes": holes}
    t = tr.triangulate(face,'p')

    #tr.compare(plt, face, t)


    return t




# verts1=np.asarray([[0,0],[1,5],[2,0],[4,3],[5,8],[4.5,13],[5.5,13],[9,0],[11,3],[1.5,9],
#                     [11.5,3.5],[11.5,4],[11,4.2],[10.5,4.5],[.9,9],[-3.5,-0.5],[-3.3,-0.8],[6,3],[6,3.04],[7,2.5],
#                     [7.7,1],[7.8,0],[-2,10],[-2,2.8],[-3,1.5],[9.5,4],[9.3,3.4],[12,1]]) #[6.5,11],[7.5,9]
# lines1=np.asarray([[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,4],[4,7],[7,8],[1,9],[8,10],[10,11],[11,12],[12,13],[9,14],
#                     [14,15],[15,16],[3,17],[17,18],[18,19],[19,20],[20,21],[0,16],[14,22],[22,23],
#                     [23,24],[13,25],[25,26],[26,27]]) #[6,20],[6,21]
# msh=InflateNetwork(verts1,lines1,.5)

# # verts1=np.asarray([[0,0],[1,5],[2,0],[0.5,4],[2,3],[0.5,9],[4,2],[2.5,8],[2.5,2],[1.5,6],[4.5,6],[0,6],[1,6]])
# # lines1=np.asarray([[0,1],[1,2],[3,4],[5,6],[7,8],[9,10],[11,12]])




# hole=[[1.8,4.7]]

# # face = {
# #   "vertices": verts1,
# #   "segments": lines1}
# # tr.compare(plt, face, face)


# face = {
#   "vertices": verts1,
#   "segments": lines1}
# t = tr.triangulate(face,'p')


# # face = {
# #   "vertices": [[0,1],[0,11],[12,9],[10,-1],[2,2],[5,7],[7,1]],
# #   "segments": [[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,4]],
# #   "holes": [[5,5]]}
# # t = tr.triangulate(face,'p')


# tr.compare(plt, face, t)
# plt.show()




# msh.export('path.stl')



# verts2=np.asarray([[0,0],[1,5],[2,0],[0.5,4],[2,3],[0.5,9],[4,2],[2.5,8],[2.5,2],[1.5,6],[4.5,6],[0,6],[1,6]])
# lines2=np.asarray([[0,1],[1,2],[3,4],[5,6],[7,8],[9,10],[11,12]])
# SplitLines(verts2,lines2)

