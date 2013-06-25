import sys
import matplotlib.pyplot as plt
from numpy import *

def main():
    
    if len(sys.argv)<1:
        print "error! nodes location  data required!"
    else:
        fn =sys.argv[1]

    data=[]

    try:
        f = open( fn, 'r' )
    except IOError:
        print "Can not open file %s for reading"%fn
        raise
    except:
        print "Unexpected exception", sys.exc_info()[0]
        raise

    for line in f.readlines():
        arr=[]
        sub=line.split()
        arr.append(float(sub[0]))
        arr.append(float(sub[1]))
        data.append(arr)
    f.close()
    plt.figure()    
    locs=array(data)

    labels=['{0}'.format(int(i)) for i in range(len(locs))]
    #labels=['A','B','C']
    x=locs[:,0]
    y=locs[:,1]
    plt.scatter(x,y)
    plt.plot(x, y, 'o',ms=18,  markerfacecolor = 'b', markeredgecolor =None)

    for label, xx, yy in zip(labels, x, y):
        plt.annotate(label, xy=(xx,yy), textcoords="offset points", ha = 'center',
        va = 'center', color='w'
        
        )
    plt.show()

if __name__ == '__main__': # If this file is run directly
    main()

    
