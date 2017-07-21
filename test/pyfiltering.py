import numpy as np
import sys

try:
    xrange
except NameError:
    xrange = range


class interval():
    def __init__(self,start,end,a,x,pval):
        self.start = start
        self.end = end
        self.a = a
        self.x = x
        self.pval = pval

    def __lt__(self,temp):
        return self.start < temp.start

    def length(self):
        return self.end - self.start + 1

class cluster():
    def __init__(self,start,end,most_sig_int):
        self.start = start
        self.end = end
        self.most_sig_int = most_sig_int
        self.n_intervals = 1 
    
    def update(self,new_int):
        if new_int.start <= self.end:
            self.end = max(new_int.end,self.end)
            self.n_intervals += 1
            if new_int.pval < self.most_sig_int.pval:
                self.most_sig_int = new_int
            elif new_int.pval == self.most_sig_int.pval and new_int.length() < self.most_sig_int.length():
                self.most_sig_int = new_int
            return 1
        else:
            return 0

if __name__ in "__main__":
    T = np.loadtxt(sys.argv[1],delimiter=',',skiprows=1,ndmin=2)
    if T.shape[0]==0:
        import sys
        sys.exit(0)
    intervals = list()
    for i in xrange(T.shape[0]):
        intervals.append(interval(T[i,1],T[i,1]+T[i,0]-1,T[i,2],T[i,3],T[i,4]))
    intervals.sort()

    clusters = list()
    tmp_cluster = cluster(intervals[0].start,intervals[0].end,intervals[0])
    for i in xrange(1,len(intervals)):
        out = tmp_cluster.update(intervals[i])
        if out==0:
            clusters.append(tmp_cluster)
            tmp_cluster = cluster(intervals[i].start,intervals[i].end,intervals[i])
    clusters.append(tmp_cluster)
    
    Tc = np.zeros((len(clusters),8))
    for i in xrange(len(clusters)):
        Tc[i,0] = clusters[i].most_sig_int.start
        Tc[i,1] = clusters[i].most_sig_int.end
        Tc[i,2] = clusters[i].most_sig_int.a
        Tc[i,3] = clusters[i].most_sig_int.x
        Tc[i,4] = clusters[i].most_sig_int.pval
        Tc[i,5] = clusters[i].n_intervals
        Tc[i,6] = clusters[i].start
        Tc[i,7] = clusters[i].end

    np.savetxt(sys.argv[2],Tc,fmt='%d %d %d %d %e %d %d %d')

