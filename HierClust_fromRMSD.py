from scipy import cluster
import numpy as np
import sys,os
from itertools import combinations

# Hierarchical clustering using scipy.cluster.hierarchical
# Input: condensed distance matrix upper triangular, including 0.0 diagonal (argv[1]);
# decoy list to map it back to the clustered points (argv[2]); and cophenetic distance cutoff
# for flat cluster generation.

def Get_n_clusters_list_of_lists(distance_file,decoy_file,nu):
    print "Starting"
    distancesfromfile = np.fromfile(distance_file,sep=" ")
    decoys = open(decoy_file,'r').readlines()
    distances = [n for n in distancesfromfile]
    print "Done reading distances"
    # Do the clustering                                                                                                                                                                                                                                                      
    linkage_mtx = cluster.hierarchy.linkage(distances,method='complete')
    print "Done clustering"
    # Generation of independent clusters, cutoff in this context if the highest RMSD that two members
    # of different clusters can have before being considered in different clusters, it's a cutoff
    # for the linkage criterium (complete in this case), it's very conservative:
    flat_clusters = cluster.hierarchy.fcluster(linkage_mtx,nu,criterion='maxclust')
    print "Done making flat clusters"
    # Get cluster centers, i.e., structures with lowest average RMSD to other in same cluster
    # Start a cluster dictionary:
    str_clust_tups = []
    clust_dict = {}
    for i,n in enumerate(decoys):
        str_clust_tups.append( (n,flat_clusters[i]) )
        if flat_clusters[i] not in clust_dict.keys():
            clust_dict[flat_clusters[i]] = []
            clust_dict[flat_clusters[i]].append(n)
        else:
            clust_dict[flat_clusters[i]].append(n)
    # Then expand the distance matrix for easy access
    n_dec = len(decoys)
    storage = np.zeros((n_dec,n_dec)) # make a big matrix full of zeros 
    i_storage = np.triu_indices(n_dec,1) # get its upper triangle indices excluding diagonal
    storage[i_storage] = distances # fill the upper triangle with the distances that were condensed
    print "Done creating cluster dictionary and distance matrix"
    output = []
    print "clust_dict keys: %i"%len(clust_dict.keys())
    print clust_dict
    for key in clust_dict.keys():
        cluster_list = clust_dict[key]
        output.append(cluster_list)
        
    return output

def Get_clusters_by_RMSD_cutoff(distance_file,decoy_file,cutoff,output_dir):
    print "Starting"
    distancesfromfile = np.fromfile(distance_file,sep=" ")
    decoys = open(decoy_file,'r').readlines()
    distances = [n for n in distancesfromfile]
    print "Done reading distances"
    # Do the clustering
    linkage_mtx = cluster.hierarchy.linkage(distances,method='complete')
    print "Done clustering"
    # Generation of independent clusters, cutoff in this context if the highest RMSD that two members
    # of different clusters can have before being considered in different clusters, it's a cutoff
    # for the linkage criterium (complete in this case), it's very conservative:
    flat_clusters = cluster.hierarchy.fcluster(linkage_mtx,cutoff,criterion='distance')
    print "Done making flat clusters"
    # Get cluster centers, i.e., structures with lowest average RMSD to other in same cluster
    # Start a cluster dictionary:
    str_clust_tups = []
    clust_dict = {}
    for i,n in enumerate(decoys):
        #print "%d %s"%( flat_clusters[i], n)
        str_clust_tups.append( (n,flat_clusters[i]) )
        if flat_clusters[i] not in clust_dict.keys():
            clust_dict[flat_clusters[i]] = []
            clust_dict[flat_clusters[i]].append(n)
        else:
            clust_dict[flat_clusters[i]].append(n)
    
    # Then expand the distance matrix for easy access
    n_dec = len(decoys)
    storage = np.zeros((n_dec,n_dec)) # make a big matrix full of zeros
    i_storage = np.triu_indices(n_dec,1) # get its upper triangle indices excluding diagonal
    storage[i_storage] = distances # fill the upper triangle with the distances that were condensed
    print "Done creating cluster dictionary and distance matrix"
    # Now get pairwise distances for all decoys in cluster and keep the min
    cluster_centers = {}
    for key in clust_dict.keys(): # for each cluster
        places = []
        #print "Printing clust_dict[key] for cluster %d"%key
        #print clust_dict[key]
        if len(clust_dict[key]) > 1:
            for decoy in clust_dict[key]: # Get the decoys in it and their indices in the original array
                places.append(decoys.index(decoy))
            combs = [sorted(i) for i in combinations(places,2)] # Make combinations of them to use as indexes. Each combination is sorted to go for the upper triangular
            av_RMSD = []
            for place in places: # for each decoy compute the average distance to the others
                vals = []
                for comb in combs:
                    if place in comb:
                        vals.append(storage[comb[0]][comb[1]])
                av_RMSD.append((place,np.average(vals)))
            #print av_RMSD
            cluster_centers[key] = decoys[min(av_RMSD,key=lambda x: x[1])[0]] # get the minimum and save it
            #print "The cluster center is %s"%(decoys[min(av_RMSD)[0]])
        else:
            cluster_centers[key] = clust_dict[key][0]
    # Print stuff
    sorted_clust_str = sorted(str_clust_tups, key=lambda x: x[1])
    out_clusters = open(output_dir+'/clusters.list','w')
    out_centers = open(output_dir+'/centers.list','w')
    for i in sorted_clust_str:
        out_clusters.write("CLUST: "+i[0][:-1]+" "+str(i[1])+'\n')
        print "CLUST: "+i[0][:-1]+" "+str(i[1])
    for key in clust_dict.keys():
        out_centers.write("CENT: "+cluster_centers[key][:-1]+" "+str(key)+'\n')
        print "CENT: "+cluster_centers[key][:-1]+" "+str(key)
    out_clusters.close()
    out_centers.close()

if __name__ == '__main__':
    Get_clusters_by_RMSD_cutoff(sys.argv[1],sys.argv[2],float(sys.argv[3]),'./')

'''
print "Starting"
distancesfromfile = np.fromfile(sys.argv[1],sep=" ")
decoy_file = open(sys.argv[2],'r')
decoys = decoy_file.readlines()
distances = [n for n in distancesfromfile]
print "Done reading distances"
# Do the clustering
linkage_mtx = cluster.hierarchy.linkage(distances,method='complete')
print "Done clustering"
cutoff = float(sys.argv[3])
# Generation of independent clusters, cutoff in this context if the highest RMSD that two members
# of different clusters can have before being considered in different clusters, it's a cutoff
# for the linkage criterium (complete in this case), it's very conservative:
flat_clusters = cluster.hierarchy.fcluster(linkage_mtx,cutoff,criterion='distance')
print "Done making flat clusters"
# Get cluster centers, i.e., structures with lowest average RMSD to other in same cluster
 # Start a cluster dictionary:
str_clust_tups = []
clust_dict = {}
for i,n in enumerate(decoys):
    str_clust_tups.append( (n,flat_clusters[i]) )
    if flat_clusters[i] not in clust_dict.keys():
        clust_dict[flat_clusters[i]] = []
        clust_dict[flat_clusters[i]].append(n)
    else:
        clust_dict[flat_clusters[i]].append(n)

 # Then expand the distance matrix for easy access
n_dec = len(decoys)
storage = np.zeros((n_dec,n_dec)) # make a big matrix full of zeros
i_storage = np.triu_indices(n_dec,1) # get its upper triangle indices excluding diagonal
storage[i_storage] = distances # fill the upper triangle with the distances that were condensed
print "Done creating cluster dictionary and distance matrix"
 # Now get pairwise distances for all decoys in cluster and keep the min
cluster_centers = {}

for key in clust_dict.keys(): # for each cluster
    places = []
    if len(clust_dict[key]) > 1:
        for decoy in clust_dict[key]: # Get the decoys in it and their indices in the original array
            places.append(decoys.index(decoy))
        combs = [sorted(i) for i in combinations(places,2)] # Make combinations of them to use as indexes. Each combination is sorted to go for the upper triangular
        av_RMSD = []
        for place in places: # for each decoy compute the average distance to the others
            vals = []
            for comb in combs:
                if place in comb:
                    vals.append(storage[comb[0]][comb[1]])
            av_RMSD.append((place,np.average(vals)))
            cluster_centers[key] = decoys[min(av_RMSD)[0]] # get the minimum and save it
    else:
        cluster_centers[key] = clust_dict[key][0]

# Print stuff

sorted_clust_str = sorted(str_clust_tups, key=lambda x: x[1])
for i in sorted_clust_str:
    print "CLUST: "+i[0][:-1]+" "+str(i[1])

for key in clust_dict.keys():
    print "CENT: "+cluster_centers[key][:-1]+" "+str(key)
'''
