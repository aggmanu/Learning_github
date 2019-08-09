import numpy as np
import itertools as it
import time
import scipy
import math
import matplotlib.pyplot as plt

import Class_persistence as Cp

#
#locations = np.random.uniform(0, 1, size=(5,2))


#new = Cp.filtration(locations)
#new.generate_simplices()

#new.evaluate_persistence_homology()
#print(new.homology_gps)



def get_b_ad(cell_types):

  beta_cells = np.where(cell_types== 2)[0].tolist() 
  alpha_cells = np.where(cell_types== 1)[0].tolist() 
  delta_cells = np.where(cell_types== 3)[0].tolist()
  ad_cells = alpha_cells + delta_cells

  return beta_cells, ad_cells


stage =  3
islet  = 12170

#stage =  1
#islet  = 240846

#stage =  2
#islet  = 203210

#stage = 2
#islet = 80368

thresh = 50
#thresh = math.inf
beta_thresh = thresh
print('Loading Stage:', stage)

#all_islets =  np.load('./Data/stage'+str(stage)+'.npy')
start = time.time()
all_islets =  np.loadtxt('./Data/stage'+str(stage)+'.tsv', delimiter='\t')
finish = time.time()
print('Time taken to load all the islets of this stage:', finish-start)
input('Press key to continue.')

islet_idxs  = np.where(all_islets[:,0] == islet)[0]

islet_info = all_islets[islet_idxs, :]
islet_info_orig = all_islets[islet_idxs, :]

beta_cells = np.where(islet_info[:,4] == 2)[0].tolist() 
alpha_cells = np.where(islet_info[:,4] == 1)[0].tolist() 
delta_cells = np.where(islet_info[:,4] == 3)[0].tolist()
ad_cells = alpha_cells + delta_cells



#ad_cells = np.where(islet_info[:,4] == 1 or islet_info[:,4] = 3 )[0] 
locations = islet_info[:, 2:4]
print('Loaded information for islet',  islet)
print('Calculating pairwise dist between all of the cells...')
start = time.time()
all_pair_dist = scipy.spatial.distance.cdist(locations, locations)
finish = time.time()

print('Time taken to calculate pairwise dist between all cells:', finish-start)

def check_intersects_beta(a1_idx, a2_idx):
  
      intersects_beta = False

      ad1_ad2 = all_pair_dist[a1_idx, a2_idx]
      dist_thresh = 2*ad1_ad2
  
      for beta_pair in it.combinations(beta_cells, 2):
          b1_idx = beta_pair[0]
          b2_idx = beta_pair[1]

          b1_ad1 = all_pair_dist[b1_idx, a1_idx]
          b1_ad2 = all_pair_dist[b1_idx, a2_idx]

          if b1_ad1 > dist_thresh and b1_ad1 > dist_thresh:
            continue

          b2_ad1 = all_pair_dist[b2_idx, a1_idx]
          b2_ad2 = all_pair_dist[b2_idx, a2_idx]

          if b2_ad1 > dist_thresh and b2_ad1 > dist_thresh:
            continue

          b1_b2 = all_pair_dist[b1_idx, b2_idx]

          if b1_b2 > ad1_ad2:
            continue


          A1 = locations[a1_idx]
          A2 = locations[a2_idx]
          
          B1 = locations[b1_idx]
          B2 = locations[b2_idx]

          C1 = A2-A1
          C2 = B1-B2

          mat = np.array([[C1[0], C2[0]]\
                          ,[C1[1], C2[1]]])
  
          C3 = B1 - A1
  
          #try:
          try:
            pars = np.linalg.solve(mat, C3)
            if max(pars) <= 1 and min(pars) >=0:
                intersects_beta = True
                break
          except:
            #Parallel lines
            continue

      return intersects_beta
              
def check_intersects_ad(b1_idx, b2_idx):
  
      intersects_ad = False

      b1_b2 = all_pair_dist[b1_idx, b2_idx]
      dist_thresh = 2*b1_b2
  
      for ad_pair in it.combinations(ad_cells, 2):
          ad1_idx = ad_pair[0]
          ad2_idx = ad_pair[1]

          ad1_b1 = all_pair_dist[b1_idx, ad1_idx]
          ad1_b2 = all_pair_dist[b2_idx, ad1_idx]

          if ad1_b1 > dist_thresh and ad1_b2 > dist_thresh:
            continue

          ad2_b1 = all_pair_dist[b1_idx, ad2_idx]
          ad2_b2 = all_pair_dist[b2_idx, ad2_idx]

          if ad2_b1 > dist_thresh and ad2_b2 > dist_thresh:
            continue
  
          ad1_ad2 = all_pair_dist[ad1_idx, ad2_idx]

          if ad1_ad2 > b1_b2:
            continue
  
          A1 = locations[ad1_idx]
          A2 = locations[ad2_idx]
          
          B1 = locations[b1_idx]
          B2 = locations[b2_idx]
  
          C1 = A2-A1
          C2 = B1-B2
          mat = np.array([[C1[0], C2[0]]\
                          ,[C1[1], C2[1]]])
  
          C3 = B1 - A1
  
          try:
            pars = np.linalg.solve(mat, C3)
  
            if max(pars) <= 1 and min(pars) >=0:
                intersects_ad = True
                break
          except:
            #Parallel lines
            continue
  

      return intersects_ad


def simplex1_check(v1, v2):
  
  if all_pair_dist[v1, v2] > thresh:
    return False
  elif (v1 in beta_cells) and (v2 in beta_cells):
    return (not check_intersects_ad(v1, v2))
  #elif (v1 in ad_cells) and (v2 in ad_cells): 
  #  return (not check_intersects_beta(v1, v2))
  else:
    return False

start = time.time()
original = Cp.filtration(locations, simplex1_check = simplex1_check)

original.generate_simplices()
original.evaluate_persistence_homology()

#new.evaluate_persistence_homology()
hom = original.homology_gps[0]

deaths_orig_beta = []
deaths_orig_ad = []
deaths_orig = []

for pt in hom:
  death = pt[0][1]
  if math.isinf(death):
    deaths_orig.append(thresh)
  else:
    deaths_orig.append(death)
finish = time.time()
print(finish - start)
input('w')


#plt.hist(deaths, bins = 100, alpha = 0.5)
#plt.show()

wmin = math.inf
permute = islet_info[:, 4]

#for i in range(10):
while wmin > 0.05:
  # Permute
  
  islet_info[:, 4] = np.random.permutation(islet_info[:, 4])

  beta_cells, ad_cells = get_b_ad(islet_info[:, 4])
  
  #beta_cells = np.where(islet_info[:,4] == 2)[0].tolist() 
  #alpha_cells = np.where(islet_info[:,4] == 1)[0].tolist() 
  #delta_cells = np.where(islet_info[:,4] == 3)[0].tolist()
  #ad_cells = alpha_cells + delta_cells
  
  new = Cp.filtration(locations, simplex1_check = simplex1_check)
  
  new.generate_simplices()
  new.evaluate_persistence_homology()
  
  #new.evaluate_persistence_homology()
  hom = new.homology_gps[0]
  
  deaths = []
  
  for pt in hom:
    death = pt[0][1]
    if math.isinf(death):
      deaths.append(thresh)
    else:
      deaths.append(death)
  
  wdist = scipy.stats.wasserstein_distance(deaths_orig, deaths)
  print(wdist)
  if wdist < wmin:
    wmin = wdist
    permute = islet_info[:, 4]
  

fig, (ax1, ax2) = plt.subplots(1, 2)

beta_cells, ad_cells = get_b_ad(islet_info_orig[:, 4])
print(beta_cells)
print(ad_cells)
ax1.scatter(locations[beta_cells][:,0], locations[beta_cells][:,1], color='g')
ax1.scatter(locations[ad_cells][:,0], locations[ad_cells][:,1], color='r')


beta_cells, ad_cells = get_b_ad(permute)
print(beta_cells)
print(ad_cells)
ax2.scatter(locations[beta_cells][:,0], locations[beta_cells][:,1], color='g')
ax2.scatter(locations[ad_cells][:,0], locations[ad_cells][:,1], color='r')

plt.show()



