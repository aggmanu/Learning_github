import numpy as np
from sklearn.metrics import pairwise_distances
import itertools as it
import math
import time
from sys import getsizeof

class filtration:

  def __init__(self, locations, simplex1_check=None, simplex2_check=None, simplex3_check=None):

    self.simplex_to_order = dict()
    self.order_to_simplex = dict()
    self.locations = locations
    self.simplex1_check = simplex1_check
    self.simplex2_check = simplex2_check
    self.simplex3_check = simplex3_check
    self.homology_gps = []
    self.order = -1

    self.parameters = []
    self.boundaries = []
    self.col_operations = []
    self.is_last_boundary = []
    self.homology_gps = [[], [], []]



  def add_simplex(self, simplex, parameter, boundaries):

    self.boundaries.append(boundaries)
    self.simplex_to_order.update({simplex:self.order})
    self.order_to_simplex.update({self.order:simplex})
    self.parameters.append(parameter)
    self.is_last_boundary.append([])
    self.col_operations.append([self.order])

    if boundaries:
      last_boundary_order = max(boundaries)
      self.is_last_boundary[last_boundary_order].append(self.order)

  def generate_simplices(self):

    data = pairwise_distances(self.locations)
    n_rows, n_cols = data.shape

    edges = []

    parameter = 0

    # Add 0-simplices
    for v in range(n_rows):
      self.order = self.order + 1
      boundary = []
      self.add_simplex(frozenset({v,}), parameter, boundary)

    # Make a list of 1-simplices

    for pair in it.combinations(list(range(n_rows)), 2):
        # Add here any restrictions on the edges

        if math.isinf(data[pair[0], pair[1]]):
            continue

        if self.simplex1_check:
            if not self.simplex1_check(pair[0], pair[1]):
                continue

        edges.append([(pair[0], pair[1]),data[pair[0], pair[1]]])

            
    # Now sort the list onthe basis of the distance
    edges.sort(key=lambda x: x[1])
    
    # The 'boundary' for 0-simplices will store its neighbors
    n_edges = len(edges)


    count = 0

    # Iterate through the edges sorted by length
    for edge, parameter in edges:
        count = count + 1
        #print(count, n_edges, end='\r')

        simplex = frozenset(edge)
        v0_set = frozenset({edge[0],})
        v1_set = frozenset({edge[1],})

        v0_order = self.simplex_to_order[v0_set]
        v1_order = self.simplex_to_order[v1_set]

        # Add 1-simplices
        self.order = self.order + 1
        boundaries = [v0_order, v1_order]

        self.add_simplex(simplex, parameter, boundaries)

        self.boundaries[v0_order].append(v1_order)

        self.boundaries[v1_order].append(v0_order)

        ## Find common neighbors

        #common_neighbors = []
        #for v in self.boundaries[v0_order]:
        #  if v in self.boundaries[v1_order]:
        #    common_neighbors.append(v)

        ## Add 2-simplices
        #for v in common_neighbors:

        #    if self.simplex2_check:
        #        if not self.simplex2_check(edge[0], edge[1], v):
        #            continue


        #    #print(order, simplex, end='\r')
        #    self.order = self.order + 1
        #    simplex = frozenset({edge[0], edge[1], v})
        #    o1 = self.simplex_to_order[frozenset({edge[0], edge[1]})]
        #    o2 = self.simplex_to_order[frozenset({edge[1], v})]
        #    o3 = self.simplex_to_order[frozenset({edge[0], v})]
        #    boundaries = [o1, o2, o3]
        #    self.add_simplex(simplex, parameter, boundaries)

            #for v2 in self.boundaries[v]:
            #  if v2 in common_neighbors:
            #      s1 = frozenset({v, v2, edge[0]})
            #      if self.simplex2_check:
            #          if not self.simplex2_check(v, v2, edge[0]):
            #              continue
            #      s2 = frozenset({v2, edge[0], edge[1]})
            #      if self.simplex2_check:
            #          if not self.simplex2_check(v2, edge[0], edge[1]):
            #              continue
            #      s4 = frozenset({v, v2, edge[1]})
            #      if self.simplex2_check:
            #          if not self.simplex2_check(v, v2, edge[1]):
            #              continue

            #      o1 = self.simplex_to_order[frozenset({v, v2, edge[0]})]
            #      o2 = self.simplex_to_order[frozenset({v2, edge[0], edge[1]})]
            #      o4 = self.simplex_to_order[frozenset({v, v2, edge[1]})]

           
        #for pair in it.combinations(common_neighbors, 2):
        #   v = pair[0]
        #   v2 = pair[1]
        #   simplex = frozenset({v, v2, edge[0], edge[1]})
        #   if simplex not in self.simplex_to_order:
        #     print(order, simplex, end='\r')
        #     self.order = self.order + 1
        #     o1 = self.simplex_to_order[frozenset({v, v2, edge[0]})]
        #     o2 = self.simplex_to_order[frozenset({v2, edge[0], edge[1]})]
        #     o3 = self.simplex_to_order[frozenset({v, edge[0], edge[1]})]
        #     o4 = self.simplex_to_order[frozenset({v, v2, edge[1]})]
        #     boundaries = [o1, o2, o3, o4]
        #     self.add_simplex(simplex, parameter, boundaries)
           

         ## Add 3-simplices
         #for v in common_neighbors:
         #   for v2 in self.boundaries[v]:
         #     if v2 in common_neighbors:
         #       simplex = frozenset({v, v2, edge[0], edge[1]})
         #       if simplex not in self.simplex_to_order:
         #         print(order, simplex, end='\r')
         #         order = order + 1
         #         o1 = self.simplex_to_order[frozenset({v, v2, edge[0]})]
         #         o2 = self.simplex_to_order[frozenset({v2, edge[0], edge[1]})]
         #         o3 = self.simplex_to_order[frozenset({v, edge[0], edge[1]})]
         #         o4 = self.simplex_to_order[frozenset({v, v2, edge[1]})]
         #         boundaries = [o1, o2, o3, o4]
         #         self.add_simplex(simplex, order, parameter, boundaries)

  def reduce_filtration(self, simplex_order):
  
          last_boundary_order = max(self.boundaries[simplex_order])
  
          first_last_boundary_order = min(self.is_last_boundary[last_boundary_order])
  
          while first_last_boundary_order < simplex_order:


            self.boundaries[simplex_order] = list(\
                            set(self.boundaries[simplex_order])\
                              .symmetric_difference(\
                                set(self.boundaries[first_last_boundary_order])))

            # Record the column operation
            self.col_operations[simplex_order] = self.col_operations[simplex_order] +\
                                              self.col_operations[first_last_boundary_order]

            # This boundary is not the last boundary of the simplex anymore
            self.is_last_boundary[last_boundary_order].remove(simplex_order)


            if self.boundaries[simplex_order]:

              # Find the new max
              last_boundary_order = max(self.boundaries[simplex_order])

              # Update the is_last_boundary of this new max
              self.is_last_boundary[last_boundary_order].append(simplex_order)

              # Find the new first_last_boundary_order
              first_last_boundary_order = min(self.is_last_boundary[last_boundary_order])

            else:
              break



  
  def evaluate_persistence_homology(self):


    for i in range(self.order+1):
        simplex = self.order_to_simplex[i]
        if len(simplex) < 2:
            continue
        self.reduce_filtration(i) 

    # The aim right now is to get only H0
    for i in range(self.order+1):
        simplex = self.order_to_simplex[i]
        if len(simplex) > 1:
          break

        if not self.is_last_boundary[i]:
          birth_order = i
          birth = self.parameters[i]
          death = math.inf
          death_order = math.inf
          self.homology_gps[0].append([(birth, death), (birth_order, death_order), simplex])
        elif len(self.is_last_boundary[i]) == 1:
          birth_order = i
          birth = self.parameters[i]
          death_order = self.is_last_boundary[i][0]
          death = self.parameters[death_order]
          self.homology_gps[0].append([(birth, death), (birth_order, death_order), simplex])
        else:
          print('Unexpected')



    #for i in is

    #for i in range(self.order+1):
    #    simplex = self.order_to_simplex[i]
    #    if len(simplex) < 2:
    #        continue

    #    if not self.boundaries[i]:
    #        if not self.is_last_boundary[i]:
    #            birth = self.parameters[i]
    #            birth_order = i 
    #            death = math.inf
    #            dim = len(simplex)
    #            self.homology_gps[dim-1].append([(birth, death), (birth_order, death_order), self.col_operations[i]])
    #            print(simplex)
    #    else:
    #        death = self.parameters[i]
    #        death_order = i
    #        birth_order = max(self.boundaries[i])
    #        birth = self.parameters[birth_order]
    #        birth_simplex = self.order_to_simplex[birth_order]
    #        dim = len(birth_simplex)

    #        if death != birth:
    #            self.homology_gps[dim-1].append([(birth, death), (birth_order, death_order), self.col_operations[i]])




