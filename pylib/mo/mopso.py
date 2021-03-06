# mopso.py

import random, math, copy
from .individual import *
import config

class mopso:
    '''
    Implementation of MOPSO algorithm.
    '''
    def __init__(self, archive_size=100, mutation_rate=0.1, c1=1.0, c2=1.0, w= 0.4):
        '''
        Constructor. 
        Parameters: archive_size, mutation rate (default value 10%) and crossover rate (default value 100%). 
        '''
        random.seed()
        self.P = []
        for i in range(config.population_size):
          self.P.append(Individual()) 
        self.output_population('init_pop')

        self.archive = []
        self.archive_size = archive_size
        self.mutation_rate = mutation_rate
        self.w = w
        self.c1 = c1
        self.c2 = c2       

    def output_population(self, name):
        obj_file = name + "_obj.csv"
        var_file = name + "_var.csv"
        obj_out = open(obj_file, 'w')
        var_out = open(var_file, 'w')
        for a in self.P:
          if a.constraint_violation >= 0.0:
            for o in a.objectives:
              obj_out.write("" + str(o) + ", ")
            obj_out.write("\n")
            var_out.write("" + str(a.variables) + "\n")
        obj_out.close()
        var_out.close()

    def output_archive(self, name):
        obj_file = name + "_obj.csv"
        var_file = name + "_var.csv"
        obj_out = open(obj_file, 'w')
        var_out = open(var_file, 'w')
        for a in self.archive:
          for o in a.objectives:
            obj_out.write("" + str(o) + ", ")
          obj_out.write("" + str(a.constraint_violation) + ", " + str(a.distance)+ "\n")
          for v in a.variables:
            var_out.write("" + str(v) + ", ")
          var_out.write("\n")
        obj_out.close()
        var_out.close()
          
    def run(self, num_generation):
      if config.eval_func == '':
	print( 'Error from MOPSO: No evaluation function is provided!')
	exit()
      self.update_archive()
      for g in range(num_generation):
        #print "============== start gen " + str(g) + " ============"
        self.update_population()
        mut_gen = num_generation*self.mutation_rate
        if g < int(mut_gen):
          for p in self.P:
            p.mutate(math.pow(1.0-float(g)/mut_gen, 1.5)) 
        for p in self.P:
          p.evaluate()
          p.update_personal_best()
        self.update_archive()
        #print "============== done gen " + str(g) + " ============"
#        self.output_population('pop_' + str(g))
        #print out only the front
        self.output_archive('ar_' + str(g))

    def check_archive_constraints_negative(self):
      for a in self.archive:
        if a.constraint_violation < 0.0:
          continue
        else:
          return False
      return True

    def update_population(self):
      gbest_index_range = int(len(self.archive)*0.3)
      for p in self.P:
        if self.check_archive_constraints_negative() or gbest_index_range == 0:
          gbest = random.choice(self.archive)
        else:
          gbest = random.choice(self.archive[0:gbest_index_range])
#        print "Update_population, gbest =  " + str(gbest.objectives) + str(gbest.variables) + str(gbest.constraint_violation)
        for i in range(p.num_variables):
          p.velocities[i] = self.w*p.velocities[i] + \
            self.c1*random.random()*(p.variables_personal_best[i]-p.variables[i]) + \
            self.c2*random.random()*(gbest.variables[i]-p.variables[i])
          p.variables[i] += p.velocities[i]
        p.maintain_variable_range()
         
    def update_archive(self):
      #print "UPDATE ARCHIVE!!!!!"
      for p in self.P:
        #print "---- pop item " + str(p.objectives) + str(p.constraint_violation)
        insert_flag = True
        if len(self.archive) == 0:
          self.archive.append(copy.deepcopy(p))
          #print "--archive is empty, add it to archive"
        else:
          remove_list = []
          for a in self.archive:
            #print "---- archive item " + str(a.objectives) + str(a.constraint_violation)
            if p.close_to(a):
#              if p == a:
              insert_flag = False
              #print "!!! archive item a and p are equal!!!"
              break 
            if a >> p:
              insert_flag = False
              #print "!!! archive item dominates p !!!"
              break 
            elif p >> a:
#              print "!!! archive item " + str(a.objectives) + str(a.variables) + str(a.constraint_violation)
#              print "!!!pop item " + str(p.objectives) + str(p.variables) + str(p.constraint_violation)
              #print "!!! archive itme  failed to dominate p, remove a from archive"
              remove_list.append(a) 
          for r in remove_list:
            self.archive.remove(r)

          if insert_flag:
            if len(self.archive) < self.archive_size:
              self.archive.append(copy.deepcopy(p))
#              print "***add item into archive: " + str(p.objectives) + str(p.variables) + str(p.constraint_violation)
            else:
              self.crowding_distance_assignment(self.archive)
              self.sort_crowding(self.archive)
#              print"!!! after cd sorting:"
#              for a in self.archive:
#                print "***: " + str(a.objectives) + str(a.variables) + str(a.constraint_violation) + ", " + str(a.distance)
              range_start = int(len(self.archive)*0.9);
              if range_start == self.archive_size-1:
                range_start = 0
              replace_index = random.choice(range(range_start, self.archive_size-1))
              a = self.archive[replace_index]
#              print "!!!archive full start = " + str(range_start)
#              print "!!!archive full! replace: " + str(a.objectives) + str(a.variables) + str(a.constraint_violation)
#              print "!!!with : " + str(p.objectives) + str(p.variables) + str(p.constraint_violation)
              self.archive[replace_index] = copy.deepcopy(p)

      self.crowding_distance_assignment(self.archive)
      self.sort_crowding(self.archive)
#      for a in self.archive:
#        print "archive: " + str(a.objectives) + str(a.variables) + str(a.constraint_violation)
#      print "========================================="
             
    def sort_objective(self, P, obj_idx):
        for i in range(len(P) - 1, -1, -1):
            for j in range(1, i + 1):
                s1 = P[j - 1]
                s2 = P[j]
                
                if s1.objectives[obj_idx] > s2.objectives[obj_idx]:
                    P[j - 1] = s2
                    P[j] = s1
                    
    def sort_crowding(self, P):
        for i in range(len(P) - 1, -1, -1):
            for j in range(1, i + 1):
                s1 = P[j - 1]
                s2 = P[j]
                
                if crowded_comparison(s1, s2) < 0:
                    P[j - 1] = s2
                    P[j] = s1
                
    def crowding_distance_assignment(self, front):
        '''
        Assign a crowding distance for each solution in the front. 
        '''
        for p in front:
            p.distance = 0
        
        for obj_index in range(front[0].num_objectives):
            self.sort_objective(front, obj_index)
            if not (front[0].objectives[obj_index] == float('inf') and front[len(front)-1].objectives[obj_index] == float('inf')): 
              front[0].distance = float('inf')
              front[len(front) - 1].distance = float('inf')
              obj_min = front[0].objectives[obj_index]
              obj_max = front[len(front)-1].objectives[obj_index]
              for i in range(1, len(front)-1):
                  if obj_max - obj_min == 0 :
                    front[i].distance += 0.0
                  else:
                    front[i].distance += (front[i+1].objectives[obj_index] - front[i-1].objectives[obj_index])/(obj_max - obj_min)
            elif front[0].objectives[obj_index] == float('inf'):
              for i in range(1, len(front)-1):
                front[i].distance = 0.0
            elif front[len(front)-1].objectives[obj_index] == float('inf'):
              front_index = len(front)-2
              while front_index > 0 and front[front_index].objectives[obj_index] == float('inf'):
                front_index -= 1
              front[0].distance = float('inf')
              front[front_index].distance = float('inf')
              obj_min = front[0].objectives[obj_index]
              obj_max = front[front_index].objectives[obj_index]
              for i in range(1, front_index):
                  if obj_max - obj_min == 0 :
                    front[i].distance += 0.0
                  else:
                    front[i].distance += (front[i+1].objectives[obj_index] - front[i-1].objectives[obj_index])/(obj_max - obj_min)
              

    def fast_nondominated_sort(self, P):
        '''
        Discover Pareto fronts in P, based on non-domination criterion. 
        '''
        fronts = {}
        
        S = {}
        n = {}
        for s in P:
            S[s] = [] # set donimated by s
            n[s] = 0  # number of elements that dominate s
            
        fronts[1] = []
        
        for p in P:
            for q in P:
                if p == q:
                    continue
                
                if p >> q:
                    S[p].append(q)
                
                elif p << q:
                    n[p] += 1
            
            if n[p] == 0:
                fronts[1].append(p)
        
        i = 1
        
        while len(fronts[i]) != 0:
            next_front = []
            
            for r in fronts[i]:
                for s in S[r]:
                    n[s] -= 1
                    if n[s] == 0:
                        next_front.append(s)
            
            i += 1
            fronts[i] = next_front
                    
        return fronts
 
