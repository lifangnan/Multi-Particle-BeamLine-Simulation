# individual.py

import sys, random
import config

def equal_to(o1, o2):
    deadband = abs(o1)*0.0001
    if abs(o1-o2) < deadband:
      return True
    else:
      return False

def less_than(o1, o2):
    if (not equal_to(o1, o2)) and o1 < o2:
      return True
    else:
      return False

def larger_than(o1, o2):
    if (not equal_to(o1, o2)) and o1 > o2:
      return True
    else:
      return False


class Solution:
    '''
    Abstract solution. To be implemented.
    '''
    
    def __init__(self):
        '''
        Constructor. Parameters: number of objectives. 
        '''
        self.num_objectives = config.num_objectives
        self.objectives = []
        for _ in range(config.num_objectives):
            self.objectives.append(None)
        self.num_constraints = config.num_constraints
        self.constraints = []
        for _ in range(config.num_constraints):
            self.constraints.append(0.0)
        self.constraint_violation = 0.0;
        self.num_variables = config.num_variables
        self.variables = []
        self.variable_max = config.variable_max
        self.variable_min = config.variable_min
        self.velocities = []
        for _ in range(config.num_variables):
          self.variables.append(None)
#          self.variable_max.append(float('inf'))    
#          self.variable_min.append(float('-inf'))    
          self.velocities.append(0.0)
        self.objectives_personal_best = []
        self.variables_personal_best = []
        self.constraint_violation_personal_best = 0.0 
        self.rank = sys.maxint
        self.distance = 0.0
        
    def evaluate(self):
        '''
        Evaluate solution, update objectives values.
        '''
        raise NotImplementedError("Solution class have to be implemented.")
    
    def mutate(self):
        '''
        Mutation operator.
        '''
        raise NotImplementedError("Solution class have to be implemented.")
    def __rshift__(self, other):
        '''
        True if this solution dominates the other (">>" operator).
        '''
        if self.constraint_violation < 0.0 and other.constraint_violation < 0.0:
          if self.constraint_violation > other.constraint_violation:
            return True
          else:
            return False 
        elif self.constraint_violation >= 0.0 and other.constraint_violation < 0.0:
          return True
        elif self.constraint_violation < 0.0 and other.constraint_violation >= 0.0:
          return False
        else:
          dominates = False
          for i in range(len(self.objectives)):
              if larger_than(self.objectives[i], other.objectives[i]):
                  return False
              elif less_than(self.objectives[i], other.objectives[i]):
                  dominates = True
          return dominates
        
    def __lshift__(self, other):
        '''
        True if this solution is dominated by the other ("<<" operator).
        '''
        return other >> self

    def close_to(self, other):
        if self.constraint_violation >= 0.0 and other.constraint_violation >= 0.0:
          for i in range(len(self.objectives)):
            if not equal_to(self.objectives[i], other.objectives[i]):
              return False
          return True
        else:
          return False
        
    def maintain_variable_range(self):
        for i in range(self.num_variables):
          if self.variables[i] < self.variable_min[i]:
            self.variables[i] = self.variable_min[i]
            self.velocities[i] *= -1.0
          if self.variables[i] > self.variable_max[i]:
            self.variables[i] = self.variable_max[i]
            self.velocities[i] *= -1.0
            
    def update_personal_best(self):
      pbest = Solution()
      pbest.variables = self.variables_personal_best
      pbest.objectives = self.objectives_personal_best
      pbest.constraint_violation = self.constraint_violation_personal_best 

      if self >> pbest:
        self.objectives_personal_best = list(self.objectives)
        self.variables_personal_best = list(self.variables)
        self.constraint_violation_personal_best = float(self.constraint_violation)
#        print "### Personal best updated old :" + str(pbest.objectives) + str(pbest.constraint_violation)
#        print "### new :" + str(self.objectives_personal_best) + str(self.constraint_violation_personal_best)
      else:
#        print "#$#$#$ Personal best cannot be updated: " + str(pbest.objectives) + str(pbest.constraint_violation)
        if random.random() > 0.5:
          self.objectives_personal_best = list(self.objectives)
          self.variables_personal_best = list(self.variables)
          self.constraint_violation_personal_best = float(self.constraint_violation)
#          print "#$#$#$ Personal best cannot be updated: " + str(self.objectives_personal_best) + str(self.constraint_violation_personal_best)


def crowded_comparison(s1, s2):
    '''
    Compare the two solutions based on crowded comparison.
    '''
    if s1.rank < s2.rank:
        return 1
        
    elif s1.rank > s2.rank:
        return -1
        
    elif s1.distance > s2.distance:
        return 1
        
    elif s1.distance < s2.distance:
        return -1
        
    else:
        return 0

class Individual(Solution):
    '''
    Solution for the T1 function.
    '''
    def __init__(self):
        '''
        Constructor.
        '''
        Solution.__init__(self)
        for i in range(self.num_variables):
          start_val = 0.5*(self.variable_max[i] + self.variable_min[i]) 
          width = min(abs(self.variable_max[i]-start_val), abs(self.variable_min[i]-start_val))
          rand = random.random()
          self.variables[i] = self.variable_min[i]+ rand*(self.variable_max[i]-self.variable_min[i])
          self.velocities[i] = (2.0*rand-1.0)*(self.variable_max[i]-self.variable_min[i])

        self.maintain_variable_range()
        self.evaluate()
        self.objectives_personal_best = list(self.objectives)
        self.variables_personal_best = list(self.variables)
        self.constraint_violation_personal_best = float(self.constraint_violation)

    def evaluate(self):
	config.eval_func(self) 

    def mutate(self, const):
      if random.random() < const:
#        i = random.randint(0, self.num_variables-1)
        i = int(random.random()*self.num_variables)
        change_range = (self.variable_max[i]-self.variable_min[i])*0.5*const 
        range_min = max(self.variable_min[i], self.variables[i]-change_range)
        range_max = min(self.variable_max[i], self.variables[i]+change_range)
        self.variables[i] = range_min + random.random()*(range_max-range_min)
        #print "---------- mutation ---------------"
        #print "var " + str(i) + ", change_range = " + str(change_range) + ", var_num = " + str(self.num_variables)
        #print "-----------------------------------"    
