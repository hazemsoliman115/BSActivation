# -*- coding: utf-8 -*-
"""
Created on Tue May 19 15:53:28 2015

@author: hazem.soliman
"""
import numpy
import math
from scipy import optimize

class Traffic_Generator(object):
    """ A class to generate traffic """
    def __init__(self, no_periods, max_ratio = 60, min_traf = 1, max_traf = 60, std_devi = 0.1, sleep_time_ratio = 0.25, ramp_up_time_ratio = 0.25, high_time_ratio = 0.25, ramp_down_time_ratio = 0.25):
        if not max_ratio:
            self.max_ratio = max_traf/min_traf
        else:
            self.max_ratio = max_ratio
        
        self.traffic = [0 for i in range(no_periods)]
        for i in range(math.ceil(sleep_time_ratio*no_periods)):
            self.traffic[i] = abs(numpy.random.normal(loc=min_traf, scale=std_devi*max_traf))
        list_offset = math.ceil(sleep_time_ratio*no_periods)
        for i in range(math.ceil(ramp_up_time_ratio*no_periods)):
            self.traffic[min(i+list_offset,no_periods-1)] = abs(numpy.random.normal(loc=min_traf+(i+1)/math.ceil(ramp_up_time_ratio*no_periods)*max_ratio, scale=std_devi*max_traf))
        list_offset = list_offset + math.ceil(ramp_up_time_ratio*no_periods)
        for i in range(math.ceil(high_time_ratio*no_periods)):
            self.traffic[min(i+list_offset,no_periods-1)] = abs(numpy.random.normal(loc=max_traf, scale=std_devi*max_traf))
        list_offset = list_offset + math.ceil(high_time_ratio*no_periods)
        for i in range(math.ceil(ramp_down_time_ratio*no_periods)):
            self.traffic[min(i+list_offset,no_periods-1)] = abs(numpy.random.normal(loc=min_traf+(1-(i+1)/math.ceil(ramp_up_time_ratio*no_periods))*max_ratio, scale=std_devi*max_traf))
            


class Sch_User(object):
    """ A general class for users in the cellular system to be scheduled"""
    
    def __init__(self, pos_x, pos_y, q):
        
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.queue_length = q
        self.BS = None
        
    def attach_BS(self, list_APs):
        """ Find the closest, not necessarily active base station """
        min_dist = 10**6        
        for l in list_APs:
            l_dist = numpy.sqrt((self.pos_x-l.pos_x)**2 - (self.pos_y-l.pos_y)**2)
            if l_dist < min_dist:
                min_dist = l_dist
                self.BS = l
        self.BS.attached_users.append(self)
        
class AccessPoint(object):
    """ A general class for a wireless access point in the scheduling system """
    
    def __init__(self, id, pos_x, pos_y, a = 0):
        self.id = id
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.active_status = a
        self.covered_cells = []
        self.attached_users = []
        self.Sec_attached_users = []
        self.total_Q = 0
        self.clustered_cells = []
        
    def associate_cells(self, intf_g):
        """ Find the set of cells within transmission range, main cell plus interfering cells """
        self.covered_cells.append(self.id)
        for l in range(len(intf_g_obj.intf_g[self.id])):
            if intf_g_obj.intf_g[self.id][l]==1:
               self.covered_cells.append(l) 
               
    def attached_Secondary_users(self, list_APs):
        """ A function that goes over interfering cells, checks whether they are active, if not it attachs its users to the current AP  """
        for l in self.covered_cells:
            if list_APs[l].active_status == 0:
                self.attached_Secondary_users = self.attached_Secondary_users + list_APs[l].attached_users
              
              
    def Find_Sum_Qs(self):
        """ A function to find the total demand as queue lengths attached to the BSs """
        self.total_Q = 0        
        for u in self.attached_users:
            self.total_Q = self.total_Q + u.queue_length
            
        
class Interference_Graph(object):
    """ A general graph to store the interference graph """
    
    def __init__(self, list_APs, intf_threshold):
        intf_g = [[0 for i in range(len(list_APs))] for i in range(len(list_APs))]
        intf_g_power = [[0 for i in range(len(list_APs))] for i in range(len(list_APs))]
        ap_index = 0
        for ap in list_APs:
            i_ap_index = 0
            for i_ap in list_APs:
                if numpy.sqrt((ap.pos_x-i_ap.pos_x)**2 + (ap.pos_y-i_ap.pos_y)**2) <= intf_threshold and ap is not i_ap:
                    intf_g[ap_index][i_ap_index] = 1
                    intf_g[i_ap_index][ap_index] = 1
                    intf_g_power[ap_index][i_ap_index] = numpy.sqrt((ap.pos_x-i_ap.pos_x)**2 + (ap.pos_y-i_ap.pos_y)**2)/intf_threshold
                    intf_g_power[i_ap_index][ap_index] = numpy.sqrt((ap.pos_x-i_ap.pos_x)**2 + (ap.pos_y-i_ap.pos_y)**2)/intf_threshold
                i_ap_index+=1
            ap_index+=1
            
        self.intf_g = intf_g
        self.intf_g_power = intf_g_power
        self.Max_f = 0
        
    def find_max_freq(self):
        """ A function to find the maximum frequency of occurrence in the graph """
        self.Max_f = 0        
        for l in self.intf_g:
            temp_max = sum(l)+1
            if self.Max_f <= temp_max:
                self.Max_f = temp_max
                
    def find_intf_g_matrix_opt(self):
        intf_g_matrix_opt = [[0 for i in range(len(self.intf_g))] for i in range(len(self.intf_g))]
        for k in range(len(intf_g_matrix_opt)):
            intf_g_matrix_opt[k][k] = -1
            for m in range(len(intf_g_matrix_opt)):
                if self.intf_g[k][m] == 1:
                    intf_g_matrix_opt[k][m] = -1
        return(intf_g_matrix_opt)
                
class AP_Cluster(object):
    """ A class for modelling clusters of APs """
    def __init__(self):
        self.clustered_AP_list = []
        
        
if __name__ == "__main__":
    No_APs = 7 # number of access points
    AP_radius = 50 # cell radius
    No_users = 100 # number of users
    No_RBs = 100 # number of PRBs
    No_bits_CSI = 2 # number of bits used to represent CSI
    No_bits_RB = 36 #No bits transmitted in each RB, assuming SINR = 10 dB, 12 subcarriers per RB. and using 16 QAM and code rate = 4/5
    Avg_Q_Size = 1000 # Average number of bits in a queue
    intf_threshold = 50
    intf_g_matrix = [[0, 0, 1, 1, 0, 0, 0],[ 0, 0, 0, 1, 1, 0, 0],[ 1, 0, 0, 0, 0, 1, 0],[ 1, 1, 0, 0, 0, 1, 1],[0, 1, 0, 0, 0, 0, 1],[ 0, 0, 1, 1, 0, 0, 0],[ 0, 0, 0, 1, 1, 0, 0]]
    intf_g_matrix_opt = [[-1, 0, -1, -1, 0, 0, 0],[ 0, -1, 0, -1, -1, 0, 0],[ -1, 0, -1, 0, 0, -1, 0],[ -1, -1, 0, -1, 0, -1, -1],[0, -1, 0, 0, -1, 0, -1],[ 0, 0, -1, -1, 0, -1, 0],[ 0, 0, 0, -1, -1, 0, -1]]
    list_users = []
    list_APs = []
    
    for i in range(No_users): # create users objects
        u = Sch_User(numpy.random.uniform(low=0.0, high=No_APs*AP_radius), numpy.random.uniform(low=-AP_radius, high=AP_radius, size=None), numpy.random.normal(loc=Avg_Q_Size, scale=1.0))
        list_users.append(u)
        
    for i in range(No_APs): # create access points objects
        ap = AccessPoint(i, numpy.random.uniform(low=0.0, high=No_APs*AP_radius), numpy.random.uniform(low=-AP_radius, high=AP_radius, size=None))
        list_APs.append(ap)
        
    for u in list_users: # attach users to default base stations
        u.attach_BS(list_APs)
        
    for ap in list_APs: # Find total queue size for all users
        ap.Find_Sum_Qs()
        
    # Create interference graph object
    intf_g_obj = Interference_Graph(list_APs, intf_threshold)
    intf_g_matrix_opt = intf_g_obj.find_intf_g_matrix_opt()
    intf_g_obj.find_max_freq()
    
    # Find interfering cells for each cell
    for ap in list_APs:
        ap.associate_cells(intf_g_obj)
        
    
    # LP parameters
    c = [1, 1, 1, 1, 1, 1, 1]
    b = [-1, -1, -1, -1, -1, -1, -1]
    bounds = ((0,1), (0,1), (0,1), (0,1), (0,1), (0,1), (0,1))
    res = optimize.linprog(c=c, A_ub=intf_g_matrix_opt, b_ub=b, A_eq=None, b_eq=None, bounds=bounds, method='simplex', callback=None, options=None)
        
    print(res)
    # Process linear programming output
    x = [1 if res.x[i] >= 1/intf_g_obj.Max_f else 0 for i in range(len(res.x))]
    for ap in range(len(list_APs)):
        list_APs[ap].active_status = x[ap]
    
    # Sort APs according to Users buffers
    Sorted_APs_Q = sorted(list_APs, key = lambda AccessPoint: AccessPoint.total_Q, reverse=True)
    
    AP_Cluster_Object = AP_Cluster()
    
    for ap in Sorted_APs_Q:
        if ap.active_status == 0:
            exist_flag = 0
            for i_ap_ind in ap.covered_cells:
                for l_cluster in range(len(AP_Cluster_Object.clustered_AP_list)):
                    if list_APs[i_ap_ind] in AP_Cluster_Object.clustered_AP_list[l_cluster]:
                        exist_flag = 1
                        exist_loc = l_cluster
            if exist_flag == 1:
                for i_ap_ind in ap.covered_cells:
                    if list_APs[i_ap_ind] not in AP_Cluster_Object.clustered_AP_list[exist_loc]:
                        AP_Cluster_Object.clustered_AP_list[exist_loc].append(list_APs[i_ap_ind])
            else:
                AP_Cluster_Object.clustered_AP_list.append([list_APs[i_ap_ind] for i_ap_ind in ap.covered_cells])
                
    T = Traffic_Generator(24)
                
    
    
        
    
    