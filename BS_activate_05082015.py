# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 11:51:11 2015

@author: hazem.soliman
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 11:11:51 2015

@author: hazem.soliman
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 16:02:14 2015

@author: hazem.soliman
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 12:38:27 2015

@author: hazem.soliman
"""

# -*- coding: utf-8 -*-
"""
Created on Thu May 21 13:11:53 2015

@author: hazem.soliman
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May 19 15:53:28 2015

@author: hazem.soliman
"""
import numpy
import math
import copy
import matplotlib.pyplot as plt
from scipy import optimize
import operator


class Traffic_Generator(object):
    """ A class to generate traffic """
    def __init__(self, no_periods, max_ratio = 60, min_traf = 1, max_traf = 60, std_devi = 0.1, sleep_time_ratio = 0.25, ramp_up_time_ratio = 0.25, high_time_ratio = 0.25, ramp_down_time_ratio = 0.25):
        if not max_ratio:
            self.max_ratio = max_traf/min_traf
        else:
            self.max_ratio = max_ratio
        
        #Divide traffic accroding to the number of subtimes, in our case it is always 4
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
            
            
class Traffic_Generator_Ramp(object):
    """ A class to generate traffic """
    def __init__(self, no_periods, max_ratio = 15, std_devi = 0.25):

        self.max_ratio = max_ratio
        
        #Divide traffic accroding to the number of subtimes, in our case it is always 4
        self.traffic = [0 for i in range(no_periods)]
        for i in range(no_periods):
            #self.traffic[i] = abs(numpy.random.normal(loc=i*max_ratio/no_periods, scale=std_devi*max_ratio))
            self.traffic[i] = (i+1)*max_ratio/no_periods
        


class Sch_User(object):
    """ A general class for users in the cellular system to be scheduled"""
    
    def __init__(self, i, pos_x, pos_y, q):
        
        self.id = i        
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.queue_length = q
        self.BS = None
        self.best_Tx = None
        self.best_Tx_Cluster = None
        self.best_Tx_dist = None
        self.worst_InTf = None
        self.worst_InTf_cluster = None
        self.worst_Intf_dist = None
        
    def attach_BS(self, list_APs):
        """ Find the closest, not necessarily active base station """
        min_dist = 10**20   
        # loop over the BSs and find the closest one
        for l in list_APs:
            l_dist = numpy.sqrt((self.pos_x-l.pos_x)**2 + (self.pos_y-l.pos_y)**2)
            if l_dist < min_dist:
                min_dist = l_dist
                self.BS = l
        self.BS.attached_users.append(self)
        
    def Find_Tx(self, list_APs, AP_Cluster_Object, AP_radius):
        """ A function to find the best transmitting active AP for this user """
        min_dist = 10**6   
        # loop over the BSs and find the best ACTIVE one
        for l in list_APs:
            l_dist = numpy.sqrt((self.pos_x-l.pos_x)**2 + (self.pos_y-l.pos_y)**2)
            if l_dist < min_dist and l.active_status == 1:
                min_dist = l_dist
                self.best_Tx = l
                self.best_Tx_dist = numpy.sqrt((self.pos_x-l.pos_x)**2 + (self.pos_y-l.pos_y)**2)/AP_radius
        # if best transmitted is different from original BS, attach the user as a secondary user
        if self.BS != self.best_Tx:
            self.best_Tx.Sec_attached_users.append(self)
        elif self not in self.best_Tx.attached_users:
            self.best_Tx.attached_users.append(self)
        # Find the cluster it belongs
        for l_ap in AP_Cluster_Object.clustered_AP_list:
            if self.best_Tx in l_ap:
                self.best_Tx_Cluster = l_ap
        
    def Find_InTf(self, list_APs, AP_Cluster_Object, AP_radius):
        """ A function to find the worst interfering access point for this user """
        min_dist = 10**6  
        # loop over the BSs and find the worst interfering active one
        for l in list_APs:
            l_dist = numpy.sqrt((self.pos_x-l.pos_x)**2 + (self.pos_y-l.pos_y)**2)
            if l_dist < min_dist and l.active_status == 1 and not l == self.best_Tx and not l in self.best_Tx_Cluster:
                min_dist = l_dist
                self.worst_InTf = l
                self.worst_Intf_dist = numpy.sqrt((self.pos_x-l.pos_x)**2 + (self.pos_y-l.pos_y)**2)/AP_radius
        if self.worst_Intf_dist is None:
            self.worst_Intf_dist = 1
        
        # Find its cluster
        for l_ap in AP_Cluster_Object.clustered_AP_list:
            if self.worst_InTf in l_ap:
                self.worst_InTf_cluster = l_ap
        
class AccessPoint(object):
    """ A general class for a wireless access point in the scheduling system """
    
    def __init__(self, i, pos_x, pos_y, a = 0):
        self.id = i
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.active_status = a
        self.covered_cells = []
        self.Intf_cells = []
        self.attached_users = []
        self.Sec_attached_users = []
        self.total_Q = None
        self.clustered_cells = []
        
    def associate_cells(self, intf_g):
        """ Find the set of cells within transmission range, main cell plus interfering cells """
        self.covered_cells.append(self.id)
        for l in range(len(intf_g_obj.intf_g[self.id])):
            if intf_g_obj.intf_g[self.id][l]==1:
               self.covered_cells.append(l)
               self.Intf_cells.append(l)
               
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
            
    def Activate_ActorCritique(self, activ_threshold):
        """ A function that determines the activation status of the base station based on the number of attached users and the threshold """
        if len(self.attached_users) >= activ_threshold:
            self.active_status = 1
        else:
            self.active_status = 0
            
            
        
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
        """ A function to find the maximum frequency of occurrence in the graph, for the rounding """
        self.Max_f = 0        
        for l in self.intf_g:
            temp_max = sum(l)+1
            if self.Max_f <= temp_max:
                self.Max_f = temp_max
                
    def find_intf_g_matrix_opt(self):
        """ Prepares the matrix to be used in the optimization problrm """
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
        
    def Continue_Activ_Cluster(self, list_APs, list_Active_APs, list_users, p_W, s_W, u_W, BW_W, No_APs, max_ratio):
        """ A function that continues the optimization after the set cover problem """
        #print("Begin")
        for ap in  list_Active_APs:
            self.clustered_AP_list.append([ap])
        list_cluster_AP = []
        for k in self.clustered_AP_list:
            if len(k) > 1:
                for x in k:
                    list_cluster_AP.append(x)
        while(1 and len(list_Active_APs)<len(list_APs) and len(list_cluster_AP)<len(list_Active_APs)):
            #print("Inside")
            #print(list_Active_APs)
#            for u in list_users:
#                u.BS = None
#                u.best_Tx = None
#                u.best_Tx_Cluster = None
#                u.best_Tx_dist = None
#                u.worst_InTf = None
#                u.worst_InTf_cluster = None
#                u.worst_Intf_dist = None


            temp_list_users_P = [copy.copy(element) for element in list_users]
            temp_list_users_C = [copy.copy(element) for element in list_users]
            #temp_list_users_P = copy.deepcopy(list_users)
            #temp_list_users_C = copy.deepcopy(list_users)


#            for u in list_users:
#                u.attach_BS(list_APs)
#            for u in temp_list_users_P: # attach users to default base stations
#                u.attach_BS(list_APs)
#            for u in temp_list_users_C: # attach users to default base stations
#                u.attach_BS(list_APs)
            total_utility_power = 0
            total_utility_cluster = 0
            # find worst ap
            for ap in list_APs: # clear attached users
                ap.attached_users = []
                ap.Sec_attached_users = []
            user_power = 0
            user_intf = 0
            user_SIR = 0
            for u in list_users:
                u.Find_Tx(list_APs, self, AP_radius)
                user_power = user_power + u.best_Tx_dist
                u.Find_InTf(list_APs, self, AP_radius)
                user_intf = user_intf + u.worst_Intf_dist
                #print("P=",u.best_Tx_dist)
                #print("I=",u.worst_Intf_dist)
                #print("SIR=",(u.worst_Intf_dist/u.best_Tx_dist)**3)
                user_SIR += (u.worst_Intf_dist/u.best_Tx_dist)**1
            #print("P=",user_power)
            #print("I=",user_intf)
            worst_ap_P = None
            worst_ap_weight_P = -1*float("inf")
            for ap in list_APs:
                temp_weight = 0
                for u in ap.attached_users:
                    temp_weight += u.worst_Intf_dist/u.best_Tx_dist
                if temp_weight > worst_ap_weight_P and ap.active_status == 0:
                    worst_ap_weight_P = temp_weight
                    worst_ap_P = ap
            # change in power
            total_utility_power = total_utility_power - p_W/No_APs
            # change in number of users
            total_utility_power = total_utility_power - (-u_W*(len(list_users)/len(list_Active_APs))/max_ratio + u_W*(len(list_users)/(len(list_Active_APs)+1))/max_ratio)
            # change in received signal power and received interference power
            #QoS = (user_intf/len(list_users))**3/(user_power/len(list_users))**3
            QoS = user_SIR/len(list_users)
            worst_ap_P.active_status = 1
            self.clustered_AP_list.append([worst_ap_P])
            for ap in list_APs: # clear attached users
                ap.attached_users = []
                ap.Sec_attached_users = []
            temp_user_power = 0
            temp_user_intf = 0
            temp_user_SIR = 0
            for u in temp_list_users_P:
                u.Find_Tx(list_APs, self, AP_radius)
                temp_user_power = temp_user_power + u.best_Tx_dist
                u.Find_InTf(list_APs, self, AP_radius)
                temp_user_intf = temp_user_intf + u.worst_Intf_dist
                temp_user_SIR += (u.worst_Intf_dist/u.best_Tx_dist)**1
            #print("TP=",temp_user_power)
            #print("TI=",temp_user_intf)
            #temp_QoS_P = (temp_user_intf/len(temp_list_users_P))**3/(temp_user_power/len(list_users))**3
            temp_QoS_P = temp_user_SIR/len(temp_list_users_P)
            
            total_utility_power = total_utility_power - s_W*QoS/10 + s_W*temp_QoS_P/10
            
            # Prepare for clustering by removing what you did
            self.clustered_AP_list.remove([worst_ap_P])
            worst_ap_P.active_status = 0
            
            
            # Check Clustering
            for ap in list_APs: # clear attached users
                ap.attached_users = []
                ap.Sec_attached_users = []
            for u in list_users:
                u.Find_Tx(list_APs, self, AP_radius)
                user_power = user_power + u.best_Tx_dist
                u.Find_InTf(list_APs, self, AP_radius)
                user_intf = user_intf + u.worst_Intf_dist
            list_InTfs = {ap : 0 for ap in list_Active_APs}
            for u in list_users:
                if u.worst_InTf:
                    list_InTfs[u.worst_InTf] += 1
            sorted_list_InTfs = sorted(list_InTfs, key = list_InTfs.get, reverse = True)
            worst_ap_C = None
            for ap in sorted_list_InTfs:
                if [ap] in self.clustered_AP_list:
                    worst_ap_C = ap
                    #print("Yes")
                    #print([worst_ap_C])
                    #print(self.clustered_AP_list)
                    break
            #worst_ap_C = sorted_list_InTfs[0]
            
            #print(worst_ap_C in list_cluster_AP)
            # change in BW
            nearest_dist = 1*float("inf")
            nearest_ap = None
            #for n_ap_ind in worst_ap_C.Intf_cells:
            for n_ap_ind in range(len(list_APs)):
                temp_dist = (ap.pos_x-list_APs[n_ap_ind].pos_x)**2 + (ap.pos_y-list_APs[n_ap_ind].pos_y)**2
                if temp_dist <= nearest_dist and list_APs[n_ap_ind].active_status == 1 and list_APs[n_ap_ind] is not worst_ap_C:
                    nearest_dist = temp_dist
                    nearest_ap = list_APs[n_ap_ind]
                    
#            print(nearest_ap)
#            print(worst_ap_C)
#            print(self.clustered_AP_list)
#            print([len(h) for h in self.clustered_AP_list])
            if len(list_Active_APs) == 1:
                BW = 0
                temp_BW = 0
            for l in self.clustered_AP_list:
                if nearest_ap in l:
                    BW = len(l)
                    #print(BW)
                    l.append(worst_ap_C)
                    #print([worst_ap_C])
                    #print(self.clustered_AP_list)
                    self.clustered_AP_list.remove([worst_ap_C])
                    temp_BW = len(l)
                    break
                        
                    
            # change in received interference
            for ap in list_APs: # clear attached users
                ap.attached_users = []
                ap.Sec_attached_users = []
            temp_user_power = 0
            temp_user_intf = 0
            temp_user_SIR = 0
            for u in temp_list_users_C:
                u.Find_Tx(list_APs, self, AP_radius)
                temp_user_power = temp_user_power + u.best_Tx_dist
                u.Find_InTf(list_APs, self, AP_radius)
                temp_user_intf = temp_user_intf + u.worst_Intf_dist
                temp_user_SIR += (u.worst_Intf_dist/u.best_Tx_dist)**1
            #print("CP=",temp_user_power)
            #print("CI=",temp_user_intf)
            #temp_QoS_C = (temp_user_intf/len(temp_list_users_C))**3/(temp_user_power/len(temp_list_users_C))**3
            temp_QoS_C = temp_user_SIR/len(temp_list_users_C)
            total_utility_cluster = total_utility_cluster - s_W*QoS/10 + s_W*temp_QoS_C/10
#            print([len(h) for h in self.clustered_AP_list])
#            print("QoS = ", QoS)
#            print("Temp QoS = ", temp_QoS_C)
            total_utility_cluster = total_utility_cluster -(-BW_W*BW/No_APs + BW_W*temp_BW/No_APs)
            self.clustered_AP_list.append([worst_ap_C])
            l.remove(worst_ap_C)
            
            #Decision
            
#            print("Change in Power = ", p_W/No_APs)
#            print("Change in number of user", u_W*(len(list_users)/len(list_Active_APs))/max_ratio)
#            print("Change in number of user", u_W*(u_W*(len(list_users)/(len(list_Active_APs)+1))/max_ratio))
#            print("Change in QoS",  s_W*QoS/10)
#            print("Change in QoS",s_W*temp_QoS_P/10)
#            print("Power utility", total_utility_power)
#            print("Change in BW", (BW_W*BW/No_APs))
#            print("Change in BW", (BW_W*temp_BW/No_APs))
#            print("Change in Cluster QoS",  s_W*QoS/10 )
#            print("Change in Cluster QoS",  s_W*temp_QoS_C/10)
#            print("Cluster Utility",total_utility_cluster)
            #print(len(l))
            #print(self.clustered_AP_list)
            #print(worst_ap_C)
            if total_utility_power<=-0.04 and total_utility_cluster <= -0.04:
                break
            elif total_utility_power > total_utility_cluster:
                worst_ap_P.active_status = 1
                self.clustered_AP_list.append([worst_ap_P])
                #print(len(list_users))
                list_users = temp_list_users_P
                #print(len(list_users))
                temp_list_users_P = None
                temp_list_users_C = None
                list_Active_APs = []
                for ap in list_APs:
                    if ap.active_status == 1:
                        list_Active_APs.append(ap)
                #print("Decision = P")
            elif total_utility_cluster > total_utility_power:
                l.append(worst_ap_C)
                self.clustered_AP_list.remove([worst_ap_C])
                #print(worst_ap_C)
                #print(l)
                #print(len(list_users))
                list_users = temp_list_users_C
                #print(len(list_users))
                temp_list_users_P = None
                temp_list_users_C = None
                list_Active_APs = []
                for ap in list_APs:
                    if ap.active_status == 1:
                        list_Active_APs.append(ap)
                list_cluster_AP = []
                for k in self.clustered_AP_list:
                    if len(k) > 1:
                        for x in k:
                            list_cluster_AP.append(x)
                #print("Decision = C")
                v=1
                
                

        
        
        
                
        
            
        
    def Create_Cluster_LinProg(self, Sorted_APs_Q):
        for ap in Sorted_APs_Q:
            if ap.active_status == 0:
                exist_flag = 0
                for i_ap_ind in ap.Intf_cells:
                    for l_cluster in range(len(self.clustered_AP_list)):
                        if list_APs[i_ap_ind] in self.clustered_AP_list[l_cluster]:
                            exist_flag = 1
                            exist_loc = l_cluster
                if exist_flag == 1:
                    for i_ap_ind in ap.Intf_cells:
                        if list_APs[i_ap_ind] not in self.clustered_AP_list[exist_loc]:
                            self.clustered_AP_list[exist_loc].append(list_APs[i_ap_ind])
                else:
                    self.clustered_AP_list.append([list_APs[i_ap_ind] for i_ap_ind in ap.Intf_cells])
        for ap in Sorted_APs_Q:
            exist_flag = 0
            if ap.active_status == 1:
                for l_cluster in range(len(self.clustered_AP_list)):
                        if ap in self.clustered_AP_list[l_cluster]:
                            exist_flag = 1
                if exist_flag == 0:
                    self.clustered_AP_list.append([ap])
                    
    def Find_Cluster(self, ap):
        """ A function to find which cluster the AP belongs to, if any, and return the index of the cluster """
        cluster_index = 100
        exist_flag = 0
        for l_cluster in range(len(self.clustered_AP_list)):
            if ap in self.clustered_AP_list[l_cluster]:
                exist_flag = 1
                cluster_index = l_cluster
                
        return(exist_flag, cluster_index)
        
    def Create_Cluster_DynamicProg(self, intf_g_obj, list_APs, InTf_policy):
        """A function to create the cluster based on the dynamic programming threshold """
        local_intf_g_obj = copy.deepcopy(intf_g_obj)
        while(1):
            max_intf_link = 0
            for x in range(len(list_APs)):
                for y in range(len(list_APs)):
                    if local_intf_g_obj.intf_g_power[x][y] > max_intf_link and local_intf_g_obj.intf_g_power[x][y] > InTf_policy:
                        max_intf_link = local_intf_g_obj.intf_g_power[x][y]
                        ap1_ind = x
                        ap2_ind = y
                        
            if max_intf_link == 0:
                break
            
            ap1_exist,ap1_cluster_index = self.Find_Cluster(list_APs[ap1_ind])
            ap2_exist,ap2_cluster_index = self.Find_Cluster(list_APs[ap2_ind])
            
            
            if list_APs[ap1_ind].active_status == 1 and  list_APs[ap2_ind].active_status == 1:
                if ap1_exist == 0 and ap2_exist == 0:
                    self.clustered_AP_list.append([list_APs[ap1_ind], list_APs[ap2_ind]])
                elif ap1_exist == 0 and ap2_exist == 1:
                    self.clustered_AP_list[ap2_cluster_index].append(list_APs[ap1_ind])
                elif ap1_exist == 1 and ap2_exist == 0:
                    self.clustered_AP_list[ap1_cluster_index].append(list_APs[ap2_ind])

            local_intf_g_obj.intf_g_power[ap1_ind][ap2_ind] = 0
            local_intf_g_obj.intf_g_power[ap2_ind][ap1_ind] = 0
            

                    
            
        
        
        
        
if __name__ == "__main__":
    No_APs = 20 # number of access points
    AP_radius = 50 # cell radius
    AP_concentration = 0.01
    #No_users = 100 # number of users
    No_RBs = 100 # number of PRBs
    No_bits_CSI = 2 # number of bits used to represent CSI
    No_bits_RB = 36 #No bits transmitted in each RB, assuming SINR = 10 dB, 12 subcarriers per RB. and using 16 QAM and code rate = 4/5
    Avg_Q_Size = 1000 # Average number of bits in a queue
    #intf_threshold = 50
    no_periods = 10 
    max_ratio = 10
    min_traf = 1
    max_traf = 16
    std_devi = 0.1
    sleep_time_ratio = 0.25
    ramp_up_time_ratio = 0.25
    high_time_ratio = 0.25
    ramp_down_time_ratio = 0.25
    Intf_policy_limit = 5
    no_sample_paths = 10000
    #intf_g_matrix = [[0, 0, 1, 1, 0, 0, 0],[ 0, 0, 0, 1, 1, 0, 0],[ 1, 0, 0, 0, 0, 1, 0],[ 1, 1, 0, 0, 0, 1, 1],[0, 1, 0, 0, 0, 0, 1],[ 0, 0, 1, 1, 0, 0, 0],[ 0, 0, 0, 1, 1, 0, 0]]
    #intf_g_matrix_opt = [[-1, 0, -1, -1, 0, 0, 0],[ 0, -1, 0, -1, -1, 0, 0],[ -1, 0, -1, 0, 0, -1, 0],[ -1, -1, 0, -1, 0, -1, -1],[0, -1, 0, 0, -1, 0, -1],[ 0, 0, -1, -1, 0, -1, 0],[ 0, 0, 0, -1, -1, 0, -1]]
    list_users = []
    list_APs = []
    
    
        
#    for i in range(No_APs): # create access points objects
#        ap = AccessPoint(i, numpy.random.uniform(low=-0.5*numpy.sqrt(AP_concentration*No_APs*AP_radius), high=0.5*numpy.sqrt(AP_concentration*No_APs*AP_radius)), numpy.random.uniform(-0.5*numpy.sqrt(AP_concentration*No_APs*AP_radius), high=0.5*numpy.sqrt(AP_concentration*No_APs*AP_radius)))       
#        list_APs.append(ap)
#        
#    
#        
#    for ap in list_APs: # Find total queue size for all users
#        ap.Find_Sum_Qs()
#        
#    # Create interference graph object
#    intf_g_obj = Interference_Graph(list_APs, intf_threshold)
#    intf_g_matrix_opt = intf_g_obj.find_intf_g_matrix_opt()
#    intf_g_obj.find_max_freq()
#    
#    # Find interfering cells for each cell
#    for ap in list_APs:
#        ap.associate_cells(intf_g_obj)
#    
#    
#    AP_Cluster_Object = AP_Cluster()            
#    
#                
#    T = Traffic_Generator(no_periods, max_ratio = 16, min_traf = 1, max_traf = 16, std_devi = 0.25, sleep_time_ratio = 0.25, ramp_up_time_ratio = 0.25, high_time_ratio = 0.25, ramp_down_time_ratio = 0.25)
    
#    intf_threshold_list = [10,15,20,25,30,35,40,45,50]
    intf_threshold_list = [10,20,30,40,50]
    Tx_Policy_Output_Power = [[0 for i in range(len(intf_threshold_list))] for state_hour in range(no_periods)]
    Tx_Policy_Output_Power_per_user = [[0 for i in range(len(intf_threshold_list))] for state_hour in range(no_periods)]
    Tx_Policy_Output_QoS = [[0 for i in range(len(intf_threshold_list))] for state_hour in range(no_periods)]
    Tx_Policy_Output_QoS_No_Cluster = [[0 for i in range(len(intf_threshold_list))] for state_hour in range(no_periods)]
    Tx_Policy_Output_AvgNoUser = [[0 for i in range(len(intf_threshold_list))] for state_hour in range(no_periods)]
    
#    p_W_list = [1]
#    s_W_list = [1]
#    u_W_list = [1]
#    BW_W_list = [1]
#    Tx_Power_Power_Weight = [[0 for ij in range(len(p_W_list))] for state_hour in range(no_periods)]
#    Tx_Power_Signal_Weight = [[0 for ij in range(len(s_W_list))] for state_hour in range(no_periods)]
#    Tx_Power_User_Weight = [[0 for ij in range(len(u_W_list))] for state_hour in range(no_periods)]
#    Tx_Power_BW_Weight = [[0 for ij in range(len(BW_W_list))] for state_hour in range(no_periods)]
#    
#    Tx_QoS_Power_Weight = [[0 for ij in range(len(p_W_list))] for state_hour in range(no_periods)]
#    Tx_QoS_Signal_Weight = [[0 for ij in range(len(s_W_list))] for state_hour in range(no_periods)]
#    Tx_QoS_User_Weight = [[0 for ij in range(len(u_W_list))] for state_hour in range(no_periods)]
#    Tx_QoS_BW_Weight = [[0 for ij in range(len(BW_W_list))] for state_hour in range(no_periods)]
#    
#    Tx_User_Power_Weight = [[0 for ij in range(len(p_W_list))] for state_hour in range(no_periods)]
#    Tx_User_Signal_Weight = [[0 for ij in range(len(s_W_list))] for state_hour in range(no_periods)]
#    Tx_User_User_Weight = [[0 for ij in range(len(u_W_list))] for state_hour in range(no_periods)]
#    Tx_User_BW_Weight = [[0 for ij in range(len(BW_W_list))] for state_hour in range(no_periods)]
    
#    for p_W_ind in range(len(p_W_list)):    
#        for s_W_ind in range(len(s_W_list)):
#            for u_W_ind in range(len(u_W_list)):
#                for BW_W_ind in range(len(BW_W_list)):
    for intf_threshold_index in range(len(intf_threshold_list)):
        for path_w in range(no_sample_paths):
            print(intf_threshold_index, path_w)
            list_users = []
            list_APs = []
            for i in range(No_APs): # create access points objects
                ap = AccessPoint(i, numpy.random.uniform(low=-0.5*numpy.sqrt(AP_concentration*No_APs*numpy.pi*AP_radius**2), high=0.5*numpy.sqrt(AP_concentration*No_APs*numpy.pi*AP_radius**2)), numpy.random.uniform(-0.5*numpy.sqrt(AP_concentration*No_APs*numpy.pi*AP_radius**2), high=0.5*numpy.sqrt(AP_concentration*No_APs*numpy.pi*AP_radius**2)))       
                list_APs.append(ap)
            
        
            
            for ap in list_APs: # Find total queue size for all users
                ap.Find_Sum_Qs()
                
            # Create interference graph object
            intf_g_obj = None
            intf_g_obj = Interference_Graph(list_APs, intf_threshold_list[intf_threshold_index])
            intf_g_matrix_opt = intf_g_obj.find_intf_g_matrix_opt()
           
            intf_g_obj.find_max_freq()
            
            # Find interfering cells for each cell
            for ap in list_APs:
                ap.associate_cells(intf_g_obj)
            
            
            AP_Cluster_Object = AP_Cluster()
            T = Traffic_Generator_Ramp(no_periods, max_ratio, std_devi)
            No_users = []        
            for i in range(no_periods): # Loop to get number of users from generated data
                No_users.append(math.ceil(T.traffic[i]*No_APs))
            for state_hour in range(no_periods):
                #print("Hour = ", state_hour)
                list_users = []
                for i in range(No_users[state_hour]): # create users objects
                    u = Sch_User(i, numpy.random.uniform(low=-0.5*numpy.sqrt(AP_concentration*No_APs*numpy.pi*AP_radius**2), high=0.5*numpy.sqrt(AP_concentration*No_APs*numpy.pi*AP_radius**2)), numpy.random.uniform(low=-0.5*numpy.sqrt(AP_concentration*No_APs*numpy.pi*AP_radius**2), high=0.5*numpy.sqrt(AP_concentration*No_APs*numpy.pi*AP_radius**2)), numpy.random.normal(loc=Avg_Q_Size, scale=1.0))
                    list_users.append(u)
                for ap in list_APs: # clear attached users
                    ap.attached_users = []
                for u in list_users: # attach users to default base stations
                    u.attach_BS(list_APs)
    
                list_Active_APs = []
                for ap in list_APs: # clear attached users
                    ap.attached_users = []
                    ap.Sec_attached_users = []
                for u in list_users: # attach users to default base stations
                    u.attach_BS(list_APs)
                
                # LP parameters
                c = [1 for i in range(len(list_APs))]
                b = [-1 for i in range(len(list_APs))]
                bounds = tuple((0,1) for i in range(len(list_APs)))
                res = optimize.linprog(c=c, A_ub=intf_g_matrix_opt, b_ub=b, A_eq=None, b_eq=None, bounds=bounds, method='simplex', callback=None, options=None)
                    
                #print(path_w)
                # Process linear programming output
                x = [1 if res.x[i] >= 1/intf_g_obj.Max_f else 0 for i in range(len(res.x))]
                
                
                for ap in range(len(list_APs)):
                    list_APs[ap].active_status = x[ap]
            
            
                for ap in list_APs:
                    if ap.active_status == 1:
                        list_Active_APs.append(ap)
                 # Sort APs according to Users buffers
                p_W = 1
                s_W = 1
                u_W = 1
                BW_W = 1
                AP_Cluster_Object = AP_Cluster()
                #try:
                AP_Cluster_Object.Continue_Activ_Cluster( list_APs, list_Active_APs, list_users, p_W, s_W, u_W, BW_W, No_APs, max_ratio)
                #except:
                    #print("Error")
                    #continue
                Sorted_APs_Q = sorted(list_APs, key = lambda AccessPoint: AccessPoint.total_Q, reverse=True)
                
                #AP_Cluster_Object.Create_Cluster_LinProg(Sorted_APs_Q)
                list_Active_APs = []
                for ap in list_APs:
                    if ap.active_status == 1:
                        list_Active_APs.append(ap)
                avg_cluster_len = 0
                for l in AP_Cluster_Object.clustered_AP_list:
                    avg_cluster_len = avg_cluster_len + len(l)
    #                    print("Intf_Policy",  Intf_policy_list[InTf_policy])
                for ap in list_APs: # clear attached users
                    ap.attached_users = []
                    ap.Sec_attached_users = []
                user_power = 0
                user_intf = 0
                user_SIR = 0
                for u in list_users:
                    u.Find_Tx(list_APs, AP_Cluster_Object, AP_radius)
                    #print(u.best_Tx_dist)
                    user_power = user_power + u.best_Tx_dist
                    u.Find_InTf(list_APs, AP_Cluster_Object, AP_radius)
                    #print(u.worst_Intf_dist)
                    user_intf = user_intf + u.worst_Intf_dist
                    user_SIR += (u.worst_Intf_dist/u.best_Tx_dist)**1
                power = len(list_Active_APs)
                #print(power)
                #QoS = (user_intf/len(list_users))**3/(user_power/len(list_users))**3
                QoS = user_SIR/len(list_users)
                # No Clustering Case
                AP_Cluster_Object.clustered_AP_list = [[u] for u in list_Active_APs]
                user_power = 0
                user_intf = 0
                user_SIR = 0
                for u in list_users:
                    u.Find_Tx(list_APs, AP_Cluster_Object, AP_radius)
                    #print(u.best_Tx_dist)
                    user_power = user_power + u.best_Tx_dist
                    u.Find_InTf(list_APs, AP_Cluster_Object, AP_radius)
                    #print(u.worst_Intf_dist)
                    user_intf = user_intf + u.worst_Intf_dist
                    user_SIR += (u.worst_Intf_dist/u.best_Tx_dist)**1
                #No_Cluster_QoS = (user_intf/len(list_users))**3/(user_power/len(list_users))**3
                No_Cluster_QoS = user_SIR/len(list_users) 
                avgnouser = 0
                for ap in list_Active_APs:
                    #print(len(ap.attached_users))
                    avgnouser = avgnouser + len(ap.attached_users) + len(ap.Sec_attached_users)
                Tx_Policy_Output_Power[state_hour][intf_threshold_index] += power
                Tx_Policy_Output_Power_per_user[state_hour][intf_threshold_index] += power/len(list_users)
#                Tx_Power_Power_Weight[state_hour][p_W_ind] += power
#                Tx_Power_Signal_Weight[state_hour][s_W_ind] += power
#                Tx_Power_User_Weight[state_hour][u_W_ind] += power
#                Tx_Power_BW_Weight[state_hour][BW_W_ind] += power
                Tx_Policy_Output_QoS[state_hour][intf_threshold_index] += QoS
                Tx_Policy_Output_QoS_No_Cluster[state_hour][intf_threshold_index] += No_Cluster_QoS
#                Tx_QoS_Power_Weight[state_hour][p_W_ind] += QoS
#                Tx_QoS_Signal_Weight[state_hour][s_W_ind] += QoS
#                Tx_QoS_User_Weight[state_hour][u_W_ind] += QoS
#                Tx_QoS_BW_Weight[state_hour][BW_W_ind] += QoS
                try:
                    Tx_Policy_Output_AvgNoUser[state_hour][intf_threshold_index] += avgnouser/len(list_Active_APs)
#                    Tx_User_Power_Weight[state_hour][p_W_ind] += avgnouser/len(list_Active_APs)
#                    Tx_User_Signal_Weight[state_hour][s_W_ind] += avgnouser/len(list_Active_APs)
#                    Tx_User_User_Weight[state_hour][u_W_ind] += avgnouser/len(list_Active_APs)
#                    Tx_User_BW_Weight[state_hour][BW_W_ind] += avgnouser/len(list_Active_APs)
                except ZeroDivisionError:
                    pass
                        
                        
for intf_threshold_index in range(len(intf_threshold_list)):             
    for state_hour in range(no_periods):
        Tx_Policy_Output_Power[state_hour][intf_threshold_index] = Tx_Policy_Output_Power[state_hour][intf_threshold_index]/len(range(no_sample_paths))
        Tx_Policy_Output_Power_per_user[state_hour][intf_threshold_index] = Tx_Policy_Output_Power_per_user[state_hour][intf_threshold_index]/len(range(no_sample_paths))       
        Tx_Policy_Output_QoS[state_hour][intf_threshold_index] = Tx_Policy_Output_QoS[state_hour][intf_threshold_index]/len(range(no_sample_paths))
        Tx_Policy_Output_AvgNoUser[state_hour][intf_threshold_index] = Tx_Policy_Output_AvgNoUser[state_hour][intf_threshold_index]/len(range(no_sample_paths))
        Tx_Policy_Output_QoS_No_Cluster[state_hour][intf_threshold_index] = Tx_Policy_Output_QoS_No_Cluster[state_hour][intf_threshold_index]/len(range(no_sample_paths))
    
    #print(intf_g_matrix_opt)   
    
#plt.figure(0)
#for intf_threshold_index in range(len(intf_threshold_list)):
#    plt.plot([Tx_Policy_Output_Power[state_hour][intf_threshold_index] for state_hour in range(no_periods)], [1/Tx_Policy_Output_Power_per_user[state_hour][intf_threshold_index] for state_hour in range(no_periods)])
#    plt.xlabel('Number of Active RRHs')
#    plt.ylabel('Number of users per RRH')
#    plt.title('Number of users per RRH at distance thresholds 10,20,30,40,50,60')
#    plt.legend(["Thr="+str(i) for i in intf_threshold_list], loc=0)
    #plt.savefig('PowerIntProg29072015.pdf', bbox_inches='tight')
    #plt.savefig('PowerIntProg29072015.eps', bbox_inches='tight')
    
plt.figure(0)
for intf_threshold_index in range(len(intf_threshold_list)):
    plt.plot([No_users[state_hour] for state_hour in range(no_periods)],[Tx_Policy_Output_Power[state_hour][intf_threshold_index]  for state_hour in range(no_periods)])
    plt.xlabel('Number of Active RRHs')
    plt.ylabel('Number of users per RRH')
    plt.title('Number of users per RRH at distance thresholds 10,20,30,40,50,60')
    plt.legend(["Thr="+str(i) for i in intf_threshold_list], loc=0)
    #plt.savefig('PowerIntProg29072015.pdf', bbox_inches='tight')
    #plt.savefig('PowerIntProg29072015.eps', bbox_inches='tight')

plt.figure(1)    
for intf_threshold_index in range(len(intf_threshold_list)):
    plt.plot([Tx_Policy_Output_Power[state_hour][intf_threshold_index] for state_hour in range(no_periods)],[Tx_Policy_Output_QoS[state_hour][intf_threshold_index] for state_hour in range(no_periods)])
    #plt.axis([0, 26, 1, 29])    
    plt.xlabel('Number of Active RRHs')
    plt.ylabel('QoS (SIR) per User')
    plt.title('QoS versus Active RRHs at distance thresholds 10 to 60')
    
    #plt.savefig('QoSIntProg.pdf', bbox_inches='tight')
plt.gca().set_color_cycle(None)    
plt.figure(1) 
#plt.hold(True)
for intf_threshold_index in range(len(intf_threshold_list)):
    plt.plot([Tx_Policy_Output_Power[state_hour][intf_threshold_index] for state_hour in range(no_periods)],[Tx_Policy_Output_QoS_No_Cluster[state_hour][intf_threshold_index] for state_hour in range(no_periods)],'--')
    #plt.xlabel('Number of Active RRHs')
    #plt.ylabel('QoS (SIR) per User')
    #plt.title('QoS versus Time at distance thresholds 10,20,30,40,50,60')
    #plt.legend(["Thr = "+str(i) for i in intf_threshold_list]+["No Cluster"], loc=0)
    #plt.legend(["Threshold = "+str(i) for i in intf_threshold_list], loc=8)
    #plt.savefig('QoSIntProg29072015.pdf', bbox_inches='tight')
    #plt.savefig('QoSIntProg29072015.eps', bbox_inches='tight')
marker_list = 'os*<>' 
plt.figure(3) 
for intf_threshold_index in range(len(intf_threshold_list)):
    plt.plot(No_users,[Tx_Policy_Output_AvgNoUser[state_hour][intf_threshold_index] for state_hour in range(no_periods)],marker = marker_list[intf_threshold_index])
    plt.xlabel('Number of Users')
    plt.ylabel('Average Number of Users per RRH')
    plt.title('Average Number of Users per RRH distance thresholds 10 to 60')
    plt.legend(["Thr= "+str(i) for i in intf_threshold_list], loc=0)
    plt.savefig('AvgNoIntProg11092015.pdf', bbox_inches='tight')
    plt.savefig('AvgNoIntProg11092015.eps', bbox_inches='tight')
    
plt.figure(4)    
for intf_threshold_index in range(len(intf_threshold_list)):
    plt.plot(No_users,[Tx_Policy_Output_QoS[state_hour][intf_threshold_index] for state_hour  in range(no_periods)])
    plt.xlabel('Number of Users per RRH')
    plt.ylabel('QoS (SIR)')
    plt.title('QoS versus Number of Users per RRH at distance thresholds 10 to 60')
    #plt.legend(["Threshold = "+str(i) for i in intf_threshold_list], loc=8)
    #plt.savefig('QoSIntProg.pdf', bbox_inches='tight')
plt.gca().set_color_cycle(None)     
#plt.figure(5)    
for intf_threshold_index in range(len(intf_threshold_list)):
    plt.plot(No_users,[Tx_Policy_Output_QoS_No_Cluster[state_hour][intf_threshold_index] for state_hour  in range(no_periods)],'-.')
    #plt.xlabel('Time')
    #plt.ylabel('QoS (SIR)')
    #plt.title('QoS versus Time at distance thresholds 10,20,30,40,50,60')
    #plt.legend(["Thr = "+str(i) for i in intf_threshold_list]+["No Cluster"], loc=0)
    #plt.savefig('QoSvsNoUsersIntProg29072015.pdf', bbox_inches='tight')
    #plt.savefig('QoSvsNoUsersIntProg29072015.eps', bbox_inches='tight')
    


plt.figure(6)    
for intf_threshold_index in range(len(intf_threshold_list)):
    plt.plot(No_users,[Tx_Policy_Output_QoS[state_hour][intf_threshold_index]/2-Tx_Policy_Output_AvgNoUser[state_hour][intf_threshold_index]/16 for state_hour  in range(no_periods)])
    plt.xlabel('Number of Users per RRH')
    plt.ylabel('Overall QoS')
    plt.title('Overall QoS versus Number of Users per RRH at distance thresholds 10 to 60')
    #plt.legend(["Threshold = "+str(i) for i in intf_threshold_list], loc=0)
    #plt.savefig('QoSIntProg.pdf', bbox_inches='tight')
plt.gca().set_color_cycle(None)     
#plt.figure(7)    
for intf_threshold_index in range(len(intf_threshold_list)):
    plt.plot(No_users,[Tx_Policy_Output_QoS_No_Cluster[state_hour][intf_threshold_index]/2-Tx_Policy_Output_AvgNoUser[state_hour][intf_threshold_index]/16 for state_hour  in range(no_periods)],'-.')
    #plt.xlabel('Time')
    #plt.ylabel('Overall QoS')
    #plt.title('QoS versus Time at distance thresholds 10,20,30,40,50,60')
    #plt.legend(["Thr = "+str(i) for i in intf_threshold_list]+["No Cluster"], loc=0)
    #plt.savefig('OverallQoSIntProg29072015.pdf', bbox_inches='tight')
    #plt.savefig('OverallQoSIntProg29072015.eps', bbox_inches='tight')
    
plt.figure(8)    
for intf_threshold_index in range(len(intf_threshold_list)):
    plt.plot([1/Tx_Policy_Output_Power_per_user[state_hour][intf_threshold_index] for state_hour in range(no_periods)],[Tx_Policy_Output_AvgNoUser[state_hour][intf_threshold_index] for state_hour  in range(no_periods)])
    plt.xlabel('Number of Users per RRH')
    plt.ylabel('QoS (SIR)')
    plt.title('QoS versus Time at distance thresholds 10,20,30,40,50,60')
    plt.legend(["Threshold = "+str(i) for i in intf_threshold_list], loc=8)
    #plt.savefig('QoSIntProg.pdf', bbox_inches='tight')
    











































    
#plt.figure(6)    
#for state_hour in range(no_periods-1):
#    plt.plot([Tx_Policy_Output_Power_per_user[state_hour][intf_threshold_index] for intf_threshold_index in range(len(intf_threshold_list))],[Tx_Policy_Output_QoS[state_hour][intf_threshold_index] for intf_threshold_index in range(len(intf_threshold_list))])
#    plt.xlabel('Power per User')
#    plt.ylabel('QoS (SIR)')
#    plt.title('QoS versus Time at distance thresholds 10,20,30,40,50,60')
#    #plt.savefig('QoSIntProg.pdf', bbox_inches='tight')
#    
#plt.figure(7)    
#for state_hour in range(no_periods-1):
#    plt.plot([Tx_Policy_Output_Power_per_user[state_hour][intf_threshold_index] for intf_threshold_index in range(len(intf_threshold_list))],[Tx_Policy_Output_QoS_No_Cluster[state_hour][intf_threshold_index] for intf_threshold_index in range(len(intf_threshold_list))],'-.')
#    plt.xlabel('Power per User')
#    plt.ylabel('QoS (SIR)')
#    plt.title('QoS versus Time at distance thresholds 10,20,30,40,50,60')
#    #plt.savefig('QoSIntProg.pdf', bbox_inches='tight')
#
#my_dict ={}        
#for i in [3]:
#    for j in range(len(Tx_Policy_Output_Power_per_user[i])):
#        my_dict[Tx_Policy_Output_Power_per_user[i][j]] = Tx_Policy_Output_QoS[i][j]
#        
#my_dict_2 ={}        
#for i in [3]:
#    for j in range(len(Tx_Policy_Output_Power_per_user[i])):
#        my_dict_2[Tx_Policy_Output_Power_per_user[i][j]] = Tx_Policy_Output_QoS_No_Cluster[i][j]
#        
#        
#x = sorted(my_dict.items(), key=operator.itemgetter(0))
#y = [c[0] for c in x]
#z = [c[1] for c in x]
#
#plt.plot(y,z)
#
#x_2 = sorted(my_dict_2.items(), key=operator.itemgetter(0))
#y_2 = [c[0] for c in x_2]
#z_2 = [c[1] for c in x_2]
#
#plt.plot(y_2,z_2,'r')

#    plt.plot(Tx_Power_Signal_Weight)
#    plt.plot(Tx_QoS_Signal_Weight)
#    plt.plot(Tx_User_Signal_Weight)
#    
#    for state_hour in [0,8,14,20]:
#        plt.plot(s_W_list,Tx_Power_Signal_Weight[state_hour])
#        
#    for state_hour in [0,8,14,20]:
#        plt.plot(s_W_list,Tx_QoS_Signal_Weight[state_hour])
#        
#    for state_hour in [0,8,14,20]:
#        plt.plot(s_W_list,Tx_User_Signal_Weight[state_hour])

        
        
    
   
    #plt.show(True)
        
                    
                    
   #print(intf_g_obj.intf_g)             
   #print(intf_g_obj.intf_g_power)                 
                

        
    