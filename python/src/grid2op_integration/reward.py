import grid2op
from grid2op.Reward import BaseReward
from grid2op.Action import BaseAction
from grid2op.Environment import BaseEnv


class HADRECReward(BaseReward):
    def __init__(self, preFaultActionPenalty, actionPenalty, invalidActionPenalty, logger=None):
        self.preFaultActionPenalty = preFaultActionPenalty
        self.actionPenalty = actionPenalty
        self.invalidActionPenalty = invalidActionPenalty
        super().__init__(logger)

    def __call__(self,
                action: BaseAction,
                env: BaseEnv,
                has_error: bool,
                is_done: bool,
                is_illegal: bool,
                is_ambiguous: bool) -> float:
        # only method really required.
        # called at each step to compute the reward.
        # this is where you need to code the "formula" of your reward
        return 1

        # -----------check whether the action is applied before fault_start_time
        # if self.current_simu_time < self.fault_start_time:
        #     reward -= self.preFaultActionPenalty

        # # ---------check whether the action is an invalid action, need check qiuhua for his implementation here
        # #if ( remain_load[iact] + actionPyAry[iact]) > 0.01: # if it is a valid action
        # if ( remain_load[iact] ) > 0.01: # if it is a valid action    
        #     reward += self.actionPenalty * actionPyAry[iact] * self.action_buses_pvalue_pu[iact]
        # else:
        #     reward -= self.invalidActionPenalty 

        # #--------compute the voltage deviation part for the reward    
        # ob_volt_tmp = self.ob_vals[-1][0:self.nobbus]
        # fault_end_time = self.fault_start_time + self.fault_duration_time
        
        # for ivoltob in range (self.nobbus):
        #     if self.fault_start_time <= self.current_simu_time < fault_end_time:
        #         volt_penalty = 0.0
                
        #     elif (fault_end_time) <= self.current_simu_time < (fault_end_time + 0.33) :
        #         volt_penalty = min(0.0, ob_volt_tmp[ivoltob]-0.7)
                
        #     elif (self.current_simu_time < (fault_end_time + 0.5)) and (self.current_simu_time >= (fault_end_time + 0.33) ):
        #         volt_penalty = min(0.0, ob_volt_tmp[ivoltob]-0.8)
                
        #     elif (self.current_simu_time < (fault_end_time + 1.5)) and (self.current_simu_time >= (fault_end_time + 0.5) ):
        #         volt_penalty = min(0.0, ob_volt_tmp[ivoltob]-0.9)
                
        #     else:
        #         volt_penalty = min(0.0, ob_volt_tmp[ivoltob]- self.minVoltRecoveryLevel)
                
        #     reward += self.observationWeight * volt_penalty
        
        # volt_rew = reward - invalidact_rew
        
        # # --------- if the voltage could not be back at minVoltRecoveryLevel after fault_end_time + maxVoltRecoveryTime, done
        # unstable_rew = 0.0
        # for ivoltob in range (self.nobbus): # check with Qiuhua for his implementation     
        #     if (self.current_simu_time > (fault_end_time + self.maxVoltRecoveryTime)) and (ob_volt_tmp[ivoltob] < self.minVoltRecoveryLevel):   
        #         reward += self.unstableReward
        #         unstable_rew += self.unstableReward 
        #         done = True 
        #         #print ("------bus %d voltage still not recoverd to %f till %f second after fault clear, simulation done-------"%(self.obs_busIDs[ivoltob], self.minVoltRecoveryLevel, self.maxVoltRecoveryTime))
        #         break
        # #---------------compute reward here-------------------    
        # if done and self.steps_beyond_done is None:
        #     self.steps_beyond_done = 0

    def initialize(self, env: BaseEnv):
        # optional
        # called once, the first time the reward is used
        pass

    def reset(self, env: BaseEnv):
        # optional
        # called by the environment each time it is "reset"
        pass

    def close(self):
        # optional called once when the environment is deleted
        pass
