from grid2op.Agent import BaseAgent, RandomAgent, DoNothingAgent

class LoadSheddingAgent(BaseAgent):
    def __init__(self, action_space):
        # define here the constructor of your agent
        # here we say our agent needs "something_else" and "and_another_something"
        # to be built just to demonstrate it does not cause any problem to extend the
        # construction of the base class BaseAgent that only takes "action_space" as a constructor
        BaseAgent.__init__(self, action_space)
        self.do_nothing = self.action_space({})

    def act(self, obs, reward, done=False):
        if ((obs.current_step >= 2) & (obs.current_step < 3)):
            new_load_p = obs.load_p * 1.1
            new_load_q = obs.load_q * 1.1
            
            # this is the only method you need to implement
            # it takes an observation obs (and a reward and a flag)
            # and should return a valid action
            dictionary_describing_the_action = {
                    "injection": {
                        "load_p": new_load_p,
                        "load_q": new_load_q
                    }
                }  # this can be anything you want that grid2op understands
            
            my_action = self.action_space(dictionary_describing_the_action)
        else:
            my_action = self.do_nothing
        return my_action