import grid2op
from grid2op.Agent import BaseAgent
from stable_baselines3 import PPO
import numpy as np

class PPOAgent(BaseAgent):
    """Custom agent that bridges Grid2Op and Stable-Baselines3"""
    def __init__(self, action_space, model):
        super().__init__(action_space)
        self.model = model
        
    def act(self, observation, reward, done=False):
        # Convert Grid2Op observation to vector for PPO
        obs_vector = observation.to_vect()
        action_vec, _ = self.model.predict(obs_vector, deterministic=False)
        return self.action_space.from_vect(action_vec)

def main():
    # 1. Create Grid2Op environment
    env = grid2op.make("rte_case14_realistic",
                      reward_class=grid2op.Reward.L2RPNReward)
    
    # 2. Initialize PPO model
    model = PPO(
        policy="MlpPolicy",
        env=env,  # Directly using Grid2Op env (no vectorization)
        learning_rate=3e-4,
        n_steps=10,  # Steps per environment update
        batch_size=64,
        n_epochs=2,
        gamma=0.99,
        verbose=1
    )
    
    # 3. Training loop
    print("Starting training...")
    model.learn(total_timesteps=10)
    
    # 4. Save the trained model
    model.save("ppo_grid2op_single")
    print("Model saved")
    
    # 5. Evaluation
    agent = PPOAgent(env.action_space, model)
    obs = env.reset()
    done = False
    
    while not done:
        action = agent.act(obs)
        obs, reward, done, _ = env.step(action)
        print(f"Step reward: {reward:.2f}")
    
    env.close()

if __name__ == "__main__":
    main()