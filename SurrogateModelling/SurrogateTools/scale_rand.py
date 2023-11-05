import numpy as np

def scale_rand(lb,ub):
  S = np.random.random_sample(len(lb))
  for i in range(len(lb)):
    c = lb[i]
    d = ub[i]
    S[i] = c + S[i] * (d-c)
  return S
