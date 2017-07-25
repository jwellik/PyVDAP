import math

def beta(Na, N, Ta, T):
    beta = (Na - N * (Ta / T)) / math.sqrt(N * (Ta / T) * (1 - (Ta / T)))

