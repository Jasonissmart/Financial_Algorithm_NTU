import numpy as np
from Arithmetic_Average_Call_Option_Pricing_by_BT_and_MT import  calc_averageCallPrice
import time

if __name__ == "__main__":
    # test for t = 0 
    # St, K, r, q, σ, t, T-t, M, n, S_avet = 50 50 0.1 0.05 0.8 0 0.25 400 100 50
    paramsList = [float(i) for i in input('# St, K, r, q, σ, t, T-t, M, n, S_avet = ').split()]
    
    start_time = time.time()
    calc_averageCallPrice(paramsList, retrun_European_price=True, searchWay='Sequential search')
    print("--- Sequential search: %s seconds ---\n" % (time.time() - start_time))

    start_time = time.time()
    calc_averageCallPrice(paramsList, retrun_European_price=True, searchWay='Binary search')
    print("--- Binary search: %s seconds ---\n" % (time.time() - start_time))

    start_time = time.time()
    calc_averageCallPrice(paramsList, retrun_European_price=True, searchWay='Linear interpolation search')
    print("--- Linear interpolation search: %s seconds ---\n" % (time.time() - start_time))
