import numpy as np
from Arithmetic_Average_Call_Option_Pricing_by_BT_and_MT import  calc_averageCallPrice

if __name__ == "__main__":
    # test for equally space Save
    # St, K, r, q, σ, t, T-t, M, n, S_avet = 50 50 0.1 0.05 0.8 0 0.25 100 100 50
    # M could be any number
    paramsList = [float(i) for i in input('# St, K, r, q, σ, t, T-t, M, n, S_avet = ').split()]
    equally_space_european_price_result = []
    for m in range(50, 450, 50):
        paramsList[7] = int(m)
        equally_space_european_price_result.append(calc_averageCallPrice(paramsList, retrun_European_price=True, savg_spaced_way='Equally'))
    print(f'Linearly equally-spaced: {equally_space_european_price_result}')

    # test for Logarithmically space Save
    logarithmically_space_european_price_result = []
    for m in range(50, 450, 50):
        paramsList[7] = int(m)
        logarithmically_space_european_price_result.append(calc_averageCallPrice(paramsList, retrun_European_price=True, savg_spaced_way='Logarithmically'))
    print(f'Logarithmically equally-spaced: {logarithmically_space_european_price_result}')