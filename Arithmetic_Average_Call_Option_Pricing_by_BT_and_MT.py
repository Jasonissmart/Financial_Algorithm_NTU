import numpy as np


class Node:
    def __init__(self, sPrice):
        self.sPrice = sPrice
        self.possibleSaveList = []
        self.possibleCallList = []


class BinomialTree:
    '''
    Build binomial tree to calculate european and american lookback put option price.
    '''
    def __init__(self, paramsList, Savg_spaced_way='Equally'):
        '''
        Initialize parameters and create an empty tree which will store node information.
            Parameters:
                Savg_spaced_way: 'Equally' or 'Logarithmically'
        '''
        self.St, self.K, self.r, self.q, self.sigma, self.t, self.deltaT, self.M, self.n, self.S_avet = paramsList
        self.M = int(self.M)
        self.n = int(self.n)
        self.T = self.deltaT + self.t

        self.previous_n = self.n * self.T / self.deltaT + 1 - self.n
        self.previous_n = int(self.previous_n)
        
        self.u = np.exp(self.sigma * ((self.deltaT / self.n) ** 0.5))
        self.d = np.exp(-self.sigma * ((self.deltaT / self.n) ** 0.5))
        self.p = (np.exp((self.r - self.q) * self.deltaT / self.n) - self.d) / (self.u - self.d)

        self.tree = np.zeros((self.n + 1, self.n + 1), dtype=object)

        self.Savg_spaced_way = Savg_spaced_way

    def calc_stockPice(self):
        '''
        Calculate stock price for every node on the tree.
        '''
        u = self.u
        d = self.d

        # calculate stock price
        for ithCol in range(0, self.n + 1):
            for jthRow in range(ithCol + 1):
                self.tree[jthRow, ithCol] = Node(self.St * (u ** (ithCol - jthRow)) * (d ** jthRow))

    def record_possible_Save(self):
        self.tree[0, 0].possibleSaveList.append(self.S_avet)
        for ithCol in range(1, self.n + 1):
            for jthRow in range(ithCol + 1):
                currentNode = self.tree[jthRow, ithCol]

                A_max = (self.S_avet * self.previous_n
                        + self.St * self.u * ((1 - self.u ** (ithCol - jthRow)) / (1 - self.u))
                        + self.St * (self.u ** (ithCol - jthRow)) * self.d * ((1 - self.d ** jthRow) / (1 - self.d)))
                A_max /= (self.previous_n + ithCol)

                A_min = (self.S_avet * self.previous_n
                        + self.St * self.d * ((1 - self.d ** jthRow) / (1 - self.d))
                        + self.St * (self.d ** jthRow) * self.u * ((1 - self.u ** (ithCol - jthRow)) / (1 - self.u)))
                A_min /= (self.previous_n + ithCol)

                if((jthRow == 0) or (jthRow == ithCol)):
                    currentNode.possibleSaveList.append(A_max)
                else:
                    if(self.Savg_spaced_way == 'Linear'):
                        for k in range(0, self.M + 1):
                            A_i_j_k = ((self.M - k) / self.M) * A_max + (k / self.M) * A_min
                            currentNode.possibleSaveList.append(A_i_j_k)
                    elif(self.Savg_spaced_way == 'Logarithmically'):
                        for k in range(0, self.M + 1):
                            A_i_j_k = np.exp(((self.M - k) / self.M) * np.log(A_max) + (k / self.M) * np.log(A_min))
                            currentNode.possibleSaveList.append(A_i_j_k)


    def calc_terminal_node_Payoff(self):
        lastColIndex = self.n
        for jthRow in range(lastColIndex + 1):
            # Map each element in possibleSmaxList to possibleCallList
            currentNode = self.tree[jthRow, lastColIndex]
            for i in range(len(currentNode.possibleSaveList)):
                currentNode.possibleCallList.append(max(currentNode.possibleSaveList[i] - self.K, 0))

    def binary_search(self, array, target, low, high):
        '''
        Using binary search to find the index that target will larger or equal to array[index]
        and smaller than array[index - 1].
        '''
        if(high >= low):
            mid = low + (high - low) // 2

            if((abs(array[mid] - target) < 10 ** -8) or ((array[mid] < target) and (array[mid - 1] > target))):
                return mid

            elif(array[mid] < target):
                return self.binary_search(array, target, low, mid - 1)
            elif(array[mid] > target):
                return self.binary_search(array, target, mid + 1, high)

    def linear_interpolation_search(self, array, target, low, high):
        '''
        Using linear interpolation search to find the index that target will larger or equal to array[index]
        and smaller than array[index - 1].
        '''
        if high >= low:
            if high == low:
                return low
            
            index = int(((array[low] - target) * high + (target - array[high]) * low) / (array[low] - array[high]))

            counts = 0
            while(~((abs(array[index] - target) < 10 ** -8) or ((array[index] < target) and (array[index - 1] > target)))):
                if((abs(array[index] - target) < 10 ** -8) or ((array[index] < target) and (array[index - 1] > target))):
                    return index
                else:
                    if(counts < 0):
                        counts -= 1
                    else:
                        counts += 1
                    counts = -counts
                    index = index + counts
            return index

    def backward_induction(self, type='European', searchWay='Sequential search'):
        p = self.p

        for ithCol in range(self.n - 1, -1, -1):
            for jthRow in range(ithCol + 1):
                # Process each node
                currentNode = self.tree[jthRow, ithCol]
                upperChildNode = self.tree[jthRow, ithCol + 1]
                lowerChildNode = self.tree[jthRow + 1, ithCol + 1]
                call_upperChildNode = -1
                call_lowerChildNode = -1

                # Match every element
                for k in range(len(currentNode.possibleSaveList)):
                    # The upper child
                    Au = (((ithCol + self.previous_n) * currentNode.possibleSaveList[k] + upperChildNode.sPrice)
                         / (self.previous_n + ithCol + 1))

                    j = 0
                    if(searchWay == 'Sequential search'):
                        while(upperChildNode.possibleSaveList[j] - Au > 10**-8):
                            j += 1
                    elif(searchWay == 'Binary search'):
                        j = self.binary_search(upperChildNode.possibleSaveList, Au, 0, len(upperChildNode.possibleSaveList) - 1)
                    elif(searchWay == 'Linear interpolation search'):
                        j = self.linear_interpolation_search(upperChildNode.possibleSaveList, Au, 0, len(upperChildNode.possibleSaveList) - 1)

                    if(abs(upperChildNode.possibleSaveList[j] - Au) <= 10**-8):
                        call_upperChildNode = upperChildNode.possibleCallList[j]
                    elif(upperChildNode.possibleSaveList[j] < Au):
                        # Inner interploation
                        Wu = (upperChildNode.possibleSaveList[j - 1] - Au) / (upperChildNode.possibleSaveList[j - 1] - upperChildNode.possibleSaveList[j])
                        call_upperChildNode = Wu * (upperChildNode.possibleCallList[j]) + (1 - Wu) * (upperChildNode.possibleCallList[j - 1])
                    
                    # The lower child
                    Ad = (((ithCol + self.previous_n) * currentNode.possibleSaveList[k] + lowerChildNode.sPrice)
                         / (self.previous_n + ithCol + 1))

                    j = 0
                    if(searchWay == 'Sequential search'):
                        while(lowerChildNode.possibleSaveList[j] - Ad > 10**-8):
                            j += 1
                    elif(searchWay == 'Binary search'):
                        j = self.binary_search(lowerChildNode.possibleSaveList, Ad, 0, len(lowerChildNode.possibleSaveList) - 1)
                    elif(searchWay == 'Linear interpolation search'):
                        j = self.linear_interpolation_search(lowerChildNode.possibleSaveList, Ad, 0, len(lowerChildNode.possibleSaveList) - 1)

                    if(abs(lowerChildNode.possibleSaveList[j] - Ad) <= 10**-8):
                        call_lowerChildNode = lowerChildNode.possibleCallList[j]
                    elif(lowerChildNode.possibleSaveList[j] < Ad):
                        # Inner interploation
                        Wd = (lowerChildNode.possibleSaveList[j - 1] - Ad) / (lowerChildNode.possibleSaveList[j - 1] - lowerChildNode.possibleSaveList[j])
                        call_lowerChildNode = Wd *  (lowerChildNode.possibleCallList[j])+ (1 - Wd) * (lowerChildNode.possibleCallList[j - 1])

                    currentNodeCall = np.exp((-1) * self.r * (self.deltaT / self.n)) * (p * call_upperChildNode + (1 - p) * call_lowerChildNode)

                    if(type == 'American'):
                        # Early exercise or not
                        if(currentNodeCall < (currentNode.possibleSaveList[k] - self.K)):
                            currentNodeCall = (currentNode.possibleSaveList[k] - self.K)
                    currentNode.possibleCallList.append(currentNodeCall)

    def getPrice(self):
        return np.round(self.tree[0, 0].possibleCallList[0], 4)


def monteCarlo_calc():
    '''
    Calculate call and put pirce by Monte Carlo simulation.
        Parameters:
            paramsList: list
        Outputs:
            c: float -> call price
            p: float -> put price
    '''
    print('Monte Carlo simulation')

    # St, K, r, q, σ, t, T-t, n, S_avet, simCnt, repCnt
    paramsList = [float(i) for i in input('# St, K, r, q, σ, t, T-t, n, S_avet, simCnt, repCnt = ').split()]
    St, K, r, q, sigma, t, deltaT, n, S_avet, simCnt, repCnt = paramsList
    n, simCnt, repCnt = int(n), int(simCnt), int(repCnt)
    previous_n_to_timet = n * (deltaT + t) / deltaT - n + 1
    total_n = previous_n_to_timet + n

    def drawSample(St, K, r, q, sigma, t, deltaT, n, S_avet, simCnt):
        mean = np.log(St) + ((r - q - sigma ** 2 / 2)*(deltaT / n))
        std = sigma*(((deltaT) / n) ** 0.5)
        lnSt_list = np.random.normal(loc=mean, scale=std, size=int(simCnt))

        St_sum_list = np.add(np.multiply(S_avet, previous_n_to_timet), np.exp(lnSt_list))
        for i in range(n-1):
            mean = lnSt_list + ((r - q - sigma ** 2 / 2) * (deltaT / n))
            lnSt_list = np.random.normal(loc=mean, scale=std)
            St_sum_list = np.add(np.exp(lnSt_list), St_sum_list)
        
        Save_list = St_sum_list / total_n
        payoff_list = np.maximum(Save_list - K, 0) * np.exp(-r * deltaT)
        return np.mean(payoff_list)

    resultList = []
    for i in range(int(repCnt)):
        p = drawSample(St, K, r, q, sigma, t, deltaT, n, S_avet, simCnt)
        resultList.append(p)
    
    print(f'  Average call price upper bound: {np.round(np.mean(resultList) + 2*np.std(resultList), 4)}')
    print(f'  Average call priceclower bound: {np.round(np.mean(resultList) - 2*np.std(resultList), 4)}', '\n')


def calc_averageCallPrice(
    paramsList,
    retrun_European_price=False,
    savg_spaced_way='Linear',
    searchWay='Sequential search'):
    '''
    Calculate european and american average call option price and print them out.
    For each kind of call, it will 
    1. initialize a BinomialTree class
    2. calculate stock price for every node on the tree
    3. record possible average stock price for every node on the tree
    4. calculate the terminal nodes' payoff
    5. implement backward induction for every node one the tree
    6. print out the payoff at tree[0, 0]
        Parameters:
            paramsList: list -> includes St, K, r, q, σ, t, T-t, M, n, S_avet
            retrun_European_price: boolean -> if true, return european price and 
            stop calculating american option price.
            savg_spaced_way: str -> 'Linear' or 'Logarithmically' equally space possible average stock price
            searchWay: str -> 'Sequential search', 'Binary search' or 'Linear interpolation search'
        Output:
            print European & American average call price
    '''
    print('Binomial Tree method')
    
    # For European
    binomialTree = BinomialTree(paramsList, Savg_spaced_way=savg_spaced_way)
    binomialTree.calc_stockPice()
    binomialTree.record_possible_Save()
    binomialTree.calc_terminal_node_Payoff()
    binomialTree.backward_induction(searchWay)
    if(retrun_European_price):
        return binomialTree.getPrice()
    print('European Average Call price: ', binomialTree.getPrice())
    
    # For American
    binomialTree = BinomialTree(paramsList, Savg_spaced_way=savg_spaced_way)
    binomialTree.calc_stockPice()
    binomialTree.record_possible_Save()
    binomialTree.calc_terminal_node_Payoff()
    binomialTree.backward_induction('American', searchWay)
    print('American Average Call price: ', binomialTree.getPrice(), '\n')
    
    
if __name__ == "__main__":
    # test for t = 0 
    # St, K, r, q, σ, t, T-t, M, n, S_avet = 50 50 0.1 0.05 0.8 0 0.25 100 100 50
    paramsList = [float(i) for i in input('# St, K, r, q, σ, t, T-t, M, n, S_avet = ').split()]
    calc_averageCallPrice(paramsList)

    # test for t = 0.25
    # St, K, r, q, σ, t, T-t, M, n, S_avet = 50 50 0.1 0.05 0.8 0.25 0.25 100 100 50
    paramsList = [float(i) for i in input('# St, K, r, q, σ, t, T-t, M, n, S_avet = ').split()]
    calc_averageCallPrice(paramsList)

    # St, K, r, q, σ, t, T-t, n, S_avet, simCnt, repCnt = 50 50 0.1 0.05 0.8 0 0.25 100 50 10000 20
    monteCarlo_calc()
