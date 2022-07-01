import csv
import math
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--algorithm', type = int, default = 1, help = '1 or 2')
parser.add_argument('--N', type = int, default = 2, help = 'should be integer >= 2')
parser.add_argument('--visualization', type = int, default = 1, help = '0: False and 1: True')
parser.add_argument('--earlystopping', type = int, default = 0, help = '0: False and 1: True')
parser.add_argument('--molecular', type = int, default = 2)
parser.add_argument('--denominator', type = int, default = 1)
args = parser.parse_args()

EARLY = True if args.earlystopping == 1 else False
ALPHA = args.molecular/args.denominator

'''
    find_vertex1: A function to compute a vertex x of a 0/1 polytope that maximizes the ratio c*(x-x')/L1_norm(x-x')
    Input: Dimension n
           Current pivoting point, x_i
           Vector c 
    Output: Return a tuple of the vertex index and the ratio where vertex_index refers to a specific vertex (e.g. if vertex_index is 3 (n = 4), then the vertext would be (0,1,1,1))
'''

def find_vertex1(n,i,c):
    ratio = None
    vertex_index = i
    for j in range(n,i,-1): # only consider the vertex x_j where j > i since the ratio would be negative if j < i
        tmp_sum = 0
        for k in range(i,j):
            tmp_sum += c[n-k-1]
        tmp_ratio = tmp_sum/(j-i) # calculate the ratio
        if ratio == None:
            ratio = tmp_ratio
            vertex_index = j
        elif tmp_ratio > ratio: # maximize the ratio
            ratio = tmp_ratio
            vertex_index = j
    return vertex_index, ratio

'''
    Algorithm 1: MRA-based Geometric Scaling
    Input: dimension n (in other words, input would be n+1 points of a 0/1 polytope P)
           vector c = (1,2,3,...,n)
    Output: Return a list storing the point that maximize the objective function, total number of iterations (halving + augmenting),
            the number of havling iterations and the number of augmenting iterations.
'''

def Algorithm1(n, c):
    mu = max(c) + 1 # not sure how to choose the value of mu; max(c) is the L-infinite norm of vector c
    i = 0 # current pivoting point, e.g. x_i is the point in R_n whose last i coordinates are equal to 1 and whose other coordinates are equal to 0.
    halving = 0 # number of halving steps
    augmenting = 0 # number of augmenting steps
    while mu >= 1/n:
        next_pivot, ratio = find_vertex1(n,i,c)
        if next_pivot == i or ratio < mu: # halving 
            mu = mu/2
            halving += 1
        else: # augmenting
            i = next_pivot
            augmenting += 1
            if EARLY and i == n: # early stoping
                return [i, halving + augmenting, augmenting, halving]
    return [i, halving + augmenting, augmenting, halving]

'''
    find_vertex2: A function to compute a vertex x of a 0/1 polytope that satisfy c*(x-x') > mu*L1_norm(x-x')
    Input: Dimension n
           Current pivoting vertex, x_i
           Vector c 
           Number mu
    Output: Return the vertex index where vertex_index refers to a specific vertex (e.g. if vertex_index is 3 (n = 4), then the vertext would be (0,1,1,1))
            If the returned value is equal to None, it means that there is no such vertex of the 0/1 polytope.
'''

def find_vertex2(n,i,c,mu):
    vertex_index = None
    for j in range(n,i,-1):
        # another way to calculate the ratio
        # tmp_ratio = (2**n)*((2**(1-i)-2**(1-j))/(j-i))
        tmp_sum = 0
        for k in range(i,j):
            tmp_sum += c[n-k-1]
        tmp_ratio = tmp_sum/(j-i) # calculate the ratio
        if tmp_ratio > mu:
            return j
    return vertex_index


'''
    Algorithm 2: Geometric Scaling via a feasibility test
    Input: dimension n (in other words, input would be n+1 points of a 0/1 polytope P)
           vector c = (1,2,3,...,n)
    Output: Return a list storing the point that maximize the objective function, total number of iterations (halving + augmenting),
            the number of havling iterations and the number of augmenting iterations.
'''

def Algorithm2(n, c):
    mu = max(c) + 1 # not sure how to choose the value of mu; max(c) is the L-infinite norm of vector c
    i = 0 # current pivoting point, e.g. x_i is the point in R_n whose last i coordinates are equal to 1 ans whose other coordinates are equal to 0.
    halving = 0 # # number of halving steps
    augmenting = 0 # number of augmenting steps
    while mu >= 1/n:
        next_pivot = find_vertex2(n,i,c,mu)
        if next_pivot: # augmenting step
            i = next_pivot
            augmenting += 1
            if EARLY and i == n: # early stoping
                return [i, halving + augmenting, augmenting, halving]
        else: # there is no such vertex, halving step
            mu = mu/ALPHA
            halving += 1
    return [i, halving + augmenting, augmenting, halving]


'''
    Run_algorithm: A function to run the algorithm 1 or 2 for n = 2,...,N
    Input: A number 'N' 
           A number 'algo' should be 1 or 2 (Algowithm 1 would be executed if 'algo' == 1 and Algorithm 2 would be executed if 'algo' == 2)
    Output: Return a dictionary storing the number of iterations (total, halving and augmenting) for n = 2, ... , N
'''
def Run_algorithm(N, algo):
    result = {}
    if algo == 1:
        for n in range(2,N+1):
            c = [i for i in range(1,n+1)] # c = (1,2,...,n)
            result[n] = Algorithm1(n,c)
        return result
    elif algo == 2:
        for n in range(2,N+1):
            c = [math.ceil(ALPHA)**i for i in range(1,n+1)] # c = (2,4,...,2^n)
            result[n] = Algorithm2(n,c)
        return result

'''
    Visualization: Plot the number of iterations (total, halving and augmenting) after running algorithm 1 or 2
    Input: A dictionary 'result'
'''

def Visualization(result):
    Dimension = [i for i in range(2,max(result.keys())+1)] # x_label
    total_iterations = [result[key][1] for key in result] # number of augmenting steps and halving steps
    augmenting_iterations = [result[key][2] for key in result] # number of augmenting steps
    halving_iterations = [result[key][3] for key in result] # number of halving steps
    plt.plot(Dimension, total_iterations, label='Total')
    plt.plot(Dimension, augmenting_iterations,label='Augmenting')
    plt.plot(Dimension, halving_iterations,label='Halving')
    plt.xlabel("Dimension")
    plt.ylabel("Number of Iterations")
    plt.legend()
    plt.show()

'''
    Save_result: Save the result in a text file
    Input: a dictionary, result
           a number, ALGO (1 or 2)
           a number, N
'''
def Save_result(result,ALGO,N):
    filename = "Algorithm_{}_N_{}_Early_{}.csv".format(ALGO,N,EARLY)
    fielfnames = ['Dimension','Augmenting','Halving','Total','Early Stopping','Alpha']
    with open(filename,'a') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(fielfnames)
        for key, iterations in result.items():
            writer.writerow([key,iterations[2], iterations[3], iterations[1],EARLY,ALPHA])


def main():
    # INPUT
    ALGO = args.algorithm
    N = args.N
    VISUALIZATION = True if args.visualization == 1 else False
    # Raise Exception
    if ALGO != 1 and ALGO != 2:
        raise Exception("The parameter 'algo' should be 1 or 2.")
    # Run the algorithm for n = 2,...,N
    result = Run_algorithm(N,ALGO)
    print(result)
    # Save the result
    Save_result(result,ALGO=ALGO,N=N)
    # Plot the result if VISUALIZATION == True
    if VISUALIZATION:
        Visualization(result)


if __name__ == '__main__':
    # freeze_support() here if program needs to be frozen
    main()  # execute this only when run directly, not when imported!