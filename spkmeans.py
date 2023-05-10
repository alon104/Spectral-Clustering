import numpy as np
import sys
from  mykmeanssp import wam, ddg, gl, jacobi, spk
import math

np.random.seed(0)

def kmeans_pp(K, dataPoints):
    indices = [i for i in range(len(dataPoints))]
    centroidsIndices = initializeCentroids(dataPoints,K,indices)
    centroids = []
    for i in range(K):
        centroids.append(dataPoints[indices.index(centroidsIndices[i])])
    final_centroid = spk(K, 300, 0, len(centroids[0]), len(dataPoints), centroids, dataPoints)
    for i in range(K):
        print(centroidsIndices[i],end='')
        if i!=K-1:
            print(",", end='')
    print()
    for i in range(K):
        for j in range(len(centroids[0])):
            print("{:.4f}".format(final_centroid[i][j]), end='')
            if j!=len(centroids[0])-1:
                print(",", end='')
        print()

def eigengap_heuristic(A):
    n = A.shape[0]
    eigen_matrix = A
    sorted_diagonal = np.sort(np.diag(eigen_matrix)) # eigen values sorted
    delta_array = np.zeros(n-1)
    for i in range(n-1):
        delta_array[i] = abs(sorted_diagonal[i] - sorted_diagonal[i+1])
    delta_array = delta_array[:math.floor(n/2)]
    return np.argmax(delta_array) + 1

def initializeCentroids(dataPoints,K,indices):
    N = len(dataPoints)
    centroidIndices = []
    centroidIndices.append(int(np.random.choice(indices)))
    while (len(centroidIndices) < K):
        distances = [0 for i in range(N)]
        for vector in dataPoints:
            minDistance = float('inf')
            for ind in centroidIndices:
                distance = calc_distance(vector,dataPoints[indices.index(ind)])
                if distance < minDistance:
                    minDistance = distance
            distances[dataPoints.index(vector)] = minDistance
        total = np.sum(distances)
        probabilty = [dist/total for dist in distances]
        chosen = (np.random.choice(indices,p=probabilty))
        centroidIndices.append(int(chosen))
    return centroidIndices

def calc_distance(list1,list2):
    assert(len(list1) == len(list2))
    return sum([(list1[ind] - list2[ind]) ** 2 for ind in range(len(list1))])**0.5

def printMatrix(matrix):
    n = len(matrix)
    m = len(matrix[0])
    for i in range(n):
        for j in range(m):
            print("{:.4f}".format(matrix[i][j]), end='')
            if j!=m-1:
                print(",", end='')
        print()

def printJacobi(A,V):
    n = len(A)
    for i in range(n):
        print("{:.4f}".format(A[i][i]), end='')
        if i!=n-1:
            print(",", end='')
    print()
    for i in range(n):
        for j in range(n):
            print("{:.4f}".format(V[i][j]), end='')
            if j!=n-1:
                print(",", end='')
        print()

def main(k, goal, matrix, n, d):
    if goal == 'wam':
        res = wam(matrix, n, d) # res is W
        printMatrix(res)
        return 0
    elif goal == 'ddg':
        res = ddg(matrix, n, d) # res is D
        printMatrix(res)
        return 0
    elif goal == 'gl':
        res = gl(matrix, n, d) # res is L
        printMatrix(res)
        return 0
    elif goal == 'jacobi': 
        A, V = jacobi(matrix, n) # V (eigenVectors) and A is a diagonal of eigenValues
        printJacobi(A, V)
        return 0
    else:
        res = gl(matrix, n, d)
        A, V = jacobi(res, n) 
        A = np.array(A)
        V = np.array(V)
        idx = np.argsort(np.diag(A))  # get indices of sorted diagonal elements
        V = V[:, idx]  # sort columns of V using the indices
        if k == 'heuristic':
            k = eigengap_heuristic(A)
        U = V[:, :k]
        U = U.tolist()
        kmeans_pp(k, U) # need to extract some of the vectors from V, so this isn't exactly V
        
    

if __name__ == "__main__":
    num_args=len(sys.argv)
    i = 2
    if num_args == 4: #K is given
        k = int(sys.argv[1])
        i += 1
    else: 
        k = 'heuristic'
    matrix = []
    with open(sys.argv[i]) as f:
        lines = f.readlines()
    for line in lines:
        asList = line.split(',')
        matrix.append(list(float(asList[i]) for i in range(len(asList))))
    main(k, sys.argv[i-1], matrix, len(matrix), len(matrix[0]))