import matplotlib.pyplot as plt
import pandas as pd
# read csv from current working directory
# format of csv:
# size | blocking | non-blocking | collective | seq |
# no header file
def read_csv():
    X = pd.read_csv("data.csv", ' ', header=None, usecols=[0, 1, 2, 3, 4]).to_numpy()
    
    plt.plot(X[:, 0], X[:, 1], label = 'Blocking P2P')
    plt.plot(X[:, 0], X[:, 2], label = 'Non-Blocking P2P')
    plt.plot(X[:, 0], X[:, 3], label = 'Collective')
    plt.plot(X[:, 0], X[:, 4], label = 'Sequential')
    plt.xlabel('N')
    plt.ylabel('Exec. Time(sec)')
    plt.legend()
    plt.savefig('data.png')

read_csv()
