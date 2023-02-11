import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

res = pd.DataFrame([[8.932, 7.586, None],
                    [4.682, None, 4.4356],
                    [None, 3.732, None],
                    [6.372, None, None],
                    [4.503, None, 5.726],
                    [None, None, None]
                    ], columns=["group of 2", "group of 3", "group of 4"], index=[6, 8, 9, 10, 12, 15])

if __name__ == "__main__":
    ypoints = np.array([36.238 / i for i in range(5, 12)])
    xpoints = np.array([i for i in range(5, 17)])
    # res.plot()
    # print(res)
    plt.plot(xpoints, ypoints)
    for c in res:
        print(res[c])
        res[c].dropna().plot()
    plt.xlabel("Number of processes")
    plt.ylabel("Wall-time")
    plt.legend()

    plt.show()
