import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

res = pd.DataFrame([[8.932, 7.586, None],
                    [4.682, None, 4.4356],
                    [None, 3.732, None],
                    [3.833, None, None],
                    [3.22, None, 2.695],
                    [None, 2.588, None],
                    [2.735, 2.197, None],
                    [4.986, None, 4.151],
                    [5.002, None, 2.654]
                    ], columns=["2 groups", "3 groups", "4 groups"], index=[6, 8, 9, 10, 12, 15, 18, 20, 24])

if __name__ == "__main__":
    ypoints = np.array([34.635 / i for i in range(2, 27)])
    xpoints = np.array([i for i in range(2, 27)])
    # res.plot()
    # print(res)
    plt.xticks(range(5, 26, 2))
    plt.plot(xpoints, ypoints, label = "T(serial)/n processes")
    plt.xlim(2,25)
    plt.fill_between(x=xpoints, y1= ypoints-0.4, y2= ypoints+0.4, alpha=.1)
    for c in res:
        print(res[c])
        res[c].dropna().plot(marker = ".")
    plt.xlabel("Number of processes")
    plt.ylabel("Wall-time")
    plt.legend()

    plt.show()
