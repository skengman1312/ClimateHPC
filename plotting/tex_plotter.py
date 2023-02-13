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
    res = pd.read_csv("../utils/HPC_final.csv", na_values="-", index_col="Numberofprocesses")
    res.drop(["Serial"], axis = 1, inplace=True)
    print(res.columns)
    ypoints = np.array([34.635 / i for i in range(2, 27)])
    xpoints = np.array([i for i in range(2, 27)])

    ypoints = res["Expcted"]
    res.drop(["Expcted"], axis = 1, inplace=True)

    xpoints = [4,6,9,12,16]

    # res.plot()
    # print(res)
    plt.xticks(range(5, 26, 2))
    plt.plot(xpoints, ypoints, label = "T(serial)/n processes")
    # plt.xlim(5,25)
    # plt.ylim(0,10)
    plt.fill_between(x=xpoints, y1= ypoints-2, y2= ypoints+2, alpha=.1)
    for c in res:
        print(res[c])
        res[c].dropna().plot(marker = ".")
    plt.xlabel("Number of processes")
    plt.ylabel("Wall-time(s)")
    plt.legend()
    plt.savefig("ssh_walltime_plot.png", dpi = 200)
    plt.show()
