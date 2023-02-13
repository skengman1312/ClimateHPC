import pandas as pd



if __name__ == "__main__":
    df = pd.read_csv("tdata.csv")
    df["speedup"] = df.iloc[0,-1] / df["computation wall-time(s)"]
    df["efficency"] = df["speedup"]/df["Number of processes"]
    print(df.to_latex(index = False, float_format="%.3f", column_format="c|c|c|c|c|c"))