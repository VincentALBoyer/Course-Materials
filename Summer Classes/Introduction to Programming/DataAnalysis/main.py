# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import numpy as np
import pandas as pd


def main() -> object:
    csvfile = "tips.csv"
    df = pd.read_csv(csvfile, on_bad_lines="warn")
    print(df.head())
    df.info()
    print(df.isna().sum())
    print(df.describe())



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
