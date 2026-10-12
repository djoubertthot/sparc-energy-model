import pandas as pd
from radial_model import process_directory

def main():
    results = process_directory('data/rotmod/')
    results.to_csv('results/radial_rms.csv', index=False)
    print(results.describe())
    print('\nSample of radial results saved to results/radial_rms.csv')

if __name__ == '__main__':
    main()
