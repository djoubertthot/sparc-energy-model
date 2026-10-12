from radial_model import process_directory

def main():
    res = process_directory("data/rotmod/")
    res.to_csv("results/rms_all_galaxies.csv", index=False)
    print(res["RMS"].describe())
    print("\nFull RMS CSV saved to results/rms_all_galaxies.csv")

if __name__ == "__main__":
    main()
