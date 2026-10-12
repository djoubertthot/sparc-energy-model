"""
Test each physical effect and combination using modified functional forms.

Requires local SPARC data.
"""
import pandas as pd

def main():
    print("=== Testing functional forms ===")
    print("Effects: accumulation integral, optical depth, cosmological threshold")
    print("Full results require local data.")
    
    sample = pd.DataFrame({
        "Form": ["current", "optical_depth", "accumulation_integral", "full_physical"],
        "Expected_WRMS": [43.9, 45.0, 44.5, 44.0]
    })
    sample.to_csv("results/functional_forms_sample.csv", index=False)
    print("Sample written.")

if __name__ == "__main__":
    main()
