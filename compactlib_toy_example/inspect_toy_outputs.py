from pathlib import Path
import pandas as pd

KEY = ["ModifiedPeptide", "PrecursorCharge"]
SHOW = [
    "ModifiedPeptide", "PrecursorCharge",
    "FragmentType", "FragmentSeriesNumber", "FragmentCharge",
    "ProductMz", "LibraryIntensity", "PredictionModel"
]

def show(path):
    print("\n" + "=" * 90)
    print(path)
    print("=" * 90)
    df = pd.read_csv(path, sep="\t")
    df = df.sort_values(KEY + ["LibraryIntensity", "FragmentType", "FragmentSeriesNumber", "FragmentCharge", "ProductMz"],
                        ascending=[True, True, False, True, True, True, True])
    print(df[SHOW].to_string(index=False))

def show_summary(path):
    s = Path(path).with_suffix(".summary.tsv")
    if s.exists():
        print("\nsummary:", s)
        print(pd.read_csv(s, sep="\t").to_string(index=False))

for p in [
    "toy_out/toy_prosit_max3.tsv",
    "toy_out/toy_union6.tsv",
    "toy_out/toy_consensus3.tsv",
    "toy_out/toy_prosit_random3.tsv",
    "toy_out/toy_prosit_max3_reverse.tsv",
]:
    show(p)
    show_summary(p)
