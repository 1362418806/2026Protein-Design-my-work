from pathlib import Path
import random
import pandas as pd

AA = "ACDEFGHIKLMNPQRSTVWY"
random.seed(1)
base = "M" + "VSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"[:238]
out = Path("/mnt/data/synbio_gfp_v3/tests/synthetic_data")
out.mkdir(parents=True, exist_ok=True)
(out / "AAseqs of 5 GFP proteins.txt").write_text(f">sfGFP\n{base}\n>avGFP\n{base}\n", encoding="utf-8")
pd.DataFrame({"Sequence": [base[::-1]], "Name": ["dummy"]}).to_csv(out / "Exclusion_List.csv", index=False)
rows=[]
for i in range(250):
    seq = list(base)
    k = random.randint(0, 8)
    for p in random.sample(range(1, len(seq)), k):
        seq[p] = random.choice(AA)
    brightness = 4.2 - 0.13*k + random.gauss(0,0.15)
    rows.append({"Sequence": "".join(seq), "Brightness": brightness})
pd.DataFrame(rows).to_excel(out / "GFP_data.xlsx", index=False)
print(out)
