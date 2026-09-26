#!/usr/bin/env python3
"""Read-level reassignment in Design 2 (nine runs) with three denominators: all classified reads, without phiX174 reads,
and without phiX174 reads and host reads assigned to the Shamonda virus S-segment node (taxid 159150), whose record ends
with a repetitive Chinese hamster motif. Also the same for the 17 prebuilt databases (fourteen runs), excluding reads
assigned to 159150 and to Proteus phage VB_PmiS-Isfahan (1969841) in either database."""
import csv
import os
K = "/home/hm/kdcr"; W2 = f"{K}/work_designs"
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
MSL = [f"msl{n}" for n in range(29, 42)]
RUNS9 = ["SRR26352217", "SRR26352207", "SRR26352216", "SRR26352206", "SRR26352195", "SRR26352212", "SRR26352205", "SRR26352208", "SRR26352202"]
RUNS14 = ["SRR26352207", "SRR26352206", "SRR26352204", "SRR26352215", "SRR26352205", "SRR26352203", "SRR26352202",
          "SRR26352217", "SRR26352216", "SRR26352195", "SRR26352184", "SRR26352212", "SRR26352209", "SRR26352208"]
OFFICIAL = ["20201202", "20210517", "20220607", "20220908", "20221209", "20230314", "20230605", "20231009", "20240112",
            "20240605", "20240904", "20241228", "20250402", "20250714", "20251015", "20260226", "20260626"]
HOST = {"159150", "1969841"}


def cls(p):
    return {l.split("\t", 3)[1]: l.split("\t", 3)[2] for l in open(p) if l[0] == "C"}


rows = []
phx = {r: set(l.strip() for l in open(f"{W2}/phix/{r}.phix_ids")) for r in RUNS9}
for a, b in zip(MSL, MSL[1:]):
    n = [0, 0, 0]; d = [0, 0, 0]
    for r in RUNS9:
        da, db = cls(f"{K}/results_fixed/{a}/{r}_output.txt"), cls(f"{K}/results_fixed/{b}/{r}_output.txt")
        for k, v in da.items():
            if k not in db: continue
            ch = v != db[k]; n[0] += 1; d[0] += ch
            if k in phx[r]: continue
            n[1] += 1; d[1] += ch
            if v in HOST or db[k] in HOST: continue
            n[2] += 1; d[2] += ch
    rows.append(["Design 2", f"{a}->{b}"] + [x for pair in zip(n, d) for x in pair] + [round(100 * d[i] / n[i], 3) for i in range(3)])
    print(rows[-1], flush=True)
for a, b in zip(OFFICIAL, OFFICIAL[1:]):
    n = [0, 0, 0]; d = [0, 0, 0]
    for r in RUNS14:
        da, db = cls(f"{W2}/official/cls/{a}/{r}_output.txt"), cls(f"{W2}/official/cls/{b}/{r}_output.txt")
        for k, v in da.items():
            if k not in db: continue
            ch = v != db[k]; n[0] += 1; d[0] += ch
            if v in HOST or db[k] in HOST: continue
            n[2] += 1; d[2] += ch
    rows.append(["prebuilt", f"{a}->{b}", n[0], d[0], "", "", n[2], d[2], round(100 * d[0] / n[0], 3), "", round(100 * d[2] / n[2], 3)])
    print(rows[-1], flush=True)
with open(f"{OUT}/reassignment_denominators.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["databases", "update", "classified_both", "discordant", "classified_both_non_phix", "discordant_non_phix",
                "classified_both_non_phix_non_host_record", "discordant_non_phix_non_host_record", "pct_all", "pct_non_phix", "pct_non_phix_non_host_record"])
    w.writerows(rows)
