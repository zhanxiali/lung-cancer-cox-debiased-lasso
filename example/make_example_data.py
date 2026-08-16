"""
Generate a simulated dataset with the same structure as the study data.

Produces:
  example_genotypes.raw   PLINK --recode A format (FID IID PAT MAT SEX PHENOTYPE + SNP columns)
  example_clinical.csv    IID, demographics, treatment, PCs, OS_month, DEAD

No participant data are used. Genotypes are drawn under Hardy-Weinberg equilibrium with
LD blocks; survival times are simulated from an exponential Cox model in which a small
number of variants carry non-null effects, so the pipeline produces a realistic mix of
signal and noise.

Usage:  python make_example_data.py [--n 400] [--snps 77] [--seed 2026]
"""
import argparse, csv, math, random

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--snps", type=int, default=77)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    rng = random.Random(a.seed)

    # SNP metadata: MAF, and LD blocks of 2-4 correlated variants
    mafs = [round(rng.uniform(0.05, 0.48), 3) for _ in range(a.snps)]
    blocks, i = [], 0
    while i < a.snps:
        size = min(rng.choice([1, 1, 2, 3, 4]), a.snps - i)
        blocks.append(list(range(i, i + size)))
        i += size

    alleles = ["A", "C", "G", "T"]
    snp_names, effect_allele = [], []
    for j in range(a.snps):
        ea, oa = rng.sample(alleles, 2)
        snp_names.append(f"rsSIM{j+1:05d}_{ea}")
        effect_allele.append((ea, oa))

    # true log-hazard ratios: 6 non-null variants
    beta = [0.0] * a.snps
    for j in rng.sample(range(a.snps), 6):
        beta[j] = rng.choice([-1, 1]) * rng.uniform(0.08, 0.20)

    geno, clin = [], []
    for k in range(a.n):
        g = [0] * a.snps
        for blk in blocks:
            # shared latent liability induces LD within a block
            z = rng.gauss(0, 1)
            for j in blk:
                p = mafs[j]
                u = 0.75 * z + 0.66 * rng.gauss(0, 1)
                thr_lo = _qnorm((1 - p) ** 2)
                thr_hi = _qnorm(1 - p * p)
                g[j] = 0 if u < thr_lo else (1 if u < thr_hi else 2)
        geno.append(g)

        age = max(35, min(90, rng.gauss(65.6, 10.6)))
        sex = 1 if rng.random() < 0.47 else 2
        smok = rng.choices([0, 1, 2], weights=[0.12, 0.57, 0.31])[0]
        stage = 1 if rng.random() < 0.52 else 0          # 1 = late
        surgery = 1 if rng.random() < 0.60 else 0
        chemo = 1 if rng.random() < 0.43 else 0
        rads = 1 if rng.random() < 0.28 else 0
        pcs = [round(rng.gauss(0, 0.02), 5) for _ in range(3)]

        lp = (0.030 * (age - 65.6) - 0.25 * (sex == 2) + 0.85 * stage
              - 0.88 * surgery + 0.10 * chemo + 0.05 * rads
              + sum(b * x for b, x in zip(beta, g)))
        t_event = rng.expovariate(0.010 * math.exp(lp))
        t_cens = rng.expovariate(1 / 160.0)
        os_month = round(min(t_event, t_cens), 1)
        dead = int(t_event <= t_cens)
        clin.append(dict(IID=f"SIM{k+1:04d}", SEX=sex, AGE=round(age, 1), smoksort=smok,
                         early_late=stage, surgery=surgery, chemotx=chemo, RADS=rads,
                         pc1=pcs[0], pc2=pcs[1], pc3=pcs[2],
                         OS_month=os_month, DEAD=dead))

    with open(f"{a.outdir}/example_genotypes.raw", "w") as f:
        f.write(" ".join(["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"] + snp_names) + "\n")
        for k, g in enumerate(geno):
            row = [f"SIM{k+1:04d}", f"SIM{k+1:04d}", "0", "0", str(clin[k]["SEX"]), "-9"]
            f.write(" ".join(row + [str(x) for x in g]) + "\n")

    cols = ["IID", "SEX", "AGE", "smoksort", "early_late", "surgery", "chemotx",
            "RADS", "pc1", "pc2", "pc3", "OS_month", "DEAD"]
    with open(f"{a.outdir}/example_clinical.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in clin:
            w.writerow(r)

    with open(f"{a.outdir}/example_truth.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["SNP", "effect_allele", "other_allele", "MAF", "true_beta"])
        for j in range(a.snps):
            w.writerow([snp_names[j].rsplit("_", 1)[0], effect_allele[j][0],
                        effect_allele[j][1], mafs[j], round(beta[j], 4)])

    d = sum(c["DEAD"] for c in clin)
    print(f"n = {a.n}, SNPs = {a.snps}, events = {d} ({100*d/a.n:.1f}%), "
          f"non-null variants = {sum(b != 0 for b in beta)}")

def _qnorm(p):
    # Acklam's inverse normal CDF approximation
    if p <= 0: return -8.0
    if p >= 1: return 8.0
    a_ = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
          1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b_ = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
          6.680131188771972e+01, -1.328068155288572e+01]
    c_ = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
          -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d_ = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
          3.754408661907416e+00]
    pl, ph = 0.02425, 1 - 0.02425
    if p < pl:
        q = math.sqrt(-2 * math.log(p))
        return (((((c_[0]*q+c_[1])*q+c_[2])*q+c_[3])*q+c_[4])*q+c_[5]) / ((((d_[0]*q+d_[1])*q+d_[2])*q+d_[3])*q+1)
    if p > ph:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c_[0]*q+c_[1])*q+c_[2])*q+c_[3])*q+c_[4])*q+c_[5]) / ((((d_[0]*q+d_[1])*q+d_[2])*q+d_[3])*q+1)
    q = p - 0.5; r = q * q
    return (((((a_[0]*r+a_[1])*r+a_[2])*r+a_[3])*r+a_[4])*r+a_[5])*q / (((((b_[0]*r+b_[1])*r+b_[2])*r+b_[3])*r+b_[4])*r+1)

if __name__ == "__main__":
    main()
