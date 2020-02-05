import json
import requests
import csv

# This script reads a vcf file to prioritize the effects of variants
# using VEP variant impact. Variants causing multiple effects are prioritize
# taking into account VEP ranking of high, moderate, low and modifier categories.
# Usage:
# python VCF_parsing.py
# Output file: output_vcf_03_02.txt
# The outout produced contains the following columns:
# chromosome, position, reference allele, alternate allele, depth of sequence coverage, 
# number of reads for alternate, alternate allele fraction, allele frequency in exac,
# variant consequence, impact, biotype of consequence

class_impact = {} #

# Reading impact table
with open("vep_impact.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        if "ablation" not in row:
            class_impact[row[0]] = row[1]

# Level of impact
dict_impact = {"HIGH": 4, "MODERATE":3, "LOW":2, "MODIFIER":1, "NA":0}

with open("Challenge_data.vcf", "r") as f1, open("output_vcf_03_02.txt", "w") as f2:
    contador = 1
    for line in f1:
        if not line.startswith("#"):
            line = line.rstrip().split()
            chr = line[0]
            pos = line[1]
            ref = line[3]
            alt = line[4].split(",")
            extra = line[7].split(";")
            DPs = extra[7].split("DP=", 1)[1].split(",")
            AOs = extra[5].split("AO=", 1)[1].split(",")
            DP = int(DPs[0])
            for sub_mut in range(len(alt)):
                AO = int(AOs[sub_mut])
                string = "-".join([chr,pos,ref,alt[sub_mut]])
                r = requests.post("http://exac.hms.harvard.edu/rest/bulk/variant", json=[string])
                contador += 1
                if r.status_code == 200:
                    parsed_json = json.loads(r.text)
                    for feature in parsed_json[string]:
                        if feature == "variant":
                            if "allele_freq" in parsed_json[string][feature]:
                                freq = [str(parsed_json[string][feature]["allele_freq"])]
                            else:
                                freq = ["NA"]
                            if "vep_annotations" in parsed_json[string][feature]:
                                vep_annot = parsed_json[string][feature]["vep_annotations"]
                                total_conseqs = [annot["major_consequence"] for annot in vep_annot]
                                total_impact_conseq = []
                                try:
                                    total_impact_conseq = [class_impact[conseq] for conseq in total_conseqs]
                                except:
                                    total_impact_conseq = ["NA"] * len(total_conseqs)
                                total_impact_value = [dict_impact[conseq] for conseq in total_impact_conseq]
                                biotype_dict = {annot["major_consequence"]:annot["BIOTYPE"] for annot in vep_annot}
                                indices = [i for i,x in enumerate(total_impact_value) if x == max(total_impact_value)]
                                sub_conseq = [total_conseqs[i] for i in indices]
                                innecesary_dict = {}
                                for i in sub_conseq:
                                    if i not in innecesary_dict:
                                        innecesary_dict[i] = 1
                                    else:
                                        innecesary_dict[i]+=1
                                top_conseq = [k for k,v in innecesary_dict.items() if v == max(innecesary_dict.values())]
                                try:
                                    impact = [class_impact[i] for i in top_conseq]
                                except:
                                    impact = ["NA"] * len(top_conseq)
                                try:
                                    biotype = [biotype_dict[i] for i in top_conseq]
                                except:
                                    biotype = ["NA"] * len(top_conseq)
                            else:
                                freq = ["NA"]
                                top_conseq = ["NA"]
                                impact = ["NA"]
                                biotype = ["NA"]
                else:
                    freq = ["NA"]
                    impact = ["NA"]
                    biotype = ["NA"]

                freq = ",".join(freq)
                if not top_conseq:
                   top_conseq = ["NA"]
                top_conseq = ",".join(top_conseq)
                impact = ",".join(impact)
                biotype = ",".join(biotype)
                f2.write(chr + "\t" + pos + "\t" + ref + "\t" + alt[sub_mut] + "\t" + str(DP) + "\t" + str(AO) + "\t" + str(float(AO/DP)) + "\t" + freq + "\t" + top_conseq + "\t" + impact + "\t" + biotype + "\n")