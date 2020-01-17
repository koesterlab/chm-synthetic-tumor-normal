from cyvcf2 import VCF, Writer

vcf_in = VCF(snakemake.input[0])
vcf_in.add_info_to_header({"ID": "SOMATIC", "Number": 0, "Type": "Flag", "Description": "Somatic variant"})
vcf_in.add_info_to_header({"ID": "VAF", "Number": 1, "Type": "Float", "Description": "Variant allele frequency"})

vcf_out = Writer(snakemake.output[0], vcf_in)

vaf = float(snakemake.wildcards.percentage) / 100

for rec in vcf_in:
    if rec.genotypes == [[1, 0, True]]:
        # if in CHM1 only, consider as somatic
        rec.INFO["SOMATIC"] = True
        rec.INFO["VAF"] = vaf
    vcf_out.write_record(rec)

vcf_out.close()
vcf_in.close()
