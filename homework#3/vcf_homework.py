import re

vcf_file = 'homework#3/307_genes_for_CP.search_strict_match.vcf'

depth_over_30 = 0
depth_per_sample_over_30 = 0
mnv = 0
homo = 0
hetero = 0
found_rs9697 = False

smplsnmbr = 0

with open(vcf_file) as f:
    for line in f:
        if line.startswith('##'):
            continue
            
        if line.startswith('#'):
            cols = line.strip().split('\t')
            info_col = cols.index('INFO')
            format_col = cols.index('FORMAT')
            ref_col = cols.index('REF')
            alt_col = cols.index('ALT')
            id_col = cols.index('ID')
            smplsnmbr = len(cols) - format_col - 1
            continue
        
        parts = line.strip().split('\t')
        
        

        # question2
        info = parts[info_col]
        dp = re.search(r'DP=(\d+)', info)
        if dp:
            total_dp = int(dp.group(1))
            if total_dp > 30:
                depth_over_30 += 1
            
            avg_dp = total_dp / smplsnmbr   #counting this out of personal curiousity
            if avg_dp > 30:
                depth_per_sample_over_30 += 1

        # question3
        ref = parts[ref_col]
        alt = parts[alt_col]
        for alt_allele in alt.split(','):
            if len(ref) > 1 and len(alt_allele) > 1 and len(ref) == len(alt_allele):
                mnv += 1
        
        # question4
        has_homo = False
        has_hetero = False
        
        for i in range(format_col + 1, len(parts)):
            gt = parts[i].split(':')[0]
            
            if '/' in gt:
                alleles = gt.split('/')
            elif '|' in gt:
                alleles = gt.split('|')
            else:
                continue
            
            if len(alleles) == 2:
                if alleles[0] == alleles[1] and alleles[0] not in ['.', '0']:
                    has_homo = True
                elif alleles[0] != alleles[1] and '.' not in alleles and '0' in alleles:
                    has_hetero = True
        
        if has_homo:
            homo += 1
        if has_hetero:
            hetero += 1
        
        # question5
        var_id = parts[id_col]
        if 'rs9697' in var_id:
            found_rs9697 = True
            print(f"\nFound rs9697 at position {parts[1]}")
            print(f"Full line: {line[:200]}...")
            
print(f"1. - Number of samples: {smplsnmbr}") 
print(f"2a. - Variants with total depth >30: {depth_over_30}")
print(f"2b. - number of depth >30 per sample: {depth_per_sample_over_30}")
print(f"3. - MNV variants: {mnv}")
print(f"4. - Homozygous variants: {homo}")
print(f"5. - Heterozygous variants: {hetero}")
print(f"6. - rs9697 found: {found_rs9697}")