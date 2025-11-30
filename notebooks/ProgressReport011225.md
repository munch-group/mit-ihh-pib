# Progress Report 01/12/2025
Author: Mit Van Bruggen 
Date: 29/11/2025 

## Rerunning the iHS analysis 
- I compared the relate results (which turned out to be incomplete) to my ihs results. This showed that the relate results didn't contain the non-PAR region, while my ihs analysis completely skips the PAR regions. Ofcourse this caused a problem 
- I realised that the problem was probably mostly in the relate results, but also wanted to identify why my ihs analysis didn't contain the PAR1 region. 
    - I digged deeper into the documentation of the software I used (hapbin) to see if the format of my input files were wrong, what the limitations are etc. 
    - I found the following: 
        - https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-h 
            - maybe my .hap files for chrX weren't made correctly because i didn't use the --chrX parameter when using impute? (I still don't know this, haven't tried it this way because that parameter does something with phasing, and my 1000G files should already be phased correctly)
            - it also says that impute2 format is where you put dummy - characters for males in non-PAR regions. The october format did something else: it just removed one of the haplotypes for males to compensate for the hemizygosity, producing the variable line count. 
        - In the source code of ihsbin: 
            ```std::getline(file,line);
            m_snpLength = (line.size()+1)/2; 
            ```
            So this assumes the first line determines the haplotype count for all lines. 
            - Hapbin reads the first line of my .hap file (a PAR1 variant with 6404 haplotypes) and sets m_snpLength=6404 for the entire file. So when it reads non-PAR lines with only 4,468 haplotypes, there's a mismatch. This could make it fail on non-PAR variants. 
            - This explains why the results of our october analysis start at position 2,781,309 (right after PAR1): hapbin is probably failing of producing invalid results for all the PAR1 variants because it expects 6404 haplotypes but the format changes at the PAR boundary. 
    - What I did next: 
        - I created a .hap file where all lines have the same number of haplotypes (impute2 format). So now we have hap files with 6404 haplotypes on every line, using - for the missing male haplotypes in non-PAR regions.

# Overview of the 4 iHS analyses 

## Confusion 
I did multiple analyses, but wasn't sure which files were used for it. 
1. October analysis: results/ihs (428,414 snps)
2. November analysis 1: results/ihs_fixed_par1 (433,705 snps)
3. November analysis 2: attempted but failed: input files with dummy - characters
4. November analysis 3: changed the hapbin source code to handle dummy characters and treat them as zero (428414 snps (same as october))

### October analysis (results/ihs)
- input HAP file: chrX.impute.hap.fixed.gz (329M, ASCII text)
- Format: Variable field counts
PAR1 lines: 6404 fields (everyone diploid)
Non-PAR lines: ~4468 fields (males haploid - only 1 copy)
- Biologically: CORRECT, males properly represented as hemizygous
- Technical problem: hapbin reads the first line to determine field count (6404) and then expects all lines to have 6404 fields 
    - when it hits non-par lines with only 4468 fields, it just stops reading early and fills the rest with zeros. 
    - This causes silent data corruption, the analysis completes without crashing but uses incomplete data. 
- Result: produced iHS scores, but they are based on partially incorrect data. 

### November analysis 1 (results/ihs_fixed_par1)
- The name is that way because the objective of the analysis was to make ihsbin work on par1 (but we now know that that is not what it fixed, because ihsbin just does not work on par1 in general)
- Input HAP file: chrX.impute.hap.gz (350M, ASCII text), the ORIGINAL file
- Format: ihsbin takes 6404 fields everywhere BUT there IS a variable field count (but we didn't know that then, I fould out later that it had 4806 in the non-PAR region)
    - Everyone treated as diploid in all regions (this is what we tried to compensate for in october by removing the second haplotype for males, making it 4468)
    - For males in non-PAR: their single haplotype is duplicated (e.g., if male has "0", file shows "0 0")
- Biologically: incorrect because males are treated as diploid when they're actually hemizygous
- Technically: worked but still we had a variable field count (same problem as in october) 
- Problem: iHH scores are artificially inflated because males contribute 2 identical haplotypes instead of 1, which increases the haplotype homozygosity measurements and leads to biased iHS scores. 
- Result: produced 5290 more SNPS than october, but those are biologically wrong 

### November analysis 2: failed (22/11)
- I tried to run reprocess_chrX_hapbin_v2. sh to create files with dummy - for males. 
- Output created: chrX.impute.hap.hapbin_format.gz (1.3 GB)
    - has 6404 fields everywhere
    - males in non-PAR: real haplotype + - (dummy character) 
    - biologically correct 
- it failed because the script stopped after creating the hapbin_format.gz file 
    - it never ran convert_to_hapbin.py to create the final hapbin format 
    - reason: i discovered in hapbin's source code that it crashes on - characters. The convert() function throws and exception for any character except '0', '1', space or newline. 
**IHS scores really differ**
Example variant X:2781309:
- October: iHH_0=9.77, iHH_1=7.17 → Lower values (but based on corrupted data)
- November (fixed_par1): iHH_0=10.58, iHH_1=10.53 → Higher values (males double-counted)
The November scores are artificially higher because males' duplicated haplotypes increase homozygosity measurements.

The November (fixed_par1) analysis treats males incorrectly:
- Doubles their contribution in non-PAR regions
- Makes extended haplotype homozygosity appear higher than it really is
- Could lead to false positive signals of selection
- The new approach should give the most accurate iHS scores for chromosome X.



### November analysis 3: worked
- New Solution
    1. Modified hapbin source code (hapbin_source/src/hapbin.hpp)
        - Added support for - character (lines 177-181)
        - Treats - as missing data (increments position but doesn't set bit)
        - Built custom ihsbin binary
    2. New analysis will use:
        - Input: chrX.impute.hap.hapbin_format.gz (already created, has dummy -)
        - Tool: Modified ihsbin that handles - characters
        - Output: New directory hapbin_with_dummy/ and results/ihs_with_dummy/
        - This approach is both biologically correct AND technically compatible
    
    When creating chrX.impute.hap.hapbin_format.gz (so the file with equal field counts because of the dummies ) with script: 
    /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/fix_chrX_hapbin_v2_13385889.err

    /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/ThisWorks/reprocess_chrX_hapbin_v2.sh

    /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/ThisWorks/fix_chrX_haplotypes_for_hapbin.py 

    Which use this as input: 
    ```
    VCF="$DATA_ROOT/raw/chr${CHR}/chr${CHR}.vcf.gz"
    HAP_IN="$WORK_DIR/chr${CHR}.impute.hap.gz"
    LEG_IN="$WORK_DIR/chr${CHR}.impute.legend.gz"
    ``` 

    I got many warnings saying that there were 4806 haps found and 6404 expected. But at this point I was under the impression that our chrX.vcf.gz file already had equal field counts (this is what I based my fixed_par1 analysis on and what was the biological mistake: we were treating males as diploid). 

 