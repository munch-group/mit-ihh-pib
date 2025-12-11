#!/bin/bash
#
# Wrapper script to:
# 1. Identify iHS candidates for all populations (with p-values)
# 2. Run pairwise comparison using p-value threshold
#

echo "================================================================================"
echo "Step 1: Identifying iHS candidates for all populations"
echo "================================================================================"
echo ""

python /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/ThisWorks/01a_identify_ihs_candidates_all_populations.py

if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Candidate identification failed!"
    exit 1
fi

echo ""
echo "================================================================================"
echo "Step 2: Running pairwise RELATE vs iHS comparison"
echo "================================================================================"
echo ""

python /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/ThisWorks/07_pairwise_relate_ihs_comparison.py

if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Pairwise comparison failed!"
    exit 1
fi

echo ""
echo "================================================================================"
echo "All steps completed successfully!"
echo "================================================================================"
