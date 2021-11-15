#!/bin/bash

for output_name in v1.0.0

do

    PBS="#!/bin/bash\n\
    #PBS -N ${output_name}\n\
    #PBS -o scripts/PBS/DEACoastlines_summary_${output_name}.out\n\
    #PBS -e scripts/PBS/DEACoastlines_summary_${output_name}.err\n\
    #PBS -l storage=gdata/v10+gdata/r78+gdata/xu18+gdata/fk4\n\
    #PBS -P r78\n\
    #PBS -q hugemem\n\
    #PBS -l walltime=06:00:00\n\
    #PBS -l mem=256GB\n\
    #PBS -l jobfs=2GB\n\
    #PBS -l ncpus=1\n\
    #PBS -l wd\n\
    module use /g/data/v10/public/modules/modulefiles\n\
    module load dea\n\
    module load otps\n\
    python3 /g/data/r78/DEACoastlines/deacoastlines_summary.py $output_name"

    echo -e ${PBS} | qsub || echo "${output_name} failed" >> log.txt
    sleep 0.2
    echo "Submitting summary $output_name"

done
