#!/bin/bash

raster_version="v1.1.0"
vector_version="v1.1.1"

# 128gb memory
# for study_area in 2083 660 2423 3436 3437 3641 3539 3538 2374 2169 3089 2320 3640 2373 2372 2271 2270 2784 1969 2377 2379 2830 1813 676 780 678 677 2891 2890 2077 2789 2380 2278 2790 2889 1192 880 670 669 2282 984 567 559 773 774 3193 1098 1090 988 1088 1091 563 461 462 1982 1981 1880 1879 1779 1778 882 1092 1097 460 3297 3298 1089 987 565 568 778 881 985 776 675 3502 3400 661 1979 3500 1472 1575 2385 3501 3294 560 459 458 3399 3296 1301 1200 1199 2381 1405 1303 862 761 760 777 569 463 2082 1980 674 468 570 571 3091 3397 2183 2080 2582 3192 2179 2079 2078 3233 2788 4150 2626 2482 3131 3029 3132 2990 2785 467 566 464 2280 466 465 2688 1404 1302 2013 2218 2217 2115 2116 1711 2787 2687 2827 2886 2786 3031 3030 2929 2928 2319 2829 2828 2726 2523 2015 1913 1912 2421 1609 1608 1811 2625 1507 1506 4819 2588 2587 2383 2384 1298 1196 564 2181 1299 1197 1710 1297 1198 1096 764 662 865 2487 664 663 2382 562 561 1677 1676 1474 1371 1295 1296 1093 986 964 966 965 1066 1271 1270 1169 1168 5237 4107 3798 3800 4209 4208 3805 4307 4104 3799 4413 4311 4514 4003 4002 3901 3900 4207 5337 5432 3807 3705 3603 3602 5233 4308 4616 3905 4309 5134 5235 3804 3702 3701 3803 3802 5543 4511 4409 4512 5441 5339 5023 4921 5440 5437 4614 5335 5433 5436 5435 5334 5333 5331 5332 5127 4923 5330 4717 5328 5226 4310 4615 3806 5131 5029 5329 5024 5031 5133 5132 5230 5128 5232 4820 4718 4006 3904 3902 4922 4005 5231 5129 8627 8626 8625 8533 8531 8530 8529 8527 8434 8433 8432 8431 8430 8429 8428 8427 8426 8423 8422 8421 8419 8418 8335 8334 8333 8332 8331 8330 8329 8234 8233 8232 8231 8228 8227 8216 8215 8214 8211 8210 8133 8132 8131 8113 8111 8108 8107 8032 8030 8008 8007 8006 7929 7928 7906 7905 7904 7900 7832 7831 7828 7800 7799 7797 7731 7729 7727 7696 7632 7631 7630 7629 7628 7593 7592 7531 7530 7529 7490 7489 7385 7384 7285 7284 7283 7282 7276 7275 7178 7177 7176 7175 7174 7172 7061 6969 6968 6965 6958 6957 6867 6866 6863 6862 6861 6860 6857 6856 6761 6760 6759 6758 6757 6756 6754 6753 6657 6656 6557 6554 6553 6551 6550 6461 6460 6459 6458 6457 6361 6360 6357 6356 6355 6354 6260 6258 6160 6057 5956 5849 5749 5647 5646 5545 6728 6414 6923 4769 4768 4458 4560 5077 5076 5481 5792 5690 5793 5691 5478 5479 5480 4048 3946 3945 5689 5587 5586 5584 5482 5585 5483 6724 7029 6927 6932 5485 5384 5180 5179 5379 5282 5375 5381 5278 4047 4150 4355 3844 3843 4970 5693 4251 4149 3742 3741 4871 4457 4456 4354 4252 4763 5795 5376 7030 5378 5275 7135 6929 7132 5796 4971 4869 5174 5073 4766 4664 4767 4665 6828 4663 4870 7124 5794 5692 5176 4867 4765 4764 5074 4973 4972 4660 4558 7123 6309 7023 7022 7225 6717 7429 6714 7033 6925 6411 7327 6103 6612 7134 7133 7031 6818 6512 6410 6920 6921 6614 6819 6001 6000 5899 5898 6308 6307 6206 6205 6725 6620 6519 6514 6412 7229 7230 6823 6722 6518 6621 6830 6515 7026 7028 7027 6931 6721 6622 7231 6625 6413 6727 7233 7131 6623 6726 6624 6829 6517 6516 6618 6619 1473 763 2988 2481 2081 3398 8528 8532 8031 7281 4661 5279 5280 6613 6716 5277 6453 6555 6549 6755 6864 7798 8213 7930 7931 7830 7730 7024 7025 7127 7130 6930 7932 6824 6926 6715 5281 5383 5484 5380 5177 5377 5276 5175 4868 3335 2727 2728 2422 2524 2014 1812 1193 673 1372 1574 2182 2180 2486 2888 4004 4412 4513 5130 6820 7032 8112

# 256gb memory
# for study_area in 5852 2584 5030 5228 4106 3903 3704 2989 5338 4411 4105 1067 2484 1195 2485 2586 2685 2686 2281 762 5126 5125 5227 1194 7901 8109 7129 671 672 864 8212 7903 5851 5544 5236 6651 6358 7902 7071 3090 2683 863 7232 779 2483 2583 2585 2684 2887 3295 3703 4206 4410 4559 4662 5075 5178 5229 5382 5850 5853 5952 5953 5954 5955 6158 6159 6259 6261 6353 6359 6454 6455 6456 6511 6552 6556 6652 6653 6654 6858 6865 6922 7070 7128 7173 7279 7280 7386 7387 7388 7695 8004 8005 8110 8316 8317 8318 8319 8420 8524 8525 8526

# Validation tiles
# for study_area in 464 560 563 567 569 570 763 777 778 880 881 882 984 985 986 987 1089 1090 1091 1096 1097 1168 1169 1193 1194 1195 1196 1197 1198 1199 1200 1302 1404 1405 1473 1474 1507 1608 1609 1710 1711 1811 1812 1913 2014 2015 2081 2115 2116 2181 2182 2183 2217 2280 2319 2384 2421 2485 2486 2523 2625 2726 2727 2828 6205 6410 6512 6613 6619 6622 6716 6721 6727 6818 6824 6829 6926 6930 6931 7028 7031 7032 7130 7132 7629 7630 7631 7729 7730 7830 7831 7930 7931 8031 8032 8133 8233 8234 8332 8333 8334 8432 8433 8434 8527 8528 8530 8531 8532 8533 4661 4765 4868 4972 5073 5074 5175 5179 5276 5278 5279 5280 5281 5377 5378 5379 5380 5381 5383 5484 5483 5485 5586 5587 5690 5691 5692 5693 5793 5795 5898 6655 8003 5025

# Experiments
for study_area in 2374
# for study_area in 880

do

    PBS="#!/bin/bash\n\
    #PBS -N ${study_area}_${vector_version}\n\
    #PBS -o PBS/DEACoastlines_${study_area}_${vector_version}.out\n\
    #PBS -e PBS/DEACoastlines_${study_area}_${vector_version}.err\n\
    #PBS -l storage=gdata/v10+gdata/r78+gdata/xu18+gdata/fk4\n\
    #PBS -P r78\n\
    #PBS -q hugemem\n\
    #PBS -l walltime=12:00:00\n\
    #PBS -l mem=256GB\n\
    #PBS -l jobfs=2GB\n\
    #PBS -l ncpus=1\n\
    #PBS -l wd\n\
    module use /g/data/v10/public/modules/modulefiles\n\
    module load dea/20200713\n\
    module load otps\n\
    python3 /g/data/r78/DEACoastlines/deacoastlines_statistics.py $study_area $raster_version $vector_version"
#         python3 /g/data/r78/DEACoastlines/deacoastlines_generation.py $study_area $raster_version $vector_version"

    echo -e ${PBS} | qsub || echo "${study_area} failed" >> log.txt
    sleep 0.2
    echo "Submitting study area $study_area $vector_version"

done
