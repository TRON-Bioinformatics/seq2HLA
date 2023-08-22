#!/usr/bin/env python

import numpy as np
import scipy.stats as st


def calculate_confidence(allele_count, allele_vector):
    paril = 1.0 - st.norm.cdf(allele_count, np.mean(allele_vector), np.std(allele_vector, ddof = 1))
    poutlier = st.binom.pmf(0, len(allele_vector), paril)
    return (1.0 - poutlier)


def main():
    val = [8498, 6256, 25164]
    types = ["A*24", "B*35", "C*04"]
    val_list = [
        [
            1244,
            5940,
            1722,
            660,
            5770,
            222,
            688,
            374,
            544,
            1890,
            636,
            4272,
            192,
            132,
            212,
            212,
            1150,
            2426,
            596,
            160
        ],
        [
            2718,
            782,
            1754,
            1284,
            5832,
            2436,
            1088,
            654,
            630,
            924,
            2086,
            420,
            344,
            1316,
            1314,
            4534,
            310,
            3386,
            1664,
            896,
            2002,
            1372,
            4500,
            3398,
            1374,
            3508,
            2194,
            3488,
            1192,
            300,
            436,
            2846,
            78,
            276,
            322
        ],
        [
            6002,
            2068,
            2332,
            6386,
            3518,
            3078,
            5418,
            3028,
            4836,
            2232,
            2034,
            2096,
            16068
        ]
    ]

    print("['[1] A*24', '[1] 8498', '[1] 1452.1', '[1] 1810.365', '[1] 0.0009938365', '[1] B*35', '[1] 6256', '[1] 1761.657', '[1] 1430.283', '[1] 0.02892281', '[1] C*04', '[1] 25164', '[1] 4545.846', '[1] 3800.572', '[1] 3.766758e-07', '']")

    
    print(calculate_confidence(val[0], val_list[0]))
              

if __name__ == "__main__":
    main()
