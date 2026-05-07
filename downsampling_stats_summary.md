# Downsampling Stats Summary

## Overview

- **Dataset**: 1,979 cells, 27,158 total genes
- **Original sequencing depth**: ~1.29M counts/cell (very deep SMART-seq)
- **Original gene detection**: ~6,440 genes/cell, 27,090 total genes detected

## Results by Target Depth

| Target Depth | Avg Genes/Cell (Original → Downsampled) | Total Genes Detected (Original → Downsampled) |
|-------------|----------------------------------------|----------------------------------------------|
| 50 | 6,440 → 46 | 27,090 → 12,400 |
| 5,000 | 6,440 → 1,795 | 27,090 → 24,200 |
| 10,000 | 6,440 → 2,455 | 27,090 → 25,300 |
| 15,000 | 6,440 → 2,843 | 27,090 → 25,800 |
| 25,000 | 6,440 → 3,300 | 27,090 → 26,270 |
| 50,000 | 6,440 → 3,840 | 27,090 → 26,675 |
| 100,000 | 6,440 → 4,295 | 27,090 → 26,855 |

## Key Observations

1. **Extreme downsampling (50 counts/cell)**: Only ~46 genes/cell detected, losing ~55% of total genes
2. **5,000 counts/cell**: ~1,795 genes/cell (28% of original), retains 89% of total genes
3. **100,000 counts/cell**: ~4,295 genes/cell (67% of original), retains 99% of total genes
4. **Replicate consistency**: Seeds 42 and 123 show very consistent results

## Conclusion

The downsampling is working correctly - gene detection drops as expected with lower sequencing depth.
