# Performance Report

## Model: Quark CKM
### 1 step - Full 

To run some check later. 

### 2 steps - Masses + Full 

Masses are ordered increasingly (so that the CKM angles are consistent). It take approximately 6s/iteration, with an average efficiency around 90 %. 	 

**Note:** Save to separate files to prevent overwriting the data accidentally. All data are aggregated in `quark_all.csv` and `quark_cache_all.csv`. 

|#| Date | Machine | Run Time | N | Hit | t/it | Save file |
|:---:| :---: | :---: |  :---: | :---: | :---: | :---: | :---: | :---: | 
|1.1| 2021/04/20 | [M]MacPro | 0:09:46 | 100 | 88 | 5.87s | `quark_01.csv`, `quark_cache_01.csv` |
|1.2| 2021/04/20 | [M]MacPro | 0:09:03 | 100 | 94 | 5.43s | `quark_01.csv`, `quark_cache_01.csv` |
|1.3| 2021/04/20 | [M]MacPro | 0:30:40 | 300 | 275 | 6.14s | `quark_01.csv`, `quark_cache_01.csv` |
|2.1| 2021/04/20 | [M]MacPro | 0:11:27 | 100 | 93 | 6.88s | `quark_02.csv`, `quark_cache_02.csv` |

Run details:  
	
- 1.1, 1.2 and 1.3 search initilizes c-parameters in [0, 1], limit [-2, 5].  
- 2.1 searches initilize c-parameters between [1, 5], limit [-3, 15]. 