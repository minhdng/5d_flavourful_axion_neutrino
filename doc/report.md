# Performance Report

## Model: Quark CKM
### 1 step - Full 

|#| Date | Machine | Run Time | N | Hit | Save file |
|:---:| :---: | :---: |  :---: | :---: | :---: | :---: | :---: | 
|1| 2021/04/17 | Minh's Mac Pro | 0:02:24 | 200 | 0 | NA |
|2| 2021/04/17 | Minh's Mac Pro | 0:20:25 | 2000 | 0 | NA |

### 2 steps - Masses + Full 
Masses are ordered increasingly (so that the CKM angles are consistent). It take approximately 11s / iteration, with an average efficiency of 3 - 4 %. 
 
|#| Date | Machine | Run Time | N | Hit | Save file |
|:---:| :---: | :---: |  :---: | :---: | :---: | :---: | :---: | 
|3| 2021/04/17 | Minh's Mac Pro | 0:21:27 | 100 | 7 | quark_04.csv |
|4| 2021/04/17 | Minh's Mac Pro | 4:32:20 | 1500 | 52 | quark_04.csv |
|5| 2021/04/18 | Minh's Mac Pro | 0:23:39 | 100 | 4 | quark_05.csv |
|5| 2021/04/18 | Minh's Mac Pro | 7:41:37 | 2500 | 99 | quark_05.csv |
|6| 2021/04/18 | Minh's Mac Pro | 3:05:20 | 1000 | 41 | quark\_06.csv, quark\_cache\_06.csv  |
|7| 2021/04/19 | Minh's Mac Pro | 1:11:44 | 498 | 41 | quark\_07.csv, quark\_cache\_07.csv |

Import updates:  
	
	#7: Fix a bug in `full` calculation
	 

**Note:** Save to separate file to prevent overwriting the data accidentally.  