# AoI-Pull-Target
The simulation code, energy trace data and the original data for the testbed experiments for the paper "Minimizing the AoI for Pull-Based Target-Level Data Collection in Energy-Harvesting IoTs" 

Files:
1. AoI-Pull-Target, the project for simulational experiments;
2. energy_data, the energy trace used in the simulation with μW as the unit. It is retrieved from an indoor EH embedding wireless testbed by ETH Zürich (https://zenodo.org/records/3363925).
   energy_data/data1119, the energy harvesting rate (With μW as the unit) of each second from 5:00 am to 16:06 am for each day starting from 2018/1/1, e.g., the file "data1119/10" contains the energy data on 2018/1/10.
   energy_data/datasummary1119, the aggregated energy harvesting rate from 5:00 am to 16:06 am of each day starting from 2018/1/1, and each data point is the average value of the energy harvesting rate of one second for consecutive 1000 seconds.
3. experimental_original_data, the original data of testbed experiments.
   experimental_original_data/results.xlsx, the experimental result for testbed experiments.
   experimental_original_data/light_intensity.xlsx, the light intensity at each node in the testbed experiment.

# Setups

1. Import the project "AoI-Pull-Target" in Eclipse.

2. Open the project, set the path for the energy trace:
   (a) aoi_pull_target.java, Line 1410 set "static String energy_dataset_by_day" as the path of folder energy_data/data1119.  
   (b) aoi_pull_target.java, Line 1411 set "energy_dataset_by_day_average" as the path of folder energy_data/datasummary1119.

3. Set simulational parameters as follows:
   
   aoi_pull_target.java, Line 1390-1408:

   static int node_num = 40; //Number of nodes

   static int target_num = 52; //Number of targets
   
   static int wireless_restriction = 10; // The maximum number of nodes that can transmit data during a timeslot
   
   static int working_duration = 1000; // The length of monitoring duration
   
   static int avg_query_interval = 60; //The average query interval
   
   static int energy_buffer = 300; //Energy capacity of each node
   
   static int evaluate_slot = 4; // The energy data to be used in the simulation, range:4-10, 4 refers to the data captured from 5:51am to 6:07am, 5 refers to the data captured from 6:08am to 6:24am, 6 refers to the data captured from 6:25am to 6:42am ... 10 refers to the data captured from 7:33 am to 7:49am.  

4. Get the results: Run "aoi_pull_target.java". All algorithms will repeat for 100 times and then output the simulation results.
   
   
