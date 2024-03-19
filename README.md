# AoI-Pull-Target
The simulation code, energy trace data and original data for the testbed experiments for the paper "Minimizing the AoI for Pull-Based Target-Level Data Collection in Energy-Harvesting IoTs" 

Files:
1. AoI-Pull-Target, the project for simulational experiments;
2. energy_data, the energy trace used in the simulation with μW as the unit.
   energy_data/data1119, the raw energy data for each day starting from 2018/1/1.
   energy_data/datasummary1119, the aggregated energy data of each day, and each data point is the average retrieved from an indoor EH embedding wireless testbed by ETH Zürich (https://zenodo.org/records/3363925).
3. experimental_original_data, data related to testbed experiments.
   experimental_original_data/results.xlsx, the experimental result for testbed experiments.
   experimental_original_data/light_intensity.xlsx, the light intensity at each node in the testbed experiment.

# Simulation Setups

1. Import the project "AoI-Pull-Target" in Eclipse.

2. Open the project, set the path for the energy trace:
   (a) aoi_pull_target.java, Line 1410 set "static String energy_dataset_by_day" as the path of folder energy_data/data1119.  
   (b) aoi_pull_target.java, Line 1411 set "energy_dataset_by_day_average" as the path of folder energy_data/datasummary1119.

3. Set simulational parameters as follows:
   
   aoi_pull_target.java, Line 1390-1406:

   static int node_num = 40; //Number of nodes
   static int target_num = 52; //number of targets of interest
   static int pkt_header = 12; //pkt header length
   static int coverage_range = 24; //coverage range of each node
   static int wireless_capacity = 100; // bandwidth constraint
   static double average_cover_num = 7.8; // The average number of nodes that can cover each target
   static int working_duration = 1000; // length of monitoring duration
   static int avg_query_interval = 60; // average query interval
   static int avg_query_length = 4; //average query length
   static int energy_buffer = 300; //energy capacity
   static int evaluate_slot = 4; // The energy data to be used, range:4-10, 4 refers to the data captured from 5:51am to 6:07am, 5 refers to the data captured from 6:08am to 6:24am, 6 refers to the data captured from 6:25am to 6:42am ... 10 refers to the data captured from 7:33 am to 7:49am.  


4. Get the results: Run "aoi_pull_target.java". All algorithms will repeat for 100 times and then output the simulation results.
   
   
