import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Random;
import java.util.Set;

class node {
	public int cover_target_num;
	public int[] eh_rate; // The predicted energy harvest rate
	public int[] actual_eh_rate; // The actual energy harvest rate
	public int[] cover_target;
	public int data_rate;
	public double avg_ehrate;
	public boolean[] active_slots;
	public int energy_capacity;

	public node(int duration_length, int energy_capacity) {
		this.cover_target_num = 0;
		this.eh_rate = new int[duration_length + 1];
		this.cover_target = new int[100];
		this.data_rate = 0;
		this.avg_ehrate = 0;
		this.active_slots = new boolean[duration_length + 1];
		for (int i = 1; i <= duration_length; i++) {
			this.active_slots[i] = false;
		}
		this.energy_capacity = energy_capacity;
		this.actual_eh_rate = new int[duration_length + 1];
	}
}

class target {
	public double position_x;
	public double position_y;
	public int cover_node_num;
	public int[] query_number;
	public int[] cover_node;
	public int weight;
	public boolean[] is_covered;
	public int[] current_aoi_value;

	public target(int duration_length, boolean has_weight) {
		this.position_x = Math.random() * 100;
		this.position_y = Math.random() * 100;
		this.cover_node_num = 0;
		this.query_number = new int[duration_length + 1];
		Random xc = new Random();
		double p = 0.1 + 1.8 * xc.nextDouble();
		this.weight = (int) Math.floor(p * 10);
		if (this.weight == 0) {
			this.weight = 1;
		}
		if (has_weight == false) {
			this.weight = 1;
		}
		for (int i = 1; i <= duration_length; i++) {
			this.query_number[i] = 0;
		}
		this.current_aoi_value = new int[duration_length + 1];
		this.cover_node = new int[100];
		this.is_covered = new boolean[duration_length + 1];
		for (int i = 1; i <= duration_length; i++) {
			this.is_covered[i] = false;
		}
	}
}

class test_compare {
	public static Map<Integer, Integer> valueUpSort(Map<Integer, Integer> map) {

		if (map == null || map.isEmpty()) {
			return null;
		}

		Map<Integer, Integer> sortedMap = new LinkedHashMap<Integer, Integer>();

		ArrayList<Map.Entry<Integer, Integer>> entryList = new ArrayList<Map.Entry<Integer, Integer>>(map.entrySet());

		Collections.sort(entryList, new Comparator<Map.Entry<Integer, Integer>>() {
			@Override
			public int compare(Map.Entry<Integer, Integer> o1, Map.Entry<Integer, Integer> o2) {
				return o2.getValue().compareTo(o1.getValue());
			}
		});

		for (Map.Entry<Integer, Integer> r : entryList) {
			sortedMap.put(r.getKey(), r.getValue());
		}

		return sortedMap;
	}

}

class MDP_policy_iteration {
	public int maximum_aoi;
	public int maximum_energy;
	public double policy_value[][][];
	public boolean policy_action[][][];

	public MDP_policy_iteration(int maximum_aoi, int maximum_energy) {
		this.maximum_aoi = maximum_aoi;
		this.maximum_energy = maximum_energy;
		this.policy_value = new double[maximum_aoi + 1][maximum_energy + 1][2];
		this.policy_action = new boolean[maximum_aoi + 1][maximum_energy + 1][2];
		for (int i = 1; i <= maximum_aoi; i++) {
			for (int j = 0; j <= maximum_energy; j++) {
				for (int k = 0; k <= 1; k++) {
					this.policy_value[i][j][k] = -1000;
					double r = Math.random();
					if (r <= 0.5) {
						this.policy_action[i][j][k] = true;
					} else {
						this.policy_action[i][j][k] = false;
					}
				}
			}
		}
	}
}

public class aoi_pull_target {

	static double average_table(int[] table, int table_length) {
		double sum = 0;
		for (int i = 0; i < table_length; i++) {
			sum += table[i];
		}
		return (sum / (double) table_length);
	}

	static int[] generate_average_number(int start, int end, int average, int number) {
		int result[] = new int[number];
		for (int i = 0; i < number; i++) {
			Random r = new Random();
			result[i] = r.nextInt(end - start + 1) + start;
		}
		while (true) {
			double mean = average_table(result, number);
			double p = mean - ((double) average);
			if (Math.abs(p) <= 0.01) {
				break;
			} else {
				if (mean < ((double) average)) {
					int q = 0;
					for (int j = 0; j < number; j++) {
						if (result[j] < average) {
							q = j;
							break;
						}
					}
					Random r = new Random();
					result[q] = r.nextInt(end - average + 1) + average;
				} else {
					int q = 0;
					for (int j = 0; j < number; j++) {
						if (result[j] > average) {
							q = j;
							break;
						}
					}
					Random r = new Random();
					result[q] = r.nextInt(average - start + 1) + start;
				}
			}
		}

		ArrayList<Integer> list = new ArrayList<Integer>();
		for (int i = 0; i < number; i++) {
			list.add(result[i]);
		}
		Collections.shuffle(list);
		for (int i = 0; i < number; i++) {
			result[i] = list.get(i);
		}
		return result;
	}

	static target[] generate_all_target_model(int target_num, int duration_length, boolean has_weight) {
		target[] all_target_model = new target[target_num];
		for (int i = 0; i < target_num; i++) {
			all_target_model[i] = new target(duration_length, has_weight);
		}
		return all_target_model;
	}

	static double compute_distance(double x, double y, double x1, double y1) {
		return Math.sqrt(Math.pow(x - x1, 2) + Math.pow(y - y1, 2));
	}

	static void deploy_node_target(int node_num, node[] all_node_model, int target_num, target[] all_target_model,
			int coverage_range, int pkt_header, double avg_cover_num) {
		while (true) {
			for (int i = 0; i < node_num; i++) {
				all_node_model[i].cover_target_num = 0;
			}
			for (int i = 0; i < target_num; i++) {
				all_target_model[i].cover_node_num = 0;
			}

			boolean has_empty_node = false;
			double s = 0;
			boolean[] is_covered = new boolean[target_num];
			for (int i = 0; i < target_num; i++) {
				is_covered[i] = false;
			}
			for (int i = 0; i < node_num; i++) {
				node this_node = all_node_model[i];
				double x = Math.random() * 100;
				double y = Math.random() * 100;
				for (int j = 0; j < target_num; j++) {
					target this_target = all_target_model[j];
					double d = compute_distance(x, y, this_target.position_x, this_target.position_y);
					if (d <= coverage_range) {
						is_covered[j] = true;
						this_node.cover_target_num++;
						this_node.cover_target[this_node.cover_target_num - 1] = j + 1;
						this_target.cover_node_num++;
						this_target.cover_node[this_target.cover_node_num - 1] = i + 1;
					}
				}
				s += (double) (this_node.cover_target_num);
				if (this_node.cover_target_num == 0) {
					has_empty_node = true;
				}
			}

			boolean all_covered = true;

			for (int i = 0; i < target_num; i++) {
				if (is_covered[i] == false) {
					all_covered = false;
					break;
				}
			}
			int maximum_cover_target = 0;
			for (int i = 0; i < node_num; i++) {
				if (all_node_model[i].cover_target_num > maximum_cover_target) {
					maximum_cover_target = all_node_model[i].cover_target_num;
				}
			}
			if (all_covered == true && maximum_cover_target <= 20) {
				if ((avg_cover_num) * node_num <= s && ((double) avg_cover_num + 1) * node_num >= s) {
					if (has_empty_node == false) {
						break;
					}
				}
			}
		}
		for (int i = 0; i < node_num; i++) {
			all_node_model[i].data_rate = 1;
		}

	}

	static void generate_query(int node_num, node[] all_node_model, int target_num, target[] all_target_model,
			int avg_query_interval, int avg_query_length, int duration_length) {

		for (int i = 0; i < target_num; i++) {
			target this_target = all_target_model[i];
			ArrayList<Integer> query_time = new ArrayList<Integer>();
			ArrayList<Integer> query_last_time = new ArrayList<Integer>();
			int q = 0;
			while (true) {
				Random r = new Random();
				int u = r.nextInt(2 * avg_query_interval - 1) + 1;
				q += u;
				if (q <= duration_length) {
					query_time.add(q);
				} else {
					break;
				}
			}
			for (int k = 0; k < query_time.size(); k++) {
				Random r = new Random();
				int u = r.nextInt(2 * avg_query_length - 1) + 1;
				query_last_time.add(u);
			}
			for (int k = 0; k < query_time.size(); k++) {
				int s = query_time.get(k);
				int t = query_last_time.get(k);
				for (int j = 1; j <= t; j++) {
					if (s + j - 1 <= duration_length) {
						this_target.query_number[s + j - 1]++;
					}
				}
			}
			for (int k = 1; k <= duration_length; k++) {
				this_target.query_number[k] = this_target.query_number[k] * this_target.weight;
			}
		}
	}

	static Set<Integer> intersect(int[] a, int size, Set<Integer> b) {
		Set<Integer> result = new HashSet<Integer>();
		for (int i = 0; i < size; i++) {
			int g = a[i];
			if (b.contains(g) == true) {
				result.add(g);
			}
		}
		return result;
	}

	static boolean is_remain_query(node this_node, int timeslot, target[] all_target_model) {
		boolean result = false;
		for (int i = 0; i < this_node.cover_target_num; i++) {
			target this_target = all_target_model[this_node.cover_target[i] - 1];
			if (this_target.is_covered[timeslot] == false && this_target.query_number[timeslot] > 0) {
				result = true;
				break;
			}
		}
		return result;
	}

	static Map<Integer, Integer> push_back(node[] all_node_model, target[] all_target_model,
			Map<Integer, Integer> active_slot_re, int[] current_capacity, int duration_length, int max_capacity) {

		Map<Integer, Integer> can_slot = new HashMap<Integer, Integer>();
		Map<Integer, Integer> return_slot = new HashMap<Integer, Integer>();
		int base_time = 9999;
		for (Integer n_num : active_slot_re.keySet()) {
			if (active_slot_re.get(n_num) < base_time) {
				base_time = active_slot_re.get(n_num);
			}
		}

		for (Integer n_num : active_slot_re.keySet()) {

			int start_time = active_slot_re.get(n_num);
			node this_node = all_node_model[n_num - 1];
			int push_back_time = -1;

			for (int j = start_time; j <= duration_length; j++) {
				if (is_remain_query(this_node, j, all_target_model) == true) {
					push_back_time = j;
					if (push_back_time > base_time + 99) {
						push_back_time = start_time;
					}
					break;
				}
			}
			if (push_back_time != -1 && push_back_time != start_time) {
				can_slot.put(n_num, push_back_time);
			} else if (push_back_time != -1 && push_back_time == start_time) {
				return_slot.put(n_num, start_time);
				this_node.active_slots[start_time] = true;
				for (int i = 0; i < this_node.cover_target_num; i++) {
					int u = this_node.cover_target[i];
					target this_target = all_target_model[u - 1];
					this_target.is_covered[start_time] = true;
				}
			} else {
				current_capacity[start_time] -= this_node.data_rate;
			}
		}
		if (can_slot.keySet().size() != 0) {
			Set<Integer> push_set = new HashSet<Integer>();
			push_set.addAll(can_slot.keySet());
			while (true) {
				Map<Integer, Integer> vary_set = new HashMap<Integer, Integer>();
				for (Integer p : push_set) {
					node this_node = all_node_model[p - 1];
					if (current_capacity[can_slot.get(p)] + this_node.data_rate <= max_capacity) {
						vary_set.put(p, can_slot.get(p));
						current_capacity[can_slot.get(p)] += this_node.data_rate;
						current_capacity[active_slot_re.get(p)] -= this_node.data_rate;
					}
				}
				if (vary_set.keySet().size() == 0) {
					for (Integer p : push_set) {
						node this_node = all_node_model[p - 1];
						int t = active_slot_re.get(p);
						return_slot.put(p, t);
						this_node.active_slots[t] = true;
						for (int i = 0; i < this_node.cover_target_num; i++) {
							int u = this_node.cover_target[i];
							target this_target = all_target_model[u - 1];
							this_target.is_covered[t] = true;
						}
					}
					break;
				} else {
					for (Integer p : vary_set.keySet()) {
						node this_node = all_node_model[p - 1];
						int t = vary_set.get(p);
						return_slot.put(p, t);
						this_node.active_slots[t] = true;
						for (int i = 0; i < this_node.cover_target_num; i++) {
							int u = this_node.cover_target[i];
							target this_target = all_target_model[u - 1];
							this_target.is_covered[t] = true;
						}
						push_set.remove(p);
					}
				}

			}
		} else {
		}
		return return_slot;
	}

	public static boolean has_remain_query(node this_node, int timeslot, target[] all_target_model) {
		boolean result = false;
		for (int i = 0; i < this_node.cover_target_num; i++) {
			int t = this_node.cover_target[i];
			target this_target = all_target_model[t - 1];
			if (this_target.is_covered[timeslot] == false && this_target.query_number[timeslot] > 0) {
				result = true;
				break;
			}
		}
		return result;
	}

	public static double check_result(int node_num, node[] all_node_model, int target_num, target[] all_target_model,
			int working_duration, int wireless_capacity) {
		int[] last_monitor = new int[target_num];
		for (int i = 0; i < target_num; i++) {
			last_monitor[i] = 0;
		}
		int[] energy = new int[node_num];
		for (int i = 0; i < node_num; i++) {
			energy[i] = 0;
		}
		double sum_aoi = 0;
		double sum_query_num = 0;

		double sum_coverage_time = 0;
		double total_coverage_time = 0;

		double plan_active = 0;
		double total_active = 0;
		for (int i = 1; i <= working_duration; i++) {
			Set<Integer> this_slot_coverage = new HashSet<Integer>();
			int current_capacity = 0;
			for (int j = 1; j <= node_num; j++) {
				node this_node = all_node_model[j - 1];
				if (this_node.active_slots[i] == true) {
					plan_active++;
					if (energy[j - 1] >= 100) {
						total_active++;
						current_capacity += this_node.data_rate;
						for (int k = 0; k < this_node.cover_target_num; k++) {
							this_slot_coverage.add(this_node.cover_target[k]);
						}
						energy[j - 1] -= 100;
						for (int k = 0; k < this_node.cover_target_num; k++) {
							int gh = this_node.cover_target[k];
							last_monitor[gh - 1] = i;
						}
					}

				}
			}
			for (Integer r : this_slot_coverage) {

				target this_t = all_target_model[r - 1];
				sum_coverage_time += (double) (this_t.query_number[i]);
			}
			for (int j = 1; j <= target_num; j++) {
				target this_t = all_target_model[j - 1];
				total_coverage_time += (double) (this_t.query_number[i]);
			}

			if (current_capacity > wireless_capacity) {
				System.out.println("Violate the node number constraint!");
			}
			for (int j = 1; j <= target_num; j++) {
				target this_target = all_target_model[j - 1];
				if (this_target.query_number[i] > 0) {
					sum_query_num += this_target.query_number[i] / this_target.weight;
					sum_aoi += this_target.query_number[i] * (i - last_monitor[j - 1]);
				}
			}
			for (int j = 1; j <= node_num; j++) {
				node this_node = all_node_model[j - 1];
				int sn = energy[j - 1];
				energy[j - 1] = Math.min(this_node.energy_capacity, sn + this_node.actual_eh_rate[i]);
			}
		}
		double result = sum_aoi / sum_query_num;
		return result;
	}

	public static double check_result_2(int node_num, node[] all_node_model, int target_num, target[] all_target_model,
			int working_duration, int wireless_capacity) {
		int[] last_monitor = new int[target_num];
		int total_active = 0;

		double sum_coverage_time = 0;
		double total_coverage_time = 0;

		for (int i = 0; i < target_num; i++) {
			last_monitor[i] = 0;
		}
		double[] energy = new double[node_num];
		for (int i = 0; i < node_num; i++) {
			energy[i] = 0;
		}
		double sum_aoi = 0;
		double sum_query_num = 0;
		for (int i = 1; i <= working_duration; i++) {
			int current_capacity = 0;
			Set<Integer> this_slot_coverage = new HashSet<Integer>();

			for (int j = 1; j <= node_num; j++) {
				node this_node = all_node_model[j - 1];
				if (energy[j - 1] >= 100 && wireless_capacity - current_capacity >= this_node.data_rate) {
					current_capacity += this_node.data_rate;
					total_active++;
					for (int k = 0; k < this_node.cover_target_num; k++) {
						this_slot_coverage.add(this_node.cover_target[k]);
					}
					energy[j - 1] -= 100;
					for (int k = 0; k < this_node.cover_target_num; k++) {
						int gh = this_node.cover_target[k];
						last_monitor[gh - 1] = i;
					}
				}
			}

			for (Integer r : this_slot_coverage) {
				target this_t = all_target_model[r - 1];
				sum_coverage_time += (double) (this_t.query_number[i]);
			}
			for (int j = 1; j <= target_num; j++) {
				target this_t = all_target_model[j - 1];
				total_coverage_time += (double) (this_t.query_number[i]);
			}

			for (int j = 1; j <= target_num; j++) {
				target this_target = all_target_model[j - 1];
				if (this_target.query_number[i] > 0) {
					sum_query_num += this_target.query_number[i] / this_target.weight;
					sum_aoi += this_target.query_number[i] * (i - last_monitor[j - 1]);
				}
			}

			for (int j = 1; j <= node_num; j++) {
				node this_node = all_node_model[j - 1];
				double sn = energy[j - 1];
				energy[j - 1] = Math.min(this_node.energy_capacity, sn + this_node.avg_ehrate);
			}

		}
		return sum_aoi / sum_query_num;
	}

	public static boolean check_has_query(int n_num, int time, node[] all_node_model, target[] all_target_model) {
		boolean result = false;
		node this_node = all_node_model[n_num - 1];
		for (int i = 0; i < this_node.cover_target_num; i++) {
			int g = this_node.cover_target[i];
			target this_target = all_target_model[g - 1];
			if (this_target.query_number[time] > 0) {
				result = true;
				break;
			}
		}
		return result;
	}

	public static boolean check_query_cover(int n_num, int time, node[] all_node_model, target[] all_target_model,
			boolean[] is_covered) {
		boolean result = false;
		node this_node = all_node_model[n_num - 1];
		for (int i = 0; i < this_node.cover_target_num; i++) {
			int g = this_node.cover_target[i];
			target this_target = all_target_model[g - 1];
			if (this_target.query_number[time] > 0 && is_covered[g - 1] == false) {
				result = true;
				break;
			}
		}
		return result;
	}

	public static int check_remain_aoi(int n_num, int time, node[] all_node_model, target[] all_target_model,
			boolean[] is_covered, int[] last_cover_time) {
		int result = 0;
		node this_node = all_node_model[n_num - 1];
		for (int i = 0; i < this_node.cover_target_num; i++) {
			int g = this_node.cover_target[i];
			target this_target = all_target_model[g - 1];
			if (this_target.query_number[time] > 0 && is_covered[g - 1] == false) {
				result += (time - last_cover_time[g - 1]);
			}
		}
		return result;
	}

	public static boolean check_query_cover2(int n_num, int time, node[] all_node_model, target[] all_target_model) {
		boolean result = false;
		node this_node = all_node_model[n_num - 1];
		for (int i = 0; i < this_node.cover_target_num; i++) {
			int g = this_node.cover_target[i];
			target this_target = all_target_model[g - 1];
			if (this_target.query_number[time] > 0) {
				result = true;
				break;
			}
		}
		return result;
	}

	public static double check_result_3_original(int node_num, node[] all_node_model, int target_num,
			target[] all_target_model, int working_duration, int wireless_capacity, int maximum_aoi, double lamuda) {

		MDP_policy_iteration[] all_node_MDP_model = new MDP_policy_iteration[node_num];

		Map<Integer, ArrayList<Integer>> plan_active_slot = new HashMap<Integer, ArrayList<Integer>>();
		for (int i = 1; i <= node_num; i++) {
			plan_active_slot.put(i, new ArrayList<Integer>());
		}

		for (int i = 0; i < node_num; i++) {
			Map<Integer, Double> energy_rate_possibility = new HashMap<Integer, Double>();
			Map<Integer, Integer> energy_rate_num = new HashMap<Integer, Integer>();
			double[] has_query_trans = new double[4];
			node this_node = all_node_model[i];
			all_node_MDP_model[i] = new MDP_policy_iteration(maximum_aoi, this_node.energy_capacity);
			MDP_policy_iteration this_MDP = all_node_MDP_model[i];
			for (int j = 1; j <= working_duration; j++) {
				int rate = this_node.eh_rate[j];
				if (energy_rate_num.containsKey(rate) == false) {
					energy_rate_num.put(rate, 1);
				} else {
					int k = energy_rate_num.get(rate);
					energy_rate_num.put(rate, k + 1);
				}
			}
			for (Integer p : energy_rate_num.keySet()) {
				energy_rate_possibility.put(p, ((double) energy_rate_num.get(p) / (double) working_duration));
			}

			boolean[] has_query = new boolean[working_duration];
			int has_q_slot_num = 0;
			int q_q = 0;
			int q_noq = 0;
			int no_q_slot_num = 0;
			int noq_q = 0;
			int noq_noq = 0;
			for (int j = 1; j <= working_duration - 1; j++) {
				if (check_has_query(i + 1, j, all_node_model, all_target_model) == true) {
					has_query[j] = true;
					has_q_slot_num++;
					if (check_has_query(i + 1, j + 1, all_node_model, all_target_model) == true) {
						q_q++;
					} else {
						q_noq++;
					}
				} else {
					has_query[j] = false;
					no_q_slot_num++;
					if (check_has_query(i + 1, j + 1, all_node_model, all_target_model) == true) {
						noq_q++;
					} else {
						noq_noq++;
					}
				}
			}
			has_query_trans[0] = ((double) q_q / (double) has_q_slot_num);
			has_query_trans[1] = ((double) q_noq / (double) has_q_slot_num);
			has_query_trans[2] = ((double) noq_q / (double) no_q_slot_num);
			has_query_trans[3] = ((double) noq_noq / (double) no_q_slot_num);

			while (true) {
				while (true) {
					double delta = 0;
					for (int p1 = 1; p1 <= maximum_aoi; p1++) {
						for (int p2 = 0; p2 <= this_node.energy_capacity; p2++) {
							for (int p3 = 0; p3 <= 1; p3++) {

								double temp = this_MDP.policy_value[p1][p2][p3];
								int next_aoi = -1;
								int next_energy_base = -1;
								double this_cost = -1;
								if (this_MDP.policy_action[p1][p2][p3] == false) {
									next_aoi = p1 + 1;
									if (next_aoi > maximum_aoi) {
										next_aoi = maximum_aoi;
									}
									next_energy_base = p2;
									if (p3 == 0)
										this_cost = 0;
									else
										this_cost = -p1;
								} else {
									if (p2 >= 100) {
										next_energy_base = p2 - 100;
										next_aoi = 1;
										this_cost = 0;
									} else {
										next_energy_base = p2;
										next_aoi = p1 + 1;
										if (next_aoi > maximum_aoi) {
											next_aoi = maximum_aoi;
										}
										if (p3 == 0)
											this_cost = 0;
										else
											this_cost = -p1;
									}
								}
								double new_value = 0;
								Map<Integer, Double> next_energy_possibility = new HashMap<Integer, Double>();
								for (Integer eh_rate : energy_rate_possibility.keySet()) {
									int next_rate = next_energy_base + eh_rate;
									if (next_rate > this_node.energy_capacity) {
										next_rate = this_node.energy_capacity;
									}
									if (next_energy_possibility.containsKey(next_rate) == false) {
										next_energy_possibility.put(next_rate, energy_rate_possibility.get(eh_rate));
									} else {
										double x = next_energy_possibility.get(next_rate);
										next_energy_possibility.put(next_rate,
												x + energy_rate_possibility.get(eh_rate));
									}
								}
								if (p3 == 0) {
									for (Integer h : next_energy_possibility.keySet()) {
										new_value += has_query_trans[2] * next_energy_possibility.get(h)
												* (this_cost + lamuda * this_MDP.policy_value[next_aoi][h][1]);
										new_value += has_query_trans[3] * next_energy_possibility.get(h)
												* (this_cost + lamuda * this_MDP.policy_value[next_aoi][h][0]);
									}
								} else {
									for (Integer h : next_energy_possibility.keySet()) {
										new_value += has_query_trans[0] * next_energy_possibility.get(h)
												* (this_cost + lamuda * this_MDP.policy_value[next_aoi][h][1]);
										new_value += has_query_trans[1] * next_energy_possibility.get(h)
												* (this_cost + lamuda * this_MDP.policy_value[next_aoi][h][0]);
									}
								}
								this_MDP.policy_value[p1][p2][p3] = new_value;
								double hj = Math.abs(temp - new_value);
								delta = Math.max(hj, delta);
							}
						}
					}
					if (delta <= 0.01) {
						break;
					}
				}

				boolean policy_stable = true;
				for (int p1 = 1; p1 <= maximum_aoi; p1++) {
					for (int p2 = 0; p2 <= this_node.energy_capacity; p2++) {
						for (int p3 = 0; p3 <= 1; p3++) {

							boolean temp = this_MDP.policy_action[p1][p2][p3];

							double v_active = 0;
							int next_aoi = -1;
							int next_energy_base = -1;
							double this_cost = -1;
							if (p2 >= 100) {
								next_energy_base = p2 - 100;
								next_aoi = 1;
								this_cost = 0;
							} else {
								next_energy_base = p2;
								next_aoi = p1 + 1;
								if (next_aoi > maximum_aoi) {
									next_aoi = maximum_aoi;
								}
								if (p3 == 0)
									this_cost = 0;
								else
									this_cost = -p1;
							}
							Map<Integer, Double> next_energy_possibility = new HashMap<Integer, Double>();
							for (Integer eh_rate : energy_rate_possibility.keySet()) {
								int next_rate = next_energy_base + eh_rate;
								if (next_rate > this_node.energy_capacity) {
									next_rate = this_node.energy_capacity;
								}
								if (next_energy_possibility.containsKey(next_rate) == false) {
									next_energy_possibility.put(next_rate, energy_rate_possibility.get(eh_rate));
								} else {
									double x = next_energy_possibility.get(next_rate);
									next_energy_possibility.put(next_rate, x + energy_rate_possibility.get(eh_rate));
								}
							}
							if (p3 == 0) {
								for (Integer h : next_energy_possibility.keySet()) {
									v_active += has_query_trans[2] * next_energy_possibility.get(h)
											* (this_cost + lamuda * this_MDP.policy_value[next_aoi][h][1]);
									v_active += has_query_trans[3] * next_energy_possibility.get(h)
											* (this_cost + lamuda * this_MDP.policy_value[next_aoi][h][0]);
								}
							} else {
								for (Integer h : next_energy_possibility.keySet()) {
									v_active += has_query_trans[0] * next_energy_possibility.get(h)
											* (this_cost + lamuda * this_MDP.policy_value[next_aoi][h][1]);
									v_active += has_query_trans[1] * next_energy_possibility.get(h)
											* (this_cost + lamuda * this_MDP.policy_value[next_aoi][h][0]);
								}
							}

							double v_idle = 0;
							next_aoi = -1;
							next_energy_base = -1;
							this_cost = -1;
							if (p2 >= 100) {
								next_energy_base = p2 - 100;
								next_aoi = 1;
								this_cost = 0;
							} else {
								next_energy_base = p2;
								next_aoi = p1 + 1;
								if (next_aoi > maximum_aoi) {
									next_aoi = maximum_aoi;
								}
								if (p3 == 0)
									this_cost = 0;
								else
									this_cost = -p1;
							}
							next_energy_possibility = new HashMap<Integer, Double>();
							for (Integer eh_rate : energy_rate_possibility.keySet()) {
								int next_rate = next_energy_base + eh_rate;
								if (next_rate > this_node.energy_capacity) {
									next_rate = this_node.energy_capacity;
								}
								if (next_energy_possibility.containsKey(next_rate) == false) {
									next_energy_possibility.put(next_rate, energy_rate_possibility.get(eh_rate));
								} else {
									double x = next_energy_possibility.get(next_rate);
									next_energy_possibility.put(next_rate, x + energy_rate_possibility.get(eh_rate));
								}
							}
							if (p3 == 0) {
								for (Integer h : next_energy_possibility.keySet()) {
									v_idle += has_query_trans[2] * next_energy_possibility.get(h)
											* (this_cost + lamuda * this_MDP.policy_value[next_aoi][h][1]);
									v_idle += has_query_trans[3] * next_energy_possibility.get(h)
											* (this_cost + lamuda * this_MDP.policy_value[next_aoi][h][0]);
								}
							} else {
								for (Integer h : next_energy_possibility.keySet()) {
									v_idle += has_query_trans[0] * next_energy_possibility.get(h)
											* (this_cost + lamuda * this_MDP.policy_value[next_aoi][h][1]);
									v_idle += has_query_trans[1] * next_energy_possibility.get(h)
											* (this_cost + lamuda * this_MDP.policy_value[next_aoi][h][0]);
								}
							}
							if (v_idle > v_active) {
								this_MDP.policy_action[p1][p2][p3] = false;
							} else {
								this_MDP.policy_action[p1][p2][p3] = true;
							}
							if (temp != this_MDP.policy_action[p1][p2][p3]) {
								policy_stable = false;
							}
						}
					}
				}
				if (policy_stable == true) {
					break;
				}

			}
		
			for (int p1 = 1; p1 <= maximum_aoi; p1++) {
				for (int p2 = 0; p2 <= this_node.energy_capacity; p2++) {
					for (int p3 = 0; p3 <= 1; p3++) {
						if (this_MDP.policy_action[p1][p2][p3] == true) {
							if (p3 == 0) {
								this_MDP.policy_action[p1][p2][p3] = false;
							}
						}
					}
				}
			}
		}

		int[] last_monitor = new int[target_num];
		int[] last_active = new int[node_num];

		for (int i = 0; i < target_num; i++) {
			last_monitor[i] = 0;
		}
		int[] energy = new int[node_num];
		for (int i = 0; i < node_num; i++) {
			energy[i] = 0;
			last_active[i] = 0;
		}

		for (int i = 1; i <= working_duration; i++) {
			int current_capacity = 0;
			for (int j = 1; j <= node_num; j++) {
				node this_node = all_node_model[j - 1];
				MDP_policy_iteration this_node_MDP_model = all_node_MDP_model[j - 1];
				if (energy[j - 1] >= 100 && wireless_capacity - current_capacity >= this_node.data_rate) {
					int b1 = energy[j - 1];
					int b2 = i - last_active[j - 1];
					if (b2 > maximum_aoi) {
						b2 = maximum_aoi;
					}
					int b3 = 0;
					if (check_has_query(j, i, all_node_model, all_target_model) == true) {
						b3 = 1;
					}
					boolean ac = this_node_MDP_model.policy_action[b2][b1][b3];
					if (ac == true) {
						plan_active_slot.get(j).add(i);
						current_capacity += this_node.data_rate;
						energy[j - 1] -= 100;
						last_active[j - 1] = i;
						for (int k = 0; k < this_node.cover_target_num; k++) {
							int gh = this_node.cover_target[k];
							last_monitor[gh - 1] = i;
						}
					}
				}
			}
			for (int j = 1; j <= node_num; j++) {
				node this_node = all_node_model[j - 1];
				int sn = energy[j - 1];
				energy[j - 1] = Math.min(this_node.energy_capacity, sn + this_node.eh_rate[i]);
			}
		}

		double sum_aoi = 0;
		double sum_query_num = 0;
		for (int i = 0; i < target_num; i++) {
			last_monitor[i] = 0;
		}
		for (int i = 0; i < node_num; i++) {
			energy[i] = 0;
			last_active[i] = 0;
		}

		boolean[] is_covered = new boolean[target_num];

		for (int i = 1; i <= working_duration; i++) {
			for (int k = 0; k < target_num; k++) {
				is_covered[k] = false;
			}
			int current_capacity = 0;
			for (int j = 1; j <= node_num; j++) {
				node this_node = all_node_model[j - 1];
				if (energy[j - 1] >= 100 && wireless_capacity - current_capacity >= this_node.data_rate) {
					if (plan_active_slot.get(j).contains(i) == true) {
						current_capacity += this_node.data_rate;
						energy[j - 1] -= 100;
						last_active[j - 1] = i;
						for (int k = 0; k < this_node.cover_target_num; k++) {
							int gh = this_node.cover_target[k];
							is_covered[gh - 1] = true;
							last_monitor[gh - 1] = i;
						}
					}
				}
			}
			for (int j = 1; j <= target_num; j++) {
				target this_target = all_target_model[j - 1];
				if (this_target.query_number[i] > 0) {
					sum_query_num += this_target.query_number[i];
					sum_aoi += this_target.query_number[i] * (i - last_monitor[j - 1]);
				}
			}
			for (int j = 1; j <= node_num; j++) {
				node this_node = all_node_model[j - 1];
				int sn = energy[j - 1];
				energy[j - 1] = Math.min(this_node.energy_capacity, sn + this_node.actual_eh_rate[i]);
			}
		}

		return sum_aoi / sum_query_num;
	}

	public static double check_result_greedy(int node_num, node[] all_node_model, int target_num,
			target[] all_target_model, int working_duration, int wireless_capacity, int maximum_aoi, double lamuda) {
		Map<Integer, ArrayList<Integer>> plan_active_slot = new HashMap<Integer, ArrayList<Integer>>();
		for (int i = 1; i <= node_num; i++) {
			plan_active_slot.put(i, new ArrayList<Integer>());
		}

		int[] last_monitor = new int[target_num];
		int[] last_active = new int[node_num];
		boolean[] is_covered = new boolean[target_num];
		for (int i = 0; i < target_num; i++) {
			last_monitor[i] = 0;
		}
		int[] energy = new int[node_num];
		for (int i = 0; i < node_num; i++) {
			energy[i] = 0;
			last_active[i] = 0;
		}

		double[] result = new double[3];

		for (int i = 1; i <= working_duration; i++) {
			for (int k = 0; k < target_num; k++) {
				is_covered[k] = false;
			}
			int current_capacity = 0;
			while (true) {
				int candidate_node_num = -1;
				int candidate_cover_aoi = 0;
				for (int j = 1; j <= node_num; j++) {
					node this_node = all_node_model[j - 1];
					if (energy[j - 1] >= 100 && wireless_capacity - current_capacity >= this_node.data_rate) {
						int ra = check_remain_aoi(j, i, all_node_model, all_target_model, is_covered, last_monitor);
						if (ra > candidate_cover_aoi) {
							candidate_node_num = j;
							candidate_cover_aoi = ra;
						}
					}
				}

				if (candidate_node_num != -1) {
					node this_node = all_node_model[candidate_node_num - 1];
					current_capacity += this_node.data_rate;
					energy[candidate_node_num - 1] -= 100;
					last_active[candidate_node_num - 1] = i;
					plan_active_slot.get(candidate_node_num).add(i);
					for (int k = 0; k < this_node.cover_target_num; k++) {
						int gh = this_node.cover_target[k];
						is_covered[gh - 1] = true;
						last_monitor[gh - 1] = i;
					}
				} else {
					break;
				}
			}
			for (int j = 1; j <= node_num; j++) {
				node this_node = all_node_model[j - 1];
				int sn = energy[j - 1];
				energy[j - 1] = Math.min(this_node.energy_capacity, sn + this_node.eh_rate[i]);
			}

		}

		double sum_aoi = 0;
		double sum_query_num = 0;
		double cover_query_num = 0;
		for (int i = 0; i < target_num; i++) {
			last_monitor[i] = 0;
		}
		for (int i = 0; i < node_num; i++) {
			energy[i] = 0;
			last_active[i] = 0;
		}
		for (int i = 1; i <= working_duration; i++) {
			for (int k = 0; k < target_num; k++) {
				is_covered[k] = false;
			}
			int current_capacity = 0;
			for (int j = 1; j <= node_num; j++) {
				node this_node = all_node_model[j - 1];
				if (energy[j - 1] >= 100 && wireless_capacity - current_capacity >= this_node.data_rate) {
					if (plan_active_slot.get(j).contains(i) == true) {
						current_capacity += this_node.data_rate;
						energy[j - 1] -= 100;
						last_active[j - 1] = i;
						for (int k = 0; k < this_node.cover_target_num; k++) {
							int gh = this_node.cover_target[k];
							is_covered[gh - 1] = true;
							last_monitor[gh - 1] = i;
						}
					}
				}
			}
			for (int j = 1; j <= target_num; j++) {
				target this_target = all_target_model[j - 1];
				if (this_target.query_number[i] > 0) {
					sum_query_num += this_target.query_number[i];
					sum_aoi += this_target.query_number[i] * (i - last_monitor[j - 1]);
					if (is_covered[j - 1] == true) {
						cover_query_num += this_target.query_number[i];
					}
				}
			}
			for (int j = 1; j <= node_num; j++) {
				node this_node = all_node_model[j - 1];
				int sn = energy[j - 1];
				energy[j - 1] = Math.min(this_node.energy_capacity, sn + this_node.actual_eh_rate[i]);
			}
		}

		return sum_aoi / sum_query_num;
	}

	static double[] read_original_data(String data_folder, int day, int slot, int duration_length) {
		double[] result = new double[duration_length];
		String data_file = data_folder + "\\" + String.valueOf(day);
		try (BufferedReader br = new BufferedReader(new FileReader(data_file))) {
			String line;
			int i = -1;
			int line_num = 0;
			while ((line = br.readLine()) != null) {
				line_num++;
				if (line_num >= 1 + (slot - 1) * duration_length && line_num <= slot * duration_length) {
					try {
						double number = Double.parseDouble(line);
						i++;
						result[i] = number;
					} catch (NumberFormatException e) {
						System.err.println("Invalid number: " + line);
					}
				}
			}
		} catch (IOException e) {
			System.err.println("Error reading file: " + e.getMessage());
		}
		return result;
	}

	static void output_network_topology(int node_num, node[] all_node_model, int target_num, target[] all_target_model,
			int duration_length, int wireless_bandwidth, int energy_buffer) throws IOException {
		File file1 = new File("D:\\expdata\\network_topology.txt");
		BufferedWriter writer = new BufferedWriter(new FileWriter(file1));
		for (int i = 0; i < node_num; i++) {
			node this_node = all_node_model[i];
			String x = String.valueOf(i + 1);
			x = x + ":";
			for (int j = 0; j < this_node.cover_target_num; j++) {
				int g = this_node.cover_target[j];
				if (j < this_node.cover_target_num - 1) {
					x = x + String.valueOf(g);
					x = x + ",";
				} else {
					x = x + String.valueOf(g);
				}
			}
			writer.write(x);
			if (i != node_num - 1) {
				writer.newLine();
			}
		}
		writer.close();

		File file2 = new File("D:\\expdata\\energy_data.txt");
		writer = new BufferedWriter(new FileWriter(file2));
		for (int i = 0; i < node_num; i++) {
			node this_node = all_node_model[i];
			String x = String.valueOf(i + 1);
			x = x + ":";
			for (int j = 1; j <= duration_length; j++) {
				int g = this_node.eh_rate[j];
				if (j < duration_length) {
					x = x + String.valueOf(g);
					x = x + ",";
				} else {
					x = x + String.valueOf(g);
				}
			}
			writer.write(x);
			if (i != node_num - 1) {
				writer.newLine();
			}
		}
		writer.close();

		File file3 = new File("D:\\expdata\\target_query.txt");
		writer = new BufferedWriter(new FileWriter(file3));
		for (int i = 0; i < target_num; i++) {
			target this_target = all_target_model[i];
			String x = String.valueOf(i + 1);
			x = x + ":";
			for (int j = 1; j <= duration_length; j++) {
				int g = this_target.query_number[j];
				if (j < duration_length) {
					x = x + String.valueOf(g);
					x = x + ",";
				} else {
					x = x + String.valueOf(g);
				}
			}
			writer.write(x);
			if (i != target_num - 1) {
				writer.newLine();
			}
		}
		writer.close();

		File file4 = new File("D:\\expdata\\target_weight.txt");
		writer = new BufferedWriter(new FileWriter(file4));
		for (int i = 0; i < target_num; i++) {
			target this_target = all_target_model[i];
			String x = String.valueOf(i + 1);
			x = x + ":";
			x = x + String.valueOf(this_target.weight);
			writer.write(x);
			if (i != target_num - 1) {
				writer.newLine();
			}
		}
		writer.close();

		File file5 = new File("D:\\expdata\\restrict.txt");
		writer = new BufferedWriter(new FileWriter(file5));
		writer.write(String.valueOf(wireless_bandwidth));
		writer.newLine();
		writer.write(String.valueOf(energy_buffer));

		writer.close();

	}

	static Map<Integer, ArrayList<Double>> read_average_energy_data(String folder_path) {
		Map<Integer, ArrayList<Double>> result = new HashMap<Integer, ArrayList<Double>>();
		File folder = new File(folder_path);
		if (folder.exists() && folder.isDirectory()) {

			File[] files = folder.listFiles();
			for (File file : files) {
				if (file.isFile()) {
					result.put(Integer.parseInt(file.getName()), new ArrayList<Double>());
					String data_file = folder_path + "\\" + file.getName();
					try (BufferedReader br = new BufferedReader(new FileReader(data_file))) {
						String line;
						while ((line = br.readLine()) != null) {
							try {
								double number = Double.parseDouble(line);
								result.get(Integer.parseInt(file.getName())).add(number);
							} catch (NumberFormatException e) {
								System.err.println("Invalid number: " + line);
							}
						}
					} catch (IOException e) {
						System.err.println("Error reading file: " + e.getMessage());
					}

				}
			}
		}
		return result;
	}

	static int find_similar_day(int day_num, int slot_num, int max_k,
			Map<Integer, ArrayList<Double>> energy_avg_value) {
		int result = -1;
		double min_diff = 999999;
		for (Integer r : energy_avg_value.keySet()) {
			if (r < day_num) {
				double sum_diff = 0;
				for (int i = slot_num - max_k; i <= slot_num - 1; i++) {
					sum_diff += Math.abs(energy_avg_value.get(r).get(i) - energy_avg_value.get(day_num).get(i));
				}
				if (sum_diff < min_diff) {
					min_diff = sum_diff;
					result = r;
				}
			}
		}

		return result;
	}

	static Map<Integer, Map<Integer, Integer>> compute_similar_day(int start_day, int end_day, int start_slot,
			int end_slot, Map<Integer, ArrayList<Double>> energy_avg_value, int max_k) {
		Map<Integer, Map<Integer, Integer>> result = new HashMap<Integer, Map<Integer, Integer>>();

		for (Integer d : energy_avg_value.keySet()) {
			if (d >= start_day && d <= end_day) {
				result.put(d, new HashMap<Integer, Integer>());
				for (int i = start_slot; i <= end_slot; i++) {
					int k = find_similar_day(d, i, max_k, energy_avg_value);
					result.get(d).put(i, k);
				}
			}
		}

		return result;
	}

	static int correspond_ehrate(double active_energy, double original_energy) {
		double x = original_energy / active_energy;
		int h = (int) (Math.floor(100 * x));
		if (h == 0) {
			h = 1;
		}
		return h;
	}

	static void give_energy_sum(node this_node, int sum_energy) {
		int remain = sum_energy;
		for (int i = 1; i <= 1000; i++) {
			this_node.eh_rate[i] = 0;
		}
		while (true) {
			if (remain >= 1000) {
				remain -= 1000;
				for (int i = 1; i <= 1000; i++) {
					this_node.eh_rate[i]++;
				}
			} else if (remain >= 500) {
				remain -= 500;
				for (int i = 1; i <= 500; i++) {
					this_node.eh_rate[i * 2]++;
				}
			} else if (remain >= 250) {
				remain -= 250;
				for (int i = 1; i <= 250; i++) {
					this_node.eh_rate[i * 4 - 1]++;
				}
			} else if (remain >= 125) {
				remain -= 125;
				for (int i = 1; i <= 125; i++) {
					this_node.eh_rate[i * 8 - 3]++;
				}
			} else if (remain >= 62) {
				remain -= 62;
				for (int i = 1; i <= 62; i++) {
					this_node.eh_rate[i * 16 - 7]++;
				}
			} else if (remain >= 31) {
				remain -= 31;
				for (int i = 1; i <= 31; i++) {
					this_node.eh_rate[i * 32 - 15]++;
				}
			} else if (remain >= 16) {
				remain -= 16;
				for (int i = 1; i <= 16; i++) {
					this_node.eh_rate[i * 64 - 31]++;
				}
			} else if (remain >= 8) {
				remain -= 8;
				for (int i = 1; i <= 8; i++) {
					this_node.eh_rate[i * 128 - 63]++;
				}
			} else if (remain >= 4) {
				remain -= 4;
				for (int i = 1; i <= 4; i++) {
					this_node.eh_rate[i * 256 - 127]++;
				}
			} else if (remain >= 2) {
				remain -= 2;
				for (int i = 1; i <= 2; i++) {
					this_node.eh_rate[i * 512 - 255]++;
				}
			} else if (remain >= 1) {
				remain -= 1;
				this_node.eh_rate[513]++;
			} else if (remain == 0) {
				break;
			}
		}
	}

	static void compute_predict_energy_value(node this_node, double ratio, int day, int slot_num, double active_energy,
			double alpha, Map<Integer, Map<Integer, Integer>> similar_day, String original_data_folder,
			int duration_length) {
		int s_day = similar_day.get(day).get(slot_num);
		double[] similar_day_original_data = read_original_data(original_data_folder, s_day, slot_num, duration_length);
		int similar_day_sum_energy = 0;
		for (int i = 0; i < duration_length; i++) {
			similar_day_sum_energy += correspond_ehrate(active_energy, ratio * similar_day_original_data[i]);
		}
		double[] this_day_former_data = read_original_data(original_data_folder, day, slot_num - 1, duration_length);
		int former_day_sum_energy = 0;
		for (int i = 0; i < duration_length; i++) {
			former_day_sum_energy += correspond_ehrate(active_energy, ratio * this_day_former_data[i]);
		}
		double predict_value = alpha * ((double) former_day_sum_energy)
				+ (1 - alpha) * ((double) similar_day_sum_energy);
		int sum_predict_value = (int) (Math.floor(predict_value));

		give_energy_sum(this_node, sum_predict_value);

		double[] day_original_data = read_original_data(original_data_folder, day, slot_num, duration_length);
		double sum_energy = 0;
		for (int i = 0; i < duration_length; i++) {
			int u = correspond_ehrate(active_energy, ratio * day_original_data[i]);
			this_node.actual_eh_rate[i + 1] = u;
			sum_energy += (double) u;
		}
		this_node.avg_ehrate = sum_energy / ((double) duration_length);
	}

	static int node_num = 40;
	static int target_num = 52;
	static int pkt_header = 12;
	static int coverage_range = 24;

	static String energy_data_file = "D:\\\\indoorEHdatabase\\\\2019_08_processed\\\\processed\\\\data1119";

	static int wireless_capacity = 10;
	static double average_cover_num = 7.8;
	static int working_duration = 1000;
	static int avg_query_interval = 60;
	static int avg_query_length = 4;
	static int energy_buffer = 300;
	static int day_start = 100;
	static int day_end = 199;
	static int slot_start = 4;
	static int slot_end = 10;

	static int evaluate_slot = 4; // The energy data of which period. 
	
	static String energy_dataset_by_day = ""; //The path of folder "data1119", the raw energy data 
	static String energy_dataset_by_day_average = ""; //The path of folder "datasummary1119", the average energy trace data of each day

	public static void main(String args[]){

		boolean has_weight = false;
		boolean show_infor = true;

		double aoi_0 = 0;
		double aoi_tmc = 0;
		double aoi_greedy = 0;
		double aoi_tcom = 0;
		
		Map<Integer, ArrayList<Double>> all_avg_original_data = read_average_energy_data(
				energy_dataset_by_day_average);
		Map<Integer, Map<Integer, Integer>> most_similar_day = compute_similar_day(day_start, day_end, slot_start,
				slot_end, all_avg_original_data, 4);

		int times = 0;
		
		/*
			Each data point is the average of the results from the simulation conducted 
			on the energy trace sampled starting from 01/04/2018 to 31/07/2018 (100 days)
		*/
		for (Integer day : all_avg_original_data.keySet()) {
			if (day >= day_start && day <= day_end) {
				System.out.println("Execute simulation on day number:" + day);
				times++;
				node[] all_node_model = new node[node_num];
				for (int i = 0; i < node_num; i++) {
					all_node_model[i] = new node(working_duration, energy_buffer);
				}
				target[] all_target_model = generate_all_target_model(target_num, working_duration, has_weight);
				deploy_node_target(node_num, all_node_model, target_num, all_target_model, coverage_range, pkt_header,
						average_cover_num);

				Random random = new Random();
				for (int i = 0; i < node_num; i++) {
					node ui = all_node_model[i];
					double energy_c = ((48 * (12 + ((double) (ui.cover_target_num))) / 250000) * 64.5 * 0.001)
							/ 0.000001;
					double energy_ratio = 0.6 + 0.8 * random.nextDouble();

					compute_predict_energy_value(ui, energy_ratio, day, evaluate_slot, energy_c, 0.2,
							most_similar_day, energy_dataset_by_day,
							working_duration);

				}

				generate_query(node_num, all_node_model, target_num, all_target_model, avg_query_interval,
						avg_query_length, working_duration);

				// output_network_topology(node_num, all_node_model, target_num,
				// all_target_model, working_duration,
				// wireless_capacity, energy_buffer);

				// Run the algorithm.
				long startTime = System.currentTimeMillis();
				int[] current_capacity = new int[working_duration + 1];
				for (int i = 1; i <= working_duration; i++) {
					current_capacity[i] = 0;
				}

				ArrayList<Integer> initial_schedule_node = new ArrayList<Integer>();
				Set<Integer> all_remain_node = new HashSet<Integer>();
				Set<Integer> all_remain_target = new HashSet<Integer>();
				for (int i = 0; i < node_num; i++) {
					all_remain_node.add(i + 1);
				}
				for (int i = 0; i < target_num; i++) {
					all_remain_target.add(i + 1);
				}
				while (true) {
					int cover_n = 0;
					int cover_node = -1;
					Set p = new HashSet<Integer>();
					for (Integer j : all_remain_node) {
						Set r = intersect(all_node_model[j - 1].cover_target, all_node_model[j - 1].cover_target_num,
								all_remain_target);
						if (r.size() > cover_n) {
							p = r;
							cover_node = j;
							cover_n = r.size();
						}
					}
					initial_schedule_node.add(cover_node);
					all_remain_node.remove(cover_node);
					all_remain_target.removeAll(p);
					if (all_remain_target.size() == 0) {
						break;
					}
				}

				int ini_time = 0;
				Map<Integer, Integer> consider_set = new HashMap<Integer, Integer>();
				while (true) {
					ini_time++;
					if (ini_time == 1) {
						int bo = 0;
						Map<Integer, Integer> active_slot_re = new HashMap<Integer, Integer>();
						for (int i = 0; i < initial_schedule_node.size(); i++) {
							int ui = initial_schedule_node.get(i);
							int ca = all_node_model[ui - 1].data_rate;
							if (bo == 0) {
								bo = 1;
								active_slot_re.put(ui, 100 + bo);
								current_capacity[100 + bo] += ca;
							} else {
								int can_box = 0;
								for (int j = 1; j <= bo; j++) {
									if (wireless_capacity - current_capacity[100 + j] >= ca) {
										can_box = j;
									}
								}
								if (can_box == 0) {
									bo++;
									can_box = bo;
									active_slot_re.put(ui, 100 + can_box);
									current_capacity[100 + can_box] += ca;
								} else {
									active_slot_re.put(ui, 100 + can_box);
									current_capacity[100 + can_box] += ca;
								}
							}
						}
						Map<Integer, Integer> return_n = push_back(all_node_model, all_target_model, active_slot_re,
								current_capacity, working_duration, wireless_capacity);
						consider_set.clear();
						for (Integer p : return_n.keySet()) {
							consider_set.put(p, return_n.get(p));
						}
						if (consider_set.size() == 0) {
							break;
						}
					} else {
						Map<Integer, Integer> active_slot_re = new HashMap<Integer, Integer>();
						for (Integer x : consider_set.keySet()) {
							node this_node = all_node_model[x - 1];
							if (consider_set.get(x) + 100 <= working_duration) {
								active_slot_re.put(x, consider_set.get(x) + 100);
								current_capacity[consider_set.get(x) + 100] += this_node.data_rate;
							}
						}
						if (active_slot_re.keySet().size() == 0) {
							break;
						} else {
							Map<Integer, Integer> return_n = push_back(all_node_model, all_target_model, active_slot_re,
									current_capacity, working_duration, wireless_capacity);
							consider_set.clear();
							for (Integer p : return_n.keySet()) {
								consider_set.put(p, return_n.get(p));
							}
							if (consider_set.size() == 0) {
								break;
							}
						}
					}
				}

				int[] target_last_monitor = new int[target_num];
				for (int i = 0; i < target_num; i++) {
					target_last_monitor[i] = 0;
				}

				for (int i = 1; i <= working_duration; i++) {
					for (int j = 0; j < node_num; j++) {
						node this_node = all_node_model[j];
						if (this_node.active_slots[i] == true) {
							for (int k = 0; k < this_node.cover_target_num; k++) {
								int h = this_node.cover_target[k];
								target_last_monitor[h - 1] = i;
							}
						}
					}
					for (int j = 0; j < target_num; j++) {
						target this_target = all_target_model[j];
						if (this_target.query_number[i] > 0) {
							int up = i - target_last_monitor[j];
							all_target_model[j].current_aoi_value[i] = up;
						}
					}
				}

				Map<Integer, Integer> schedule_sequence = new HashMap<Integer, Integer>();
				for (int i = 0; i < node_num; i++) {
					node this_node = all_node_model[i];
					schedule_sequence.put(i + 1, this_node.cover_target_num);
				}
				schedule_sequence = test_compare.valueUpSort(schedule_sequence);

				for (Integer i : schedule_sequence.keySet()) {

					node this_node = all_node_model[i - 1];
					int current_maximum_age = -1;
					for (int j = 0; j < this_node.cover_target_num; j++) {
						int gp = this_node.cover_target[j];
						target this_target = all_target_model[gp - 1];
						for (int k = 1; k <= working_duration; k++) {
							if (this_target.query_number[k] > 0) {
								if (this_target.current_aoi_value[k] > current_maximum_age) {
									current_maximum_age = this_target.current_aoi_value[k];
								}
							}
						}
					}
					if (current_maximum_age > 0) {
						int[] last_active_distance = new int[working_duration + 1];
						int most_recent_active = 0;
						for (int j = 1; j <= working_duration; j++) {
							last_active_distance[j] = j - most_recent_active;
							if (this_node.active_slots[j] == true) {
								most_recent_active = j;
							}
						}
						ArrayList<target> cover_target_add = new ArrayList<target>();
						for (int aj = 0; aj < this_node.cover_target_num; aj++) {
							int pk = this_node.cover_target[aj];
							cover_target_add.add(all_target_model[pk - 1]);
						}
						int[][][] cost = new int[working_duration + 1][current_maximum_age
								+ 1][this_node.energy_capacity + 1];
						boolean[][][] action = new boolean[working_duration + 1][current_maximum_age
								+ 1][this_node.energy_capacity + 1];

						for (int w = working_duration; w >= 1; w--) {
							int remain_capacity = wireless_capacity - current_capacity[w];
							int max_inter = Math.min(last_active_distance[w], current_maximum_age);

							if (w == working_duration) {
								for (int x = 1; x <= max_inter; x++) {
									int this_cost = 0;
									for (target th : cover_target_add) {
										if (th.query_number[w] > 0) {
											if (th.current_aoi_value[w] > x) {
												this_cost += th.query_number[w] * x;
											} else {
												this_cost += th.query_number[w] * th.current_aoi_value[w];
											}
										} else {

										}
									}
									for (int y = 0; y <= this_node.energy_capacity; y++) {
										if (this_node.active_slots[w] == true) {
											if (y < 100) {
												cost[w][x][y] = -1;
											} else {
												cost[w][x][y] = 0;
												action[w][x][y] = true;
											}
										} else {
											if (y < 100 || this_node.data_rate > remain_capacity) {
												cost[w][x][y] = this_cost;
												action[w][x][y] = false;
											} else {
												if (this_cost == 0) {
													cost[w][x][y] = this_cost;
													action[w][x][y] = false;
												} else {
													cost[w][x][y] = 0;
													action[w][x][y] = true;
												}
											}
										}

									}
								}
								for (int x = max_inter + 1; x <= current_maximum_age; x++) {
									for (int y = 0; y <= this_node.energy_capacity; y++) {
										cost[w][x][y] = -1;
									}
								}
							} else {
								int[] next_energy_idle = new int[this_node.energy_capacity + 1];
								int[] next_energy_active = new int[this_node.energy_capacity + 1];
								for (int aj = 0; aj <= this_node.energy_capacity; aj++) {
									next_energy_idle[aj] = Math.min(aj + this_node.eh_rate[w],
											this_node.energy_capacity);

									if (aj >= 100) {
										next_energy_active[aj] = Math.min(aj - 100 + this_node.eh_rate[w],
												this_node.energy_capacity);
									}
								}
								for (int x = 1; x <= max_inter; x++) {
									int next_interval = -1;
									if (x < current_maximum_age) {
										next_interval = x + 1;
									} else {
										next_interval = current_maximum_age;
									}
									int this_cost = 0;
									for (target th : cover_target_add) {
										if (th.query_number[w] > 0) {
											if (th.current_aoi_value[w] > x) {
												this_cost += th.query_number[w] * x;
											} else {
												this_cost += th.query_number[w] * th.current_aoi_value[w];
											}
										} else {

										}
									}
									for (int y = 0; y <= this_node.energy_capacity; y++) {
										if (this_node.active_slots[w] == true) {
											if (y < 100) {
												cost[w][x][y] = -1;
											} else {
												int next_energy = next_energy_active[y];
												if (cost[w + 1][1][next_energy] == -1) {
													cost[w][x][y] = -1;
												} else {
													cost[w][x][y] = cost[w + 1][1][next_energy];
													action[w][x][y] = true;
												}
											}
										} else {
											if (y < 100 || this_node.data_rate > remain_capacity) {
												int next_energy = next_energy_idle[y];
												if (cost[w + 1][next_interval][next_energy] == -1) {
													cost[w][x][y] = -1;
												} else {
													cost[w][x][y] = this_cost + cost[w + 1][next_interval][next_energy];
													action[w][x][y] = false;
												}
											} else {
												int cost_idle;
												int cost_active;
												int next_energy = next_energy_idle[y];
												if (cost[w + 1][next_interval][next_energy] == -1) {
													cost_idle = -1;
												} else {
													cost_idle = this_cost + cost[w + 1][next_interval][next_energy];
												}
												next_energy = next_energy_active[y];
												if (cost[w + 1][1][next_energy] == -1) {
													cost_active = -1;
												} else {
													cost_active = cost[w + 1][1][next_energy];
												}
												if (cost_idle != -1 && cost_active != -1) {
													if (cost_active < cost_idle) {
														cost[w][x][y] = cost_active;
														action[w][x][y] = true;
													} else {
														cost[w][x][y] = cost_idle;
														action[w][x][y] = false;
													}
												} else if (cost_idle == -1 && cost_active != -1) {
													cost[w][x][y] = cost_active;
													action[w][x][y] = true;
												} else if (cost_idle != -1 && cost_active == -1) {
													cost[w][x][y] = cost_idle;
													action[w][x][y] = false;
												} else {
													cost[w][x][y] = -1;
												}
											}
										}
									}
								}
								for (int x = max_inter + 1; x <= current_maximum_age; x++) {
									for (int y = 0; y <= this_node.energy_capacity; y++) {
										cost[w][x][y] = -1;
									}
								}
							}
						}

						ArrayList<Integer> new_active_slot = new ArrayList<Integer>();
						int last_interval = 1;
						int current_energy = 0;

						for (int af = 1; af <= working_duration; af++) {
							if (action[af][last_interval][current_energy] == true) {
								last_interval = 1;
								int g = current_energy;
								current_energy = Math.min(g - 100 + this_node.eh_rate[af], this_node.energy_capacity);
								if (this_node.active_slots[af] == false) {
									new_active_slot.add(af);
								}
							} else {
								if (last_interval == current_maximum_age) {
									last_interval = current_maximum_age;
								} else {
									last_interval++;
								}
								int g = current_energy;
								current_energy = Math.min(g + this_node.eh_rate[af], this_node.energy_capacity);
							}
						}

						for (int af = 0; af < new_active_slot.size(); af++) {
							int fg = new_active_slot.get(af);
							current_capacity[fg] += this_node.data_rate;
							this_node.active_slots[fg] = true;
							for (target q : cover_target_add) {
								q.is_covered[fg] = true;
							}
						}
						int ne = 0;
						for (int af = 1; af <= working_duration; af++) {
							if (this_node.active_slots[af] == true) {
								ne = af;
							}
							for (target q : cover_target_add) {
								if (q.query_number[af] > 0) {
									if (q.current_aoi_value[af] > af - ne) {
										q.current_aoi_value[af] = af - ne;
									}
								}
							}
						}
					}
				}

				double r_02 = check_result(node_num, all_node_model, target_num, all_target_model, working_duration,
						wireless_capacity);

				double r_tmc = check_result_2(node_num, all_node_model, target_num, all_target_model, working_duration,
						wireless_capacity);

				double r_greedy = check_result_greedy(node_num, all_node_model, target_num, all_target_model,
						working_duration, wireless_capacity, 20, 0.95);

				double r_tcom = check_result_3_original(node_num, all_node_model, target_num, all_target_model,
						working_duration, wireless_capacity, 20, 0.95);

				aoi_0 += r_02;
				aoi_tmc += r_tmc;
				aoi_greedy += r_greedy;
				aoi_tcom += r_tcom;

			}
		}
		System.out.println("The result: ");
		System.out.println();
		System.out.println("The average AoI of Yao's algorithm:" + aoi_tmc / times);
		System.out.println("The average AoI of greedy algorithm:" + aoi_greedy / times);
		System.out.println("The average AoI of Chiariotti's algorithm:" + aoi_tcom / times);
		System.out.println("The average AoI of our algorithm:" + aoi_0 / times);
		System.out.println();
	}
}
