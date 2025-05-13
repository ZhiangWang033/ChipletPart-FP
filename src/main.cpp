#include <iostream>
#include <vector>
#include <thread>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include "floorplan.h"
#include <omp.h>


void runSA(SACore* sa) {
  sa->run();
}

void parseFile(std::string filename, 
               float seperation,
               float reach,
               std::vector<Chiplet>& chiplets, 
               std::vector<BundledNet>& bundled_nets) {
  std::ifstream file_input(filename);
  if (!file_input.is_open()) {
    std::cout << "Error: cannot open file " << filename << std::endl;
    exit(1);
  }
  
  std::string cur_line;
  std::getline(file_input, cur_line);
  std::istringstream cur_line_buf(cur_line);
  std::vector<int> stats{std::istream_iterator<int>(cur_line_buf), 
                         std::istream_iterator<int>()};
  int num_chiplets = stats[1];
  int num_nets = stats[0];
  
  float halo_width = seperation; // This value is the seperation between two chiplets

  halo_width = 0.0;
  std::cout << "num_chiplets: " << num_chiplets << std::endl;
  std::cout << "num_nets: " << num_nets << std::endl;

  for (int i = 0; i < num_nets; i++) {
    std::getline(file_input, cur_line);
    std::istringstream cur_line_buf(cur_line);
    std::vector<float> hvec{std::istream_iterator<float>(cur_line_buf),
                            std::istream_iterator<float>()};
    int weight = hvec[0];
    //float reach = reach;
    int term_a = hvec[2] - 1;
    int term_b = hvec[3] - 1;
    bundled_nets.push_back(
        BundledNet(std::pair<int, int>(term_a, term_b), weight, reach));                
  }

  std::cout << "Finish parsing nets" << std::endl;

  float chiplet_area = 0.0;
  float chiplet_width = 0.0;
  float chiplet_height = 0.0;
  for (int i = 0; i < num_chiplets; i++) {
    //file_input >> chiplet_area;
    float x = 0.0;
    float y = 0.0;
    //float width = std::sqrt(chiplet_area);
    //float height = chiplet_area / width;
    //float min_area = chiplet_area;
    file_input >> chiplet_width;
    file_input >> chiplet_height;
    float width = chiplet_width;
    float height = chiplet_height;
    float min_area = width * height;
    chiplet_area = min_area;
    std::cout << "macro_id = " << i << " area = " << chiplet_area << "  "
              << "width = " << width << " height = " << height << std::endl;
    chiplets.push_back(Chiplet(x, y, width, height, min_area, halo_width));
  }

  std::cout << "check number of chiplets = " << chiplets.size() << std::endl;
  std::cout << "check number of nets = " << bundled_nets.size() << std::endl;
  for (auto& net : bundled_nets) {
    if (net.terminals.first < num_chiplets && net.terminals.second < num_chiplets) {
      continue;
    } else {
      std::cout << "Error: net terminals are out of range" << std::endl;
      std::cout << "terminals: " << net.terminals.first << " " << net.terminals.second << std::endl;
      exit(1);
    }
  }

  std::cout << "Finish parsing chiplets" << std::endl;

  file_input.close();  
}



float runReachSweep(std::string chiplet_file, float seperation, float reach, bool output_flag)
{
  std::vector<Chiplet> chiplets;
  std::vector<BundledNet> bundled_nets;
  parseFile(chiplet_file, seperation, reach, chiplets, bundled_nets);

  // define the parameters here
  float area_penalty_weight = 1000.0;
  float package_penalty_weight = 1.0;
  float net_penalty_weight = 100.0;

  float pos_swap_prob = 0.2;
  float neg_swap_prob = 0.2;
  float double_swap_prob = 0.2;
  float resize_prob = 0.2;
  float expand_prob = 0.2;

  int max_num_step = 2000;
  int num_perturb_per_step = 500;
  int num_threads = 10;
  unsigned init_seed = 0;
  int num_workers = 20;


  float max_temp = 100.0;  
  float min_temp = 1e-12;  // for sensity check only
  float cooling_rate = 0.95; // gemoetric cooling schedule
  int max_iter = 5000;  // number of iterations for a single SA worker
  int move_budget = 1000000; // total number of moves for each SA worker
  int num_move_per_iter = 200; // ( = max_budget / max_iter )

  // Setting for GWTW
  int gwtw_iter = 2;  
  float max_temp_derate_factor = 1.0;
  int top_k = 2; // each time copy the first top_k_ walkers
  std::vector<float> top_k_ratio { 0.5, 0.5}; // the copy ratio for each picked worker: sum(top_k_ratio_) = 1.0
  float sync_freq = 0.1; // sync up betwen workers for every max_iter * sync_freq iterations

  float max_cooling_rate = 0.99;
  float min_cooling_rate = 0.9;

  std::cout << "\n\n" << std::endl;
  std::cout << std::string(80, '*') << std::endl;

  std::cout << "[INFO] Reach = " << reach << std::endl;

  std::vector<std::unique_ptr<SACore> > workers;
  // try to workers with different cooling rate
  float delta_cooling_rate = (0.995 - cooling_rate) / (num_workers / 2 + 1);
  float best_cost = std::numeric_limits<float>::max();
  for (int worker_id = 0; worker_id < num_workers; worker_id++) {
    std::unique_ptr<SACore> sa = std::make_unique<SACore>(
      worker_id,
      chiplets,
      bundled_nets,
      area_penalty_weight,
      package_penalty_weight,
      net_penalty_weight,
      pos_swap_prob,
      neg_swap_prob,
      double_swap_prob,
      resize_prob,
      expand_prob,
      static_cast<int>(max_iter * sync_freq),
      num_perturb_per_step,
      cooling_rate + (worker_id - num_workers / 2) * delta_cooling_rate,
      init_seed + worker_id);
    workers.push_back(std::move(sa));
  }

  std::cout << "finish creating all the SA runs" << std::endl;

  workers[0]->initialize();
  float norm_area_penalty = workers[0]->getNormAreaPenalty();
  float norm_package_penalty = workers[0]->getNormPackagePenalty();
  float norm_net_penalty = workers[0]->getNormNetPenalty();

  std::cout << "norm_area_penalty: " << norm_area_penalty << std::endl;
  std::cout << "norm_package_penalty: " << norm_package_penalty << std::endl;
  std::cout << "norm_net_penalty: " << norm_net_penalty << std::endl;
 
  for (auto& worker : workers) {
    worker->setNormAreaPenalty(norm_area_penalty);
    worker->setNormPackagePenalty(norm_package_penalty);
    worker->setNormNetPenalty(norm_net_penalty);
  }


  omp_set_num_threads(std::max(8, num_threads));
  
  std::vector<int> copied_cnt;
  int copied_cnt_sum = 0;
  for (int i = 0; i < top_k; i++) {
    int cnt = static_cast<int>(std::floor(num_workers * top_k_ratio[i]));
    if (cnt == 0) {
      std::cout << "[INFO] The {} element of top_k_ratio is too samll (cnt = 0)!" << i << std::endl;
    }
    copied_cnt_sum += cnt;
    copied_cnt.push_back(cnt);
  }

  copied_cnt[0] += num_workers - copied_cnt_sum;
  std::string copied_vec_string("[INFO] The copied cnt for top vector = [ ");
  for (int i = 0; i < top_k; i++) {
    copied_vec_string += std::to_string(copied_cnt[i]) + " ";
  }
  copied_vec_string += "]";
  std::cout << copied_vec_string << std::endl;

  for (int gwtw_iter_id = 0; gwtw_iter_id < gwtw_iter; gwtw_iter_id++ ) {
    std::cout << std::string(80, '*') << std::endl;
    std::cout << "[INFO] GWTW Iteration " << gwtw_iter_id << " starts ..." << std::endl;
    if (gwtw_iter_id > 0) {
      max_temp *= max_temp_derate_factor;
    }

    std::cout << "[INFO] Set the max temperature to " << max_temp << std::endl;
    for (auto& worker : workers) {
      worker->setTemp(max_temp);
    }
    
    // in each GWTW iteration, there are multiple sync up
    int sync_iter = static_cast<int>(max_iter * sync_freq);
    for (int iter = 0; iter < max_iter; iter += sync_iter) {
      #pragma omp parallel for schedule(dynamic)
      for (int j = 0; j < workers.size(); j++) {
        workers[j]->run();
      }

      // Get the top K replicas based on HPWL
      std::sort(workers.begin(), workers.end(),
        [](const std::unique_ptr<SACore>& a, const std::unique_ptr<SACore>& b) {
          return a->getCost() < b->getCost();
        });

      std::cout << "********************** the current results after iteration " << iter << " *********************" << std::endl;
      for (auto& worker : workers) {
        std::cout << "[INFO] worker_id = " << worker->getWorkerId() << ", cost = " << worker->getCost() << std::endl;
      }


      int worker_cnt = top_k; 
      for (int j = 0; j < top_k; j++) {
        std::vector<Chiplet> best_chiplets;
        workers[j]->getMacros(best_chiplets);
        std::vector<int> pos_seq;
        workers[j]->getPosSeq(pos_seq);
        std::vector<int> neg_seq;
        workers[j]->getNegSeq(neg_seq);
        for (int i = 1; i < copied_cnt[j]; i++) {
          workers[worker_cnt]->setMacros(best_chiplets);
          workers[worker_cnt]->setPosSeq(pos_seq);
          workers[worker_cnt]->setNegSeq(neg_seq);
          worker_cnt++;      
        }      
      }
    }
  }

  SACore* best_sa = nullptr;
  for (auto& worker : workers) {
    if (worker->checkViolation() == false) {
      best_sa = worker.get();
    }
  }

  if (best_sa == nullptr) {
    std::cout << "Cannot find any valid solution" << std::endl;
    return 0.0;
  } else {
    std::cout << "Successfully find a valid solution" << std::endl;
  }


  std::vector<Chiplet> best_chiplets;
  best_sa->getMacros(best_chiplets);

  if (output_flag == true) {
    std::string output_file = chiplet_file + + ".reach_" + std::to_string(reach) + ".output";
    std::ofstream file_output(output_file);
    for (auto& chiplet : best_chiplets) {
      file_output << chiplet.getRealX()  << " " 
                  << chiplet.getRealY()  << " " 
                  << chiplet.getRealWidth() << " " 
                  << chiplet.getRealHeight() << std::endl;
    }
    file_output.close();
  }

  return best_sa->getPackageSize();
}


int main(int argc, char** argv) {
  std::string chiplet_file = std::string(argv[1]);
  float seperation = 1.0;
  std::vector<float> reaches;

  for (int i = 10; i <= 1000; i += 10) {
    reaches.push_back(i);
  }


  std::string reach_file = "reach_sweep.txt";
  std::ofstream file_output;
  file_output.open(reach_file, std::ios::out); // Fix: Change "w" to std::ios::out
  file_output.close();

  std::vector<float> package_areas;
  for (auto& reach : reaches) {
    float package_area = runReachSweep(chiplet_file, seperation, reach, false);
    package_areas.push_back(package_area);
    file_output.open(reach_file, std::ios::app); // Fix: Change "w" to std::ios::app
    file_output << reach << " " << package_area << std::endl;
    file_output.close();
  }

  return 0;
}


