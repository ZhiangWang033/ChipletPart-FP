#include <iostream>
#include <vector>
#include <thread>
#include <string>
#include <sstream>
#include <fstream>
#include <iterator>
#include "floorplan.h"


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
  float max_cooling_rate = 0.99;
  float min_cooling_rate = 0.9;

  std::cout << "\n\n" << std::endl;
  std::cout << std::string(80, '*') << std::endl;

  std::cout << "[INFO] Reach = " << reach << std::endl;

  // we first sweep the cooling rate 
  SACore* best_sa = nullptr;
  std::vector<SACore*> sa_containers;  // store all the SA runs to avoid memory leakage
  float best_cost = std::numeric_limits<float>::max();
  for (int i = 0; i < num_threads; i++) {
    float cooling_rate = max_cooling_rate - (max_cooling_rate - min_cooling_rate) * i / num_threads;     
    SACore* sa = new SACore(chiplets,
                            bundled_nets,
                            area_penalty_weight,
                            package_penalty_weight,
                            net_penalty_weight,
                            pos_swap_prob,
                            neg_swap_prob,
                            double_swap_prob,
                            resize_prob,
                            expand_prob,
                            max_num_step,
                            num_perturb_per_step,
                            cooling_rate,
                            init_seed);
    sa_containers.push_back(sa);
  }

  std::cout << "finish creating all the SA runs" << std::endl;

  sa_containers[0]->initialize();
  float norm_area_penalty = sa_containers[0]->getNormAreaPenalty();
  float norm_package_penalty = sa_containers[0]->getNormPackagePenalty();
  float norm_net_penalty = sa_containers[0]->getNormNetPenalty();

  std::cout << "norm_area_penalty: " << norm_area_penalty << std::endl;
  std::cout << "norm_package_penalty: " << norm_package_penalty << std::endl;
  std::cout << "norm_net_penalty: " << norm_net_penalty << std::endl;

  for (auto sa : sa_containers) {
    sa->setNormAreaPenalty(norm_area_penalty);
    sa->setNormPackagePenalty(norm_package_penalty);
    sa->setNormNetPenalty(norm_net_penalty);
  }

  std::vector<std::thread> threads;
  threads.reserve(sa_containers.size());
  for (auto& sa : sa_containers) {
    threads.emplace_back(runSA, sa);
  }
  
  for (auto& th : threads) {
    th.join();
  }

  int min_cost_sa_id = -1;
  int sa_id = 0;
  for (auto& sa : sa_containers) {
    if (sa->isValid() && sa->getCost() < best_cost) {
      best_cost = sa->getCost();
      best_sa = sa;
    }

    if (sa->getCost() < best_cost) {
      min_cost_sa_id = sa_id;
    }

    sa_id++;
  }

  if (best_sa == nullptr) {
    std::cout << "Cannot find any valid solution" << std::endl;
    best_sa = sa_containers[min_cost_sa_id];
  } else {
    std::cout << "Successfully find a valid solution" << std::endl;
  }
  
  std::vector<Chiplet> best_chiplets;
  best_sa->getMacros(best_chiplets);

  std::vector<int> pos_seq;
  best_sa->getPosSeq(pos_seq);
  std::vector<int> neg_seq;
  best_sa->getNegSeq(neg_seq);
  float best_cooling_rate = best_sa->getCoolingRate();
  best_sa = nullptr;

  // clean up 
  for (auto& sa : sa_containers) {
    delete sa;
  }

  sa_containers.clear();
   
  std::cout << "finish sweeping the cooling rate" << std::endl;
  
  best_cost = std::numeric_limits<float>::max();
  for (int i = 0; i < num_threads; i++) {   
    int seed = init_seed + i;
    SACore* sa = new SACore(chiplets,
                            bundled_nets,
                            area_penalty_weight,
                            package_penalty_weight,
                            net_penalty_weight,
                            pos_swap_prob,
                            neg_swap_prob,
                            double_swap_prob,
                            resize_prob,
                            expand_prob,
                            max_num_step,
                            num_perturb_per_step,
                            best_cooling_rate,
                            seed);
    sa_containers.push_back(sa);
  }

  std::cout << "finish initialize SA runs" << std::endl;


  //sa_containers[0]->initialize();
  //norm_area_penalty = sa_containers[0]->getNormAreaPenalty();
  //norm_package_penalty = sa_containers[0]->getNormPackagePenalty();
  //norm_net_penalty = sa_containers[0]->getNormNetPenalty();

  for (auto sa : sa_containers) {
    sa->setNormAreaPenalty(norm_area_penalty);
    sa->setNormPackagePenalty(norm_package_penalty);
    sa->setNormNetPenalty(norm_net_penalty);
    sa->setPosSeq(pos_seq);
    sa->setNegSeq(neg_seq);
  }

  std::cout << "finish setting the penalty parameters" << std::endl;

  threads.clear();
  threads.reserve(sa_containers.size());
  for (auto& sa : sa_containers) {
    threads.emplace_back(runSA, sa);
  }
  
  for (auto& th : threads) {
    th.join();
  }

  std::cout << "finish running SA" << std::endl;  

  min_cost_sa_id = -1;
  sa_id = 0;
  for (auto& sa : sa_containers) {
    if (sa->isValid() && sa->getCost() < best_cost) {
      best_cost = sa->getCost();
      best_sa = sa;
    }

    if (sa->getCost() < best_cost) {
      min_cost_sa_id = sa_id;
    }

    sa_id++;
  }

  if (best_sa == nullptr) {
    std::cout << "Cannot find any valid solution" << std::endl;
    best_sa = sa_containers[min_cost_sa_id];
    best_sa->checkViolation();
    return 0.0;
  } else {
    std::cout << "Successfully find a valid solution" << std::endl;
  }
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


