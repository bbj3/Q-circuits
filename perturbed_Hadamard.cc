#include "qudos.h"
using namespace std;
#include "gate.h"
#include "circuit.h"
using namespace qudos;
#include <iomanip>      // std::setprecision
#include <assert.h>
#include <random>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>



string random_Pauli_gate(std::mt19937 & mygenerator, std::uniform_real_distribution<double> & mydis);

void print_state(psi_t psi){
	int cnt=0;
	for(auto psival : psi) {
		std::cout<<cnt << ":  "<<psival<<std::endl;
		cnt++;
	}
}


psi_t get_equal_superposition_state(int Hilbertsize){
	psi_t psi;
	psi.resize(Hilbertsize,0);
	for (int i=0; i<Hilbertsize; i++) {
		psi[i] = 1/std::sqrt(double(Hilbertsize));
	}
	return psi;
}


double CalcMHWScore(vector<long double> scores)
{
  double median;
  size_t size = scores.size();

  sort(scores.begin(), scores.end());

  if (size  % 2 == 0)
  {
      median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
  }
  else 
  {
      median = scores[size / 2];
  }

  return median;
}




void get_state_from_txt_file_with_relative(std::string relative_path, psi_t & mystate){

	std::ifstream myfile (relative_path);
	std::string line;
    // run through all the gates
	int cnt = 0;
	std::complex<double> element;
	while (std::getline(myfile, line))
	{   
      // get info on gates in particular line
		std::istringstream iss(line);
		iss >> element;
		mystate[cnt] = element;
		cnt++;
	}
	myfile.close();
}


// for two  states of same size
long double get_overlap( psi_t & psi, psi_t & rho){
	assert(psi.size()==rho.size());
	std::complex<long double> sum=0;
	long double  fidelity = 0;
	long double  psi_norm = 0;
	long double  rho_norm = 0;
	for (int S=0; S<psi.size(); S++){
		sum+=psi[S]*std::conj(rho[S]);
		rho_norm += abs(rho[S])*abs(rho[S]);
		psi_norm += abs(psi[S])*abs(psi[S]);
	} 
    fidelity = std::real(sum*std::conj(sum)); // just taking real beacause it's a real number and want to repr. it as a double
    fidelity = std::sqrt(fidelity)/std::sqrt(psi_norm*rho_norm);
    return fidelity;
} 




// writes a  QFT circuit for N_qubits in the file qft.in
void Apply_noisy_Hadamard_circuit(int N_qubits, double r, int iterations, psi_t input_state, double & mean, double & variance, double & std_dev, double & overlap_median, double & mean_of_log, double & log_variance, double & std_dev_log, double & stndrd_err_of_mean){
        // erase everything  from circuit to begin with
	qudos::GatesLibrary gl;
	double rnd;
	long double overlap;
	std::mt19937 mygenerator (time(NULL));
	std::uniform_real_distribution<double> mydis(0.0, 1.0);

	int tmp_iter=0;
	string random_gate_string;
	vector<long double> overlap_vector(iterations, 0.0);
	vector<long double> log_overlap_vector(iterations, 0.0);
	vector<string> gate_vector(iterations);
	vector<vector<int> > qubit_vector(iterations, vector<int>(0));  // here we record which qubit we apply the gates to in each iteration (circuit)


	Circuit pure_H_circuit;
	vector<int> qubits(1);

	for (int q=N_qubits-1; q>=0; q--) {
		random_gate_string = "H";
		qubits[0] = q;
		pure_H_circuit.Push(Gate(gl(random_gate_string), qubits, random_gate_string));
	}
	psi_t H_state = input_state;
	pure_H_circuit.Apply(H_state);
	vector<char> qubit_vec(iterations);
    
    int count_random_gates = 0;

	for (int t=0; t<iterations; t++){
        count_random_gates = 0;
		psi_t output_state = input_state;
		Circuit perturbed_H_circuit;
		string circuit_string = "";

		int gate_cnt = 0;
		for (int q=N_qubits-1; q>=0; q--){
			rnd = mydis(mygenerator);
			cout <<" rnd: " <<rnd <<endl;
			random_gate_string = "H";
			circuit_string.append(random_gate_string);
			qubit_vector[t].push_back(q);
			gate_cnt++;
			if (rnd<r){
                random_gate_string = random_Pauli_gate(mygenerator, mydis);
				count_random_gates++;
				circuit_string.append(random_gate_string);
				qubit_vector[t].push_back(q);
				gate_cnt++;
			}
            // every qubit gets a Hadamard gate
			//perturbed_H_circuit.Push(Gate(gl(random_gate_string), qubits, random_gate_string));
		}
		gate_vector[t] = circuit_string;

		if (count_random_gates>0){  
		    // if there was a random pauli gate  then we create the circuit and apply and calculate fidelity
			vector<int> qubits(1);
			int q =0;
            qubits[0] = qubit_vector[t][q];
			for (char &c : gate_vector[t]){
				random_gate_string = c;
				perturbed_H_circuit.Push(Gate(gl(random_gate_string), qubits, random_gate_string));
				std::cout <<c<<" " << qubit_vector[t][q]<<" ";
				q++;
				qubits[0] = qubit_vector[t][q];
			}
			std::cout <<" "<<std::endl;
		    perturbed_H_circuit.Apply(output_state);
		    //print_state(output_state);
		    overlap = get_overlap(output_state, H_state);
		    overlap_vector[tmp_iter] = overlap;
		    log_overlap_vector[tmp_iter] = std::log(overlap);
		}
		else{
			overlap_vector[tmp_iter] = 1.0;
		}

		tmp_iter++;
	}
    

    // calculate mean and std dev and variance:
    
    for (int i=0; i<iterations; i++){
    	mean += overlap_vector[i];
    	mean_of_log += log_overlap_vector[i];
    	cout << mean << " "<< overlap_vector[i]<<"  " <<log_overlap_vector[i]<<" "<<endl;
    }
    mean = mean/iterations;
    mean_of_log = mean_of_log/iterations;


    for (int i=0; i<iterations; i++){
    	variance = variance + (overlap_vector[i]-mean)*(overlap_vector[i]-mean);
        log_variance = log_variance + (log_overlap_vector[i]-mean_of_log)*(log_overlap_vector[i]-mean_of_log);
    }
    variance = variance/(iterations-1);
    stndrd_err_of_mean = variance/sqrt(iterations);
    log_variance = log_variance/(iterations-1);

    std_dev = sqrt(variance);
    std_dev_log = sqrt(log_variance);


    overlap_median = CalcMHWScore(overlap_vector);

	std::string circuit_file_name= "perturberd_Hadamard_r"+to_string(r)+"_N_" + std::to_string(N_qubits) + ".in";
	std::ofstream my_circuit_stream(circuit_file_name);
	int tmpcnt=1;
	for (int i=0; i<iterations; i++) {
            // every qubit gets a Hadamard gate
		my_circuit_stream << tmpcnt << "  " << std::setprecision(10) << overlap_vector[i] <<"  "<< gate_vector[i]<<" ";
		for (int j=0; j<qubit_vector[i].size(); j++){
			my_circuit_stream << qubit_vector[i][j]<< " ";
		}
		my_circuit_stream <<endl;
		tmpcnt++;
	}
}



string random_Pauli_gate(std::mt19937 & mygenerator, std::uniform_real_distribution<double> & mydis){

	double rnd = mydis(mygenerator);
	double equal_prob = 1.0/3.0;
	if (rnd<equal_prob){
		return "X";
	}
	else if ( equal_prob<=rnd and rnd<2*equal_prob){
		return "Y";
	}
	else {
		return "Z";
	}
}




int main (){
	int N_qubits = 20;//-----------              NOTE CHANGE THIS -----------------
	string path = "TFIM_ground_states_LM";
	string relative_path = "TFIM_ground_states_LM/TFIM_ground_state_20spins.txt";

	psi_t psi_input(exp2(20));
	cout << " "<<relative_path<<endl;
	int Hilbertsize = exp2(N_qubits);
	psi_input = get_equal_superposition_state(Hilbertsize);
	get_state_from_txt_file_with_relative(relative_path,  psi_input);
	for (int i=0; i<20; i++){
		cout << psi_input[i]<<endl;
	}
    
    double mean,variance, std_dev, overlap_median, mean_of_log, std_dev_log, log_variance, stndrd_err_of_mean;

     std::ofstream out;

    // std::ios::app is the open mode "append" meaning
    // new data will be written to the end of the file.
    out.open("stats_noisy_Hada_5_feb"+std::to_string(N_qubits)+".txt", std::ios::app);

    
    out <<"r " <<"iterations "<<"mean "<<"std_dev "<<"variance "<<" standard_error_of_the_mean "<< "overlap_median " <<" mean_of_log" << "variance of log "<< "std_dev_log"<<endl;
	for (double r=0.001; r<=0.6; r=r*2.0){
		mean=0;
		variance=0;
		std_dev=0;
		mean_of_log =0;
		std_dev_log = 0;
		log_variance = 0;
        stndrd_err_of_mean = 0;
		int iteration = int (ceil(6.0/r))+1000;
		Apply_noisy_Hadamard_circuit(N_qubits, r, iteration, psi_input, mean,variance, std_dev, overlap_median, mean_of_log, log_variance, std_dev_log, stndrd_err_of_mean);
		out<<r<<" "<< iteration<< " "<< mean <<" "<< std_dev <<" " <<variance<<" " <<stndrd_err_of_mean <<" "<<overlap_median<< " " << mean_of_log <<" "<< log_variance<<" "<<std_dev_log<<" "<<endl;
	}
}
