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
#include "generate_exact_QFT_TFIM.cc"

//hallo


string random_Pauli_gate(std::uniform_real_distribution<double> & mydis);
string random_Pauli_pair(std::uniform_int_distribution<int> & mydistribution);

std::random_device device;
std::mt19937 mygenerator(device());
//mygenerator.discard(700000);

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


// for two  states of same size








// writes a  QFT circuit for N_qubits in the file qft.in
// skips phase gates where the phase applied is: e^(i*2*pi/2^{m+1}), e^(i*2*pi/2^{m+2})...e^(i*2*pi/2^{N})
void create_noisy_truncated_QFT_circuit(int N_qubits, int m, double r, int & cnt_rnd_gates, string & circuit_string, std::uniform_int_distribution<int> &  mypairdistribution){
    assert(m <= N_qubits);
    // erase everything  from circuit to begin with
    
	std::uniform_real_distribution<double> mydis(0.0, 1.0);
	double rnd = mydis(mygenerator);

    std::string circuit_file_name= "qft_noisy_truncated_m_"+std::to_string(m)+"_" + std::to_string(N_qubits) + ".in";
    std::ofstream my_circuit_stream;
    my_circuit_stream.open(circuit_file_name, std::ofstream::out | std::ofstream::trunc);

    string random_gate_1, random_gate_2;
    cnt_rnd_gates = 0;

    //cnt just labels the gates in their respective order
    int cnt=1;
    for (int i=N_qubits-1; i>=0; i--) {
        // every qubit gets a Hadamard gate
        my_circuit_stream << cnt << " " << 1 << " " << i <<" "<< "H"<<endl;
        circuit_string.append("H-"+to_string(i)+"-");
        rnd = mydis(mygenerator); 
        cnt++;
        if (rnd<r){
        	cnt_rnd_gates++;
        	random_gate_1 = random_Pauli_gate(mydis);
        	my_circuit_stream << cnt << " " << 1 << " " << i <<" "<< random_gate_1<<endl;
            circuit_string.append(random_gate_1+"-"+to_string(i)+"-");
        	cnt++;
        }

        for (int j=i-1; j>=0; j--) {
                // now we apply the controlled phase gates
            int k_phase = i-j+1;
            if(k_phase>m){
                break;
            }
            my_circuit_stream << cnt << " " << 2 << " " << j <<" "<< i << " " << "R"+std::to_string(k_phase) << endl;
            circuit_string.append("R"+std::to_string(k_phase)+"-"+to_string(j)+"-"+to_string(i)+"-");
            cnt++;
            rnd = mydis(mygenerator); 
            if (rnd<r){
            	cnt_rnd_gates++;
        	    random_gate_1 = random_Pauli_pair(mypairdistribution);   //random_Pauli_pair returns a string with 2 letters representing each random gate
        	    random_gate_2 = random_gate_1[1];
        	    random_gate_1 = random_gate_1[0];

        	    my_circuit_stream << cnt << " " << 1 << " " << j <<" "<<  random_gate_1 << endl;
                circuit_string.append(random_gate_1+"-"+to_string(j)+"-");
        	    cnt++;
        	    my_circuit_stream << cnt << " " << 1 << " " << i <<" "<<  random_gate_2 << endl;
                circuit_string.append(random_gate_2+"-"+to_string(i)+"-");
        	    cnt++;
            }
        }
    }
        // now swapping the qubits, the upper and lower indices refer to the way qubits are ordered
        // in the circuit, if we have a 3 qubit circuit diagram, then the uppermost qubit is number 0 and
        // the lowest one is number 2
        int upper = N_qubits-1; //N=3
        int lower = 0;
        if (N_qubits % 2 == 1) {
            //if we have an odd number of qubits we swap (N_qubits-1)/2 times
            for (int k=1; k<=((N_qubits-1)/2); k++) {
                my_circuit_stream << cnt << " " << 2 << " " << upper <<" "<< lower << " "<< "SWAP"<<endl;
                cnt++;
                rnd = mydis(mygenerator); 
                if (rnd<r){
                	cnt_rnd_gates++;
        	        random_gate_1 = random_Pauli_pair(mypairdistribution);
        	        random_gate_2 = random_gate_1[1];
        	        random_gate_1 = random_gate_1[0];

        	        my_circuit_stream << cnt << " " << 1 << " " << upper <<" "<<  random_gate_1 << endl;
                    circuit_string.append(random_gate_1+"-"+to_string(upper)+"-");
        	        cnt++;
        	        my_circuit_stream << cnt << " " << 1 << " " << lower <<" "<<  random_gate_2 << endl;
                    circuit_string.append(random_gate_2+"-"+to_string(lower)+"-");
        	        cnt++;
                }
                upper--;
                lower++;
            }
        }
        else {
            //when we have even number of qubits we swap (N_qubits)/2 times
            for (int k=1; k<=((N_qubits)/2); k++){
                my_circuit_stream << cnt << " " << 2 << " " << upper <<" "<< lower <<" "<<"SWAP"<< endl;
                cnt++;
                rnd = mydis(mygenerator); 
                if (rnd<r){
                	cnt_rnd_gates++;
        	        random_gate_1 = random_Pauli_pair(mypairdistribution);
        	        random_gate_2 = random_gate_1[1];
        	        random_gate_1 = random_gate_1[0];

        	        my_circuit_stream << cnt << " " << 1 << " " << upper <<" "<<  random_gate_1 << endl;
                    circuit_string.append(random_gate_1+"-"+to_string(upper)+"-");
        	        cnt++;
        	        my_circuit_stream << cnt << " " << 1 << " " << lower <<" "<<  random_gate_2 << endl;
                    circuit_string.append(random_gate_2+"-"+to_string(lower)+"-");
        	        cnt++;
                }
                upper--;
                lower++;
            }
        }
        my_circuit_stream.close();
}



void Apply_noisy_truncated_QFT(int N_qubits, double r, int iterations, psi_t input_state, double & mean, double & variance, double & std_dev, double & overlap_median, double & mean_of_log, double & log_variance, double & std_dev_log, double & stndrd_err_of_mean){
         // erase everything  from circuit to begin with
	qudos::GatesLibrary gl;
	double rnd;
	long double overlap;
	std::uniform_real_distribution<double> mydis(0.0, 1.0);

	int tmp_iter=0;
	string random_gate_string;
	vector<long double> overlap_vector(iterations, 0.0);
	vector<long double> log_overlap_vector(iterations, 0.0);
	vector<string> gate_vector(iterations);
    
    std::uniform_int_distribution<int> mypairdistribution(0,14);


    int m = 3; // truncation of QFT circuit
    write_truncated_QFT_circuit(N_qubits, m);

    // here we calculate the QFT from 
    std::string circuit_file_name= "qft_truncated_m_"+std::to_string(m)+"_" + std::to_string(N_qubits) + ".in";
    Circuit my_truncated_qft_circuit(circuit_file_name);
	
	psi_t qft_state = input_state;
	my_truncated_qft_circuit.Apply(qft_state);
    
    int count_random_gates = 0;

    std::string circuit_noisy_file_name = "qft_noisy_truncated_m_"+std::to_string(m)+"_" + std::to_string(N_qubits) + ".in";



    std::string circuit_string;
    for (int t=0; t<iterations; t++){
        circuit_string = "";
    	psi_t qft_noisy = input_state;
        count_random_gates = 0;
        create_noisy_truncated_QFT_circuit(N_qubits,  m,  r,  count_random_gates, circuit_string, mypairdistribution);
        // APPLY IT
        if (count_random_gates>0){
        	Circuit my_noisy_qft_circuit(circuit_noisy_file_name);
            my_noisy_qft_circuit.Apply(qft_noisy);
            overlap = get_overlap(qft_noisy, qft_state);
		    overlap_vector[tmp_iter] = overlap;
		    log_overlap_vector[tmp_iter] = std::log(overlap);
        }
        else{
            overlap_vector[tmp_iter] = 1.0;
            log_overlap_vector[tmp_iter] = 0.0;
        }

        tmp_iter++;
        gate_vector[t] = circuit_string;
    }


    // calculate mean and std dev and variance:
    
    for (int i=0; i<iterations; i++){
    	mean += overlap_vector[i];
    	mean_of_log += log_overlap_vector[i];
    }
    mean = mean/iterations;
    mean_of_log = mean_of_log/iterations;


    for (int i=0; i<iterations; i++){
    	variance = variance + (overlap_vector[i]-mean)*(overlap_vector[i]-mean);
        log_variance = log_variance + (log_overlap_vector[i]-mean_of_log)*(log_overlap_vector[i]-mean_of_log);
    }
    variance = variance/(iterations-1);
    log_variance = log_variance/(iterations-1);

    std_dev = sqrt(variance);
    stndrd_err_of_mean = std_dev/sqrt(iterations);
    std_dev_log = sqrt(log_variance);


    overlap_median = CalcMHWScore(overlap_vector);

	std::string aggregate_info = "aggregate_noisy_QFT_r"+to_string(r)+"_N_" + std::to_string(N_qubits) + ".in";
	std::ofstream my_circuit_stream(aggregate_info);
	int tmpcnt=1;
	for (int i=0; i<iterations; i++) {
            // every qubit gets a Hadamard gate
		my_circuit_stream << tmpcnt << "  " << std::setprecision(10) << overlap_vector[i] <<"  "<< gate_vector[i]<<" ";
		my_circuit_stream <<endl;
		tmpcnt++;
	}


}










// writes a  QFT circuit for N_qubits in the file qft.in
void Apply_noisy_Hadamard_circuit(int N_qubits, double r, int iterations, psi_t input_state, double & mean, double & variance, double & std_dev, double & overlap_median, double & mean_of_log, double & log_variance, double & std_dev_log, double & stndrd_err_of_mean){
        // erase everything  from circuit to begin with
	qudos::GatesLibrary gl;
	double rnd;
	long double overlap;
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
                random_gate_string = random_Pauli_gate(mydis);
				count_random_gates++;
				circuit_string.append(random_gate_string);
				qubit_vector[t].push_back(q);
				gate_cnt++;
			}
            
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
			log_overlap_vector[tmp_iter] = 0;
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
    log_variance = log_variance/(iterations-1);

    std_dev = sqrt(variance);
    stndrd_err_of_mean = std_dev/sqrt(iterations);
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



string random_Pauli_gate(std::uniform_real_distribution<double> & mydis){

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


string random_Pauli_pair(std::uniform_int_distribution<int> & mydistribution ){
	vector<std::string> random_gates = {"XX", "XY", "YX", "YY", "XZ", "ZX", "YZ", "ZY", "ZZ", "ZI", "IZ", "IY", "YI", "IX", "XI"};
    int number = mydistribution(mygenerator);
    std::string random_gate = random_gates[number];
	return random_gate;
}




int main (){
	int N_qubits = 25;//-----------              NOTE CHANGE THIS -----------------
	//string path = "TFIM_ground_states_LM";
	string relative_path = "TFIM_ground_states_2d_25qubits/TFIM_ground_state_2d_x_5_y_5_spins_gamma_3_total_spins_25.txt";
    //string relative_path = "TFIM_ground_states_LM/TFIM_ground_state_20spins.txt";

	psi_t psi_input(exp2(N_qubits));
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
		int iteration = (ceil(6.0/r))+1000;
		//Apply_noisy_Hadamard_circuit(N_qubits, r, iteration, psi_input, mean,variance, std_dev, overlap_median, mean_of_log, log_variance, std_dev_log, stndrd_err_of_mean);
		Apply_noisy_truncated_QFT(N_qubits, r, iteration, psi_input, mean, variance, std_dev, overlap_median, mean_of_log, log_variance, std_dev_log, stndrd_err_of_mean);
		out<<r<<" "<< iteration<< " "<< mean <<" "<< std_dev <<" " <<variance<<" " <<stndrd_err_of_mean <<" "<<overlap_median<< " " << mean_of_log <<" "<< log_variance<<" "<<std_dev_log<<" "<<endl;
	}
}
