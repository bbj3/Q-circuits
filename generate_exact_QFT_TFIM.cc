#include "qudos.h"
using namespace std;
#include "gate.h"
#include "circuit.h"
using namespace qudos;
#include <iomanip>      // std::setprecision
#include <assert.h>


    //calculating the transformation of a state explicitly using equation 5.3
    //from Chuang and Nielsen
psi_t calculate_QFT_output_explicitly(psi_t psi){
    psi_t psi_transformed;
    int N = psi.size();
    psi_transformed.resize(N,0);

    const std::complex<double> I_(0.,1.);
    std::complex<double>  sum_for_new_coeff = 0;

    for (int k=0; k<N; k++) {
        sum_for_new_coeff = 0;
        for (int j=0; j<N; j++) {
            sum_for_new_coeff = sum_for_new_coeff + psi[j]*std::exp((2.*pi*I_*double(j)*double(k))/double(N));
        }
        psi_transformed[k] = sum_for_new_coeff/std::sqrt(N);
    }
    return psi_transformed;
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



// writes a  QFT circuit for N_qubits in the file qft.in
// skips phase gates where the phase applied is: e^(i*2*pi/2^{m+1}), e^(i*2*pi/2^{m+2})...e^(i*2*pi/2^{N})
void write_truncated_QFT_circuit(int N_qubits, int m){
    assert(m <= N_qubits);
    // erase everything  from circuit to begin with


    std::string circuit_file_name= "qft_truncated_m_"+std::to_string(m)+"_" + std::to_string(N_qubits) + ".in";
    std::ofstream my_circuit_stream (circuit_file_name);

        //cnt just labels the gates in their respective order
    int cnt=1;
    for (int i=N_qubits-1; i>=0; i--) {
            // every qubit gets a Hadamard gate
        my_circuit_stream << cnt << " " << 1 << " " << i <<" "<< "H"<<endl;

        cnt++;
        for (int j=i-1; j>=0; j--) {
                // now we apply the controlled phase gates
            int k_phase = i-j+1;
            if(k_phase>m){
                break;
            }
            my_circuit_stream << cnt << " " << 2 << " " << j <<" "<< i << " " << "R"+std::to_string(k_phase) << endl;
            cnt++;
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
                upper--;
                lower++;
            }
        }
        else {
            //when we have even number of qubits we swap (N_qubits)/2 times
            for (int k=1; k<=((N_qubits)/2); k++){
                my_circuit_stream << cnt << " " << 2 << " " << upper <<" "<< lower <<" "<<"SWAP"<< endl;
                cnt++;
                upper--;
                lower++;
            }
        }
        my_circuit_stream.close();
}


// for two  states of same size
  double get_overlap( psi_t & psi, psi_t & rho){
    assert(psi.size()==rho.size());
    std::complex<double> sum=0;
    double fidelity = 0;
    double psi_norm = 0;
    double rho_norm = 0;
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
void write_QFT_circuit(int N_qubits){
        // erase everything  from circuit to begin with


    std::string circuit_file_name= "qft_" + std::to_string(N_qubits) + ".in";
    std::ofstream my_circuit_stream (circuit_file_name);

        //cnt just labels the gates in their respective order
    int cnt=1;
    for (int i=N_qubits-1; i>=0; i--) {
            // every qubit gets a Hadamard gate
        my_circuit_stream << cnt << " " << 1 << " " << i <<" "<< "H"<<endl;

        cnt++;
        for (int j=i-1; j>=0; j--) {
                // now we apply the controlled phase gates
            int k_phase = i-j+1;

            my_circuit_stream << cnt << " " << 2 << " " << j <<" "<< i << " " << "R"+std::to_string(k_phase) << endl;
            cnt++;
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
                upper--;
                lower++;
            }
        }
        else {
            //when we have even number of qubits we swap (N_qubits)/2 times
            for (int k=1; k<=((N_qubits)/2); k++){
                my_circuit_stream << cnt << " " << 2 << " " << upper <<" "<< lower <<" "<<"SWAP"<< endl;
                cnt++;
                upper--;
                lower++;
            }
        }
        my_circuit_stream.close();
    }



// writes a  QFT circuit for N_qubits in the file qft.in
void write_Hadamard_circuit(int N_qubits){
        // erase everything  from circuit to begin with


    std::string circuit_file_name= "Hadamard_circuit_" + std::to_string(N_qubits) + ".in";
    std::ofstream my_circuit_stream (circuit_file_name);
    int cnt=1;
    for (int i=N_qubits-1; i>=0; i--) {
            // every qubit gets a Hadamard gate
        my_circuit_stream << cnt << " " << 1 << " " << i <<" "<< "H"<<endl;
        cnt++;
    }
}



// writes a  QFT circuit for N_qubits in the file qft.in
void write_Hadamard_gate(int qubit_index){
        // erase everything  from circuit to begin with


    std::string circuit_file_name= "Hadamard_gate_acts_on_qubit" + std::to_string(qubit_index) + ".in";
    std::ofstream my_circuit_stream (circuit_file_name);
    // every qubit gets a Hadamard gate
    my_circuit_stream << 1 << " " << 1 << " " << qubit_index <<" "<< "H"<<endl;
}

/*
int main (){
    int qubit = 25;
    int m =3;
    write_truncated_QFT_circuit(qubit, m);
        int spin_x = 5;
    int spin_y = 5;
    std::string spin_x_str = std::to_string(spin_x);
    std::string spin_y_str = std::to_string(spin_y);
    std::string N_string = std::to_string(spin_x*spin_y);
    int N = 20;
    //spin_x*spin_y; // number of qubits
    int Hilbertsize = exp2(N);
    psi_t psi_input(Hilbertsize, 0.0);

    //std::string circuit_name = "qft_"+std::to_string(N)+".in";
    write_QFT_circuit(N);

    string relative_path = "TFIM_ground_states_LM/TFIM_ground_state_20spins.txt";
    my_mutex.lock();
    get_state_from_txt_file_with_relative(relative_path,  psi_input);
    my_mutex.unlock();

    // here we calculate the QFT from 
    write_QFT_circuit(N);
    Circuit my_qft_circuit(circuit_name);

    for (auto element : psi_input){
        std::cout <<element<<std::endl;
    }


    psi_t correct_output_state_from_qft_circut = psi_input;
    my_qft_circuit.Apply(correct_output_state_from_qft_circut); // calculate QFT output using matrices

    string relative_path_for_QFT_output = "TFIM_ground_states_2d_QFT_output/TFIM_ground_states_2d_QFT_output_x_"+spin_x_str+"_y_"+spin_y_str+"_gamma3_total_spins_"+N_string+".txt";
    std::ofstream tmp_QFT_state_stream (relative_path_for_QFT_output);
    //tmp_QFT_state_stream <<"i "<<"Re " <<"Im "<<endl;
    for (int i=0; i<Hilbertsize; i++){
        tmp_QFT_state_stream <<std::setprecision(10) <<correct_output_state_from_qft_circut[i].real() <<" "<< correct_output_state_from_qft_circut[i].imag()<<endl;
    }
    tmp_QFT_state_stream.close();

   

    // here we calculate the Hadamard circuit output:


    psi_t correct_output_QFT = psi_input;
    std::string circuit_name = "qft_"+std::to_string(N)+".in";
    write_QFT_circuit(N);
    Circuit my_qft_circuit(circuit_name);
    my_qft_circuit.Apply(correct_output_QFT);
    
    std::vector<double> overlap_vector(N,0), m_vector(N,0);

    cout <<"N: "<<N<<endl;
    cout <<" m vector size: "<<m_vector.size()<<endl;

    psi_t QFT_trunc_state;
    double overlap;
    int mycnt = 0;
    for (int m=1; m<=N; m++){
        QFT_trunc_state = psi_input;

        string trunc_circuit_name = "qft_truncated_m_"+std::to_string(m)+"_" + std::to_string(N) + ".in";
        write_truncated_QFT_circuit(N, m);
        Circuit my_trunc_qft_circuit(trunc_circuit_name);
        my_trunc_qft_circuit.Apply(QFT_trunc_state);
        overlap = get_overlap(QFT_trunc_state, correct_output_QFT);
        m_vector[mycnt] = m;
        overlap_vector[mycnt] = overlap;
        mycnt++;
    }
    

    string path_m_overlap = "m_vs_fidelity.txt";
    std::ofstream m_overlap_stream(path_m_overlap);
    //tmp_QFT_state_stream <<"i "<<"Re " <<"Im "<<endl;
    for (int i=0; i<m_vector.size(); i++){
        m_overlap_stream <<m_vector[i]<<" "<< overlap_vector[i]<<endl;
        cout << "m: "<< m_vector[i]<< " overlap: "<<std::setprecision(10)<<overlap_vector[i]<<endl;
    }
    m_overlap_stream.close();



    string relative_path_for_Hadamard_output = "TFIM_ground_states_2d_Hadamard_output/TFIM_ground_state_hadamard_output_spin_per_dim_x_"+spin_x_str+"_y_"+spin_y_str+"_gamma3_total_spins_"+N_string+".txt";
    std::ofstream hadamard_state_stream (relative_path_for_Hadamard_output);
    //hadamard_state_stream <<"i "<<"Re " <<"Im "<<endl;
    for (int i=0; i<Hilbertsize; i++){
        hadamard_state_stream <<std::setprecision(10) <<correct_output_state_from_hadamard_circut[i].real() <<" "<<correct_output_state_from_hadamard_circut[i].imag()<<endl;
    }
    hadamard_state_stream.close();


    psi_t tmp_state = psi_input;
    

    for (int gatecnt=0; gatecnt<N; gatecnt++){
        my_hadamard_circuit.Apply_one_gate_at_a_time(tmp_state,  gatecnt);
        string relative_path_for_Hadamard_inter = "TFIM_ground_states_2d_Hadamard_output/H_intermediate_gate"+to_string(gatecnt)+"_spin_"+spin_x_str+"_y_"+spin_y_str+"_gamma3.txt";
        std::ofstream hadamard_interstate_stream(relative_path_for_Hadamard_inter);
        
    }



    // here we calculate the TRUNCATED - QFT from 
    int m = 2;
    string trunc_circuit_name = "qft_truncated_m_"+std::to_string(m)+"_" + std::to_string(N) + ".in";
    write_truncated_QFT_circuit(N, m);
    Circuit my_trunc_qft_circuit(trunc_circuit_name);

    for (auto element : psi_input){
        std::cout <<element<<std::endl;
    }


    psi_t correct_output_state_trunc_qft_circut = psi_input;
    my_trunc_qft_circuit.Apply(correct_output_state_trunc_qft_circut); // calculate QFT output using matrices

    string relative_path_for_QFT_trunc_output = "TFIM_ground_states_2d_QFT_output/TFIM_ground_2d_QFTtruncm_x_"+spin_x_str+"_y_"+spin_y_str+"_gamma3_.txt";
    std::ofstream tmp_QFTtrunc_state_stream (relative_path_for_QFT_trunc_output);
    //tmp_QFT_state_stream <<"i "<<"Re " <<"Im "<<endl;
    for (int i=0; i<Hilbertsize; i++){
        tmp_QFTtrunc_state_stream <<std::setprecision(10) <<correct_output_state_trunc_qft_circut[i].real() <<" "<< correct_output_state_trunc_qft_circut[i].imag()<<endl;
    }
    tmp_QFTtrunc_state_stream.close();



    psi_t experiment_state(4,0);
    experiment_state[1] = 1;

    psi_t tmp_experiment_0 = experiment_state;
    psi_t tmp_experiment_1 = experiment_state;
    
    experiment_state[1] = 1;
    tmp_experiment_0[1] = 1;


    int qubit_index= 0;
    std::string circuit_file_name= "Hadamard_gate_acts_on_qubit" + std::to_string(qubit_index) + ".in";
    write_Hadamard_gate(qubit_index);
    Circuit my_hadamard_gate_0(circuit_file_name);
    my_hadamard_gate_0.Apply(tmp_experiment_0);
    std::cout <<"hadamard gate acts on 0: " <<endl;
    for (int i=0; i<4; i++){
        std::cout <<tmp_experiment_0[i] <<" "<<endl;
    }
      
    qubit_index= 1;
    circuit_file_name= "Hadamard_gate_acts_on_qubit" + std::to_string(qubit_index) + ".in";
    write_Hadamard_gate(qubit_index);
    Circuit my_hadamard_gate_1(circuit_file_name);
    my_hadamard_gate_1.Apply(tmp_experiment_1);
    std::cout <<"hadamard gate acts on 1: " <<endl;
    for (int i=0; i<4; i++){
        std::cout <<tmp_experiment_1[i] <<" "<<endl;
    }



}
*/