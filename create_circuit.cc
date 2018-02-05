#include "qudos.h"
using namespace std;
#include "gate.h"
#include "circuit.h"
using namespace qudos;
#include <assert.h>



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




int main(){
    //*
    int N = 16;
    int m = 5;
    write_truncated_QFT_circuit(N, m);
    //std::string circuit_name = "my_circuit_name";

    write_QFT_circuit(N);
    //Circuit my_qft_circuit(circuit_name);


}