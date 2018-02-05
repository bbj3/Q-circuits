#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
//#include <mpi.h>

namespace qudos{
	using namespace std;

	class Circuit{

		vector<qudos::Gate> gates_;
		int nqubits_;
		int mynode_;

	public:

		explicit Circuit(){
			// MPI_Comm_rank(MPI_COMM_WORLD, &mynode_);
			nqubits_=0;
		}

		explicit Circuit(const qudos::Gate & gate){
			// MPI_Comm_rank(MPI_COMM_WORLD, &mynode_);
			nqubits_=0;
			Push(gate);
		}

		explicit Circuit(const vector<qudos::Gate> & gates){
			// MPI_Comm_rank(MPI_COMM_WORLD, &mynode_);
			nqubits_=0;
			for(const auto & gate : gates){
				Push(gate);
			}
		}
		explicit Circuit(string filename){
			// MPI_Comm_rank(MPI_COMM_WORLD, &mynode_);
			nqubits_=0;
			ReadFromFile(filename);
		}
        /*
		void pop_gate(){
			gates_.erase(gates_.begin());
		}
        */


		void ReadFromFile(string filename){
			qudos::GatesLibrary gl;

			ifstream filec;
			filec.open(filename.c_str());
			if(!filec.good()&&mynode_==0){
				std::cerr<<"#Cannot open file"<<filename<<std::endl;
				std::abort();
			}

			while(filec.good()){
				int cycle;
				filec>>cycle;
			
				int qubit_order;
				filec>>qubit_order;

				vector<int> qubits(qubit_order);
				for(auto & qubit : qubits){
					filec>>qubit;
				}

				string qubitname;
				filec>>qubitname;
				if(!filec.eof()){
					Push(qudos::Gate(gl(qubitname),qubits,qubitname));
				}
			}

			if(mynode_==0){
				std::cout<<"#Circuit loaded from "<<filename<<std::endl;
				std::cout<<"#Number of qubits : "<<nqubits_<<std::endl;
				std::cout<<"#Number of gates : "<<Ngates()<<std::endl;
			}
		}

		void Push(const qudos::Gate & gate){
			gates_.push_back(gate);

			vector<int> qubits=gate.Qubits();
			const int maxqubit=*(std::max_element(qubits.begin(),qubits.end()));

			if(maxqubit+1>nqubits_){
				nqubits_=maxqubit+1;
			}
		}

		void Apply(psi_t & psi)const{
			psi_t psiout(psi);

			for(const auto & gate : gates_){
				gate.Apply(psi,psiout);
				psi.swap(psiout);
			}
		}
        
        void Apply_one_gate_at_a_time(psi_t & psi, int gate_number)const{
			psi_t psiout(psi);
            qudos::Gate gate = gates_[gate_number];
            cout<<"gate name: " << gate.get_name() <<endl;
			gate.Apply(psi,psiout);
			psi.swap(psiout);
			//pop_gate();
		}



		inline int Ngates()const{
			return gates_.size();
		}

		inline int Nqubits()const{
			return nqubits_;
		}

		inline int Hilbsize()const{
			return 1<<Nqubits();
		}

		inline void InitState(psi_t & psi)const{
			psi.resize(Hilbsize(),1.0/std::sqrt(Hilbsize()));
		}

	};

}
