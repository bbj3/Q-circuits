#include <iostream>
#include <vector>
#include <complex>
#include <limits>

namespace qudos{
  using namespace std;

  class Gate{

    //Diagonal part of the gate
    vector<complex<double> > diag_;

    //Off-Diagonal part
    vector<vector<qudos::connector_t> > offdiag_;

    //List of qubits on which the gate acts
    const vector<int> qubits_;

    vector<int> maskloc_;
    vector<conf_t> bitmask_;
    conf_t masksites_; //is 1 on bits on which the gate operates

    bool has_diagonal_;
    bool diagonal_identity;
    bool has_offdiagonal_;

    string name_;

  public:
    explicit Gate(const qudos::mat_t & gate,const vector<int> & qubits,string name=""):qubits_(qubits),name_(name){
      Init(gate);
    }

    void Init(const qudos::mat_t & gate){

      const double epsilon=std::numeric_limits<double>::epsilon();
      // const double epsilon=1.0e-10;

      diag_.reserve(gate.size());

      has_diagonal_=false;
      for(int i=0;i<gate.size();i++){
        diag_.push_back(gate[i][i]);
        if(std::abs(gate[i][i])>epsilon){
          has_diagonal_=true;
        }
      }

      diagonal_identity=false;
      if(has_diagonal_){
        diagonal_identity=true;
        for(int i=1;i<diag_.size();i++){
          if(std::abs(diag_[i]-diag_[0])>epsilon){
            diagonal_identity=false;
            break;
          }
        }
      }

      offdiag_.resize(gate.size());

      has_offdiagonal_=false;
      for(int i=0;i<gate.size();i++){
        for(int j=0;j<gate[i].size();j++){

          if(i!=j && std::abs(gate[i][j])>epsilon){
            const qudos::connector_t connector(j,gate[i][j]);
            offdiag_[i].push_back(connector);
            has_offdiagonal_=true;
          }

        }
      }


      maskloc_.resize(qubits_.size());
      masksites_=0;
      for(int loc=0;loc<qubits_.size();loc++){
        maskloc_[loc]=(1<<qubits_[loc]);
        SetBit(masksites_,qubits_[loc]);
      }

      const int localhilbsize=1<<qubits_.size();

      bitmask_.resize(localhilbsize,0);

      for(int lc=0;lc<localhilbsize;lc++){
        for(int loc=0;loc<qubits_.size();loc++){
          if(GetBit(lc,loc)){
            SetBit(bitmask_[lc],qubits_[loc]);
          }
        }
      }

    }

    std::string get_name(){
      return name_;
    }

    void Apply(const psi_t & psi, psi_t & psiout)const{
      psi_t::const_iterator psi_el=psi.begin();
      psi_t::iterator psiout_el=psiout.begin();

      if(has_diagonal_){
        if(diagonal_identity){
          #pragma omp parallel for
          for(conf_t conf=0;conf<psi.size();conf++){

            //Diagonal part of the qubit
            // *psiout_el=diag_[0]*(*psi_el);

            // psi_el++;
            // psiout_el++;
            psiout[conf]=diag_[0]*psi[conf];

          }
        }
        else{
          #pragma omp parallel for
          for(conf_t conf=0;conf<psi.size();conf++){
            const conf_t locqubit=LocalConf(conf);

            //Diagonal part of the qubit
            // *psiout_el=diag_[locqubit]*(*psi_el);
            //
            // psi_el++;
            // psiout_el++;

            psiout[conf]=diag_[locqubit]*psi[conf];
          }
        }
      }

      if(has_offdiagonal_){

        if(!has_diagonal_){
          #pragma omp parallel for
          for(int i=0;i<psiout.size();i++){
            psiout[i]=0;
          }
        }

        // psiout_el=psiout.begin();

        #pragma omp parallel for
        for(conf_t conf=0;conf<psi.size();conf++){
          const conf_t locqubit=LocalConf(conf);

          //Offdiagonal part, looping over all non-zero matrix elements
          for(const connector_t & c : offdiag_[locqubit]){
            const std::complex<double> mel=c.second;
            // *psiout_el+=psi[GlobalConf(conf,c.first)]*mel;
            psiout[conf]+=psi[GlobalConf(conf,c.first)]*mel;
          }

          // psiout_el++;
        }
      }



    }

    inline int LocalConf(conf_t conf)const{

      int state=0;
      int loc=0;
      for(const int & mask : maskloc_){
        if(conf & mask){
          state|=(1<<loc);
        }
        loc++;
      }

      return state;
    }

    //Returns the global configuration where the spins in the locl configuration lc
    //have been set globally
    inline conf_t GlobalConf(conf_t conf,int lc)const{
      return conf^((conf^bitmask_[lc])&masksites_);
    }

    std::vector<int> Qubits()const{
      return qubits_;
    }

    inline void SetBit(conf_t & c,conf_t bit)const{
      c|=(1<<bit);
    }
    inline conf_t GetBit(conf_t c,conf_t bit)const{
      return (c>>bit)&1;
    }

    inline string Name()const{
      return name_;
    }

  };

}
