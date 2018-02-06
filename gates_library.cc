#include <string>
#include <map>
#include <complex>
#include <cmath>
#include <iostream>

namespace qudos{

  class GatesLibrary{

    std::map<std::string,qudos::mat_t> gateslib_;

  public:

    explicit GatesLibrary(){
      InitDefaultGates();
    }

    void InitDefaultGates(){

      const std::complex<double> I_(0.,1.);
      const double pi(std::acos(-1.));

      qudos::mat_t Oneq,Twoq;

      Oneq.resize(2,std::vector<std::complex<double> > (2,0.));
      Oneq[0][1]=1.0;
      Oneq[1][0]=1.0;
      gateslib_["sx"]=Oneq;
      gateslib_["NOT"]=Oneq;
      gateslib_["X"]=Oneq;

      SetMatZero(Oneq);
      Oneq[0][0]=1.0;
      Oneq[1][1]=-1.0;
      gateslib_["sz"]=Oneq;
      gateslib_["Z"]=Oneq;

      SetMatZero(Oneq);
      Oneq[0][1]=-I_;
      Oneq[1][0]=I_;
      gateslib_["sy"]=Oneq;
      gateslib_["Y"]=Oneq;

      SetMatZero(Oneq);
      Oneq[0][0]=1.0;
      Oneq[1][1]=std::exp(I_*pi/4.);
      gateslib_["T"]=Oneq;

      SetMatZero(Oneq);
      Oneq[0][0]=1.0;
      Oneq[1][1]=1.0;
      gateslib_["ID1"]=Oneq;


      SetMatZero(Oneq);
      Oneq[0][0]=1.0;
      Oneq[1][1]=1.0;
      gateslib_["I"]=Oneq;

      SetMatZero(Oneq);
      Oneq[0][0]=1.0;
      Oneq[1][1]=-1.0;
      Oneq[0][1]=1.0;
      Oneq[1][0]=1.0;
      MatMultiply(Oneq,1./std::sqrt(2.));
      gateslib_["H"]=Oneq;

      SetMatZero(Oneq);
      Oneq[0][0]=1.0;
      Oneq[1][1]=1.0;
      Oneq[0][1]=-I_;
      Oneq[1][0]=-I_;
      MatMultiply(Oneq,(1.+I_)/2.);
      gateslib_["XHALF"]=Oneq;
      gateslib_["NOTHALF"]=Oneq;

      SetMatZero(Oneq);
      Oneq[0][0]=1.0;
      Oneq[1][1]=1.0;
      Oneq[0][1]=-1;
      Oneq[1][0]=1;
      MatMultiply(Oneq,(1.+I_)/2.);
      gateslib_["YHALF"]=Oneq;

      Twoq.resize(4,std::vector<std::complex<double> > (4,0.));
      Twoq[0][0]=1.0;
      Twoq[1][1]=-1.0;
      Twoq[2][2]=-1.0;
      Twoq[3][3]=1.0;
      gateslib_["sz*sz"]=Twoq;

      SetMatZero(Twoq);
      Twoq[0][0]=1.0;
      Twoq[1][2]=1.0;
      Twoq[2][1]=1.0;
      Twoq[3][3]=1.0;
      gateslib_["SWAP"]=Twoq;

      SetMatZero(Twoq);
      Twoq[0][0]=1.0;
      Twoq[1][1]=1.0;
      Twoq[2][2]=1.0;
      Twoq[3][3]=-1.0;
      gateslib_["CZ"]=Twoq;

      SetMatZero(Twoq);
      Twoq[0][0]=1.0;
      Twoq[1][1]=1.0;
      Twoq[2][3]=1.0;
      Twoq[3][2]=1.0;
      gateslib_["CX"]=Twoq;
      gateslib_["CNOT"]=Twoq;

      SetMatZero(Twoq);
      Twoq[0][0]=1.0;
      Twoq[1][1]=1.0;
      Twoq[2][3]=-I_;
      Twoq[3][2]=I_;
      gateslib_["CY"]=Twoq;

      for(int k=0;k<30;k++){
        SetMatZero(Twoq);
        Twoq[0][0]=1.0;
        Twoq[1][1]=1.0;
        Twoq[2][2]=1.0;
        Twoq[3][3]=std::exp(2.*pi*I_/double(1<<k));
        gateslib_["R"+std::to_string(k)]=Twoq;
      }

    }

    qudos::mat_t operator()(std::string gatename){
      auto gate=gateslib_.find(gatename);

      if(gate!=gateslib_.end()){
        return gate->second;
      }
      else{
        std::string delimiter = "*";

        qudos::mat_t gatet=gateslib_["ID1"];
        qudos::mat_t gatett=gatet;

        if(gatename.find(delimiter) != std::string::npos){

          std::string s=gatename;

          size_t pos = 0;
          std::string token;

          while ((pos = s.find(delimiter)) != std::string::npos) {
            token = s.substr(0, pos);

            auto gate=gateslib_.find(token);

            if(gate!=gateslib_.end()){
              gatett=MatMultiply(gatet,gate->second);
              gatet=gatett;
              s.erase(0, pos + delimiter.length());
            }
            else{
              std::cerr<<"Gate "<<gatename<<" not found in the library"<<std::endl;
              std::abort();
            }

          }

          auto gate=gateslib_.find(s);

          if(gate!=gateslib_.end()){
            gatett=MatMultiply(gatet,gate->second);
            gatet=gatett;
            s.erase(0, pos + delimiter.length());
          }
          else{
            std::cerr<<"Gate "<<gatename<<" not found in the library"<<std::endl;
            std::abort();
          }

          return gatet;
        }

        std::cerr<<"Gate "<<gatename<<" not found in the library"<<std::endl;
        std::abort();
      }
    }

    void SetMatZero(qudos::mat_t & Mat){
      for(int i=0;i<Mat.size();i++){
        for(int j=0;j<Mat[i].size();j++){
          Mat[i][j]=0.;
        }
      }
    }

    void MatMultiply(qudos::mat_t & Mat,std::complex<double> number){
      for(int i=0;i<Mat.size();i++){
        for(int j=0;j<Mat[i].size();j++){
          Mat[i][j]*=number;
        }
      }
    }

    void MatTranspose(qudos::mat_t & Mat){
      mat_t Mat1=Mat;

      for(int i=0;i<Mat.size();i++){
        for(int j=0;j<Mat[i].size();j++){
          Mat1[i][j]=Mat[j][i];
        }
      }
      Mat=Mat1;
    }

    qudos::mat_t MatMultiply(const qudos::mat_t & Mat1,const qudos::mat_t & Mat2){
      mat_t Matr(Mat1.size(),std::vector<std::complex<double> > (Mat2[0].size()));


      for(int i=0;i<Mat1.size();i++){

        for(int j=0;j<Mat1[i].size();j++){
          Matr[i][j]=0.;

          for(int k=0;k<Mat1[i].size();k++){
            Matr[i][j]+=Mat1[i][k]*Mat2[k][j];
          }

        }
      }
      return Matr;
    }

    void PrintGate(const qudos::mat_t & Mat){
      for(int i=0;i<Mat.size();i++){
        for(int j=0;j<Mat[i].size();j++){
          std::cout<<Mat[i][j]<<" ";
        }
        std::cout<<std::endl;
      }
    }

  };


}
