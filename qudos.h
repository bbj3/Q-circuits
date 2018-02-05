#ifndef QUDOS_H_
#define QUDOS_H_

#include <vector>
#include <complex>
#include <mutex>  


namespace qudos {

typedef std::vector<std::vector<std::complex<double> > > mat_t;
typedef std::vector<std::complex<double> > psi_t;
typedef unsigned int conf_t;
typedef std::pair<int,std::complex<double> > connector_t;
const std::complex<double> I_(0.,1.);
const double pi(std::acos(-1.));
std::mutex my_mutex;
}

#include "gate.h"
#include "circuit.h"
#endif
