

class SGD{

	int n_size;
	double scalar_adadelta_E_g2,  scalar_adadelta_E_delta_x2,  adadelta_rho,  adadelta_eps, delta_x, delta_x_2;

	vector<double> adadelta_E_g2, adadelta_E_delta_x2;
	vector<double> adadelta_delta_x;


        //adagrad 
	double adagrad_eps;
	vector<double> grad_squared;
	double adagrad_gamma;


        //Adam and adamax:
	double alpha, beta_1, beta_2,adam_eps;
	vector<double> m, v, m_hat, v_hat;
	int time;

        //RMSProp
	vector<double> rms_E_g2;

	double infinity_norm=0;


        // NOTE: only use one method at a time otherwise everything might go to shit
	    // that is if you want to use two different SGD methods then initialize two different objects
public:
       // good example params are:  double my_Adagrad_eps=1.0e-10,  my_adagrad_gamma=0.005,  my_adam_alpha=0.001,  my_adam_beta_1=0.9,  my_adam_beta_2=0.999,  my_adam_eps=1.0e-8; my_adadelta_rho = 0.97  my_adadelta_eps = 10^-9
	   // methods that have worked well: step_gradient_element_wise_adadelta and step_gradient_Adam
	SGD(int my_n_size, double my_adadelta_eps, double my_adadelta_rho, double my_Adagrad_eps=0, double my_adagrad_gamma=0, double my_adam_alpha=0, double my_adam_beta_1=0, double my_adam_beta_2=0, double my_adam_eps=0){
		n_size = my_n_size;
		adadelta_eps = my_adadelta_eps;
		adadelta_rho = my_adadelta_rho;
		scalar_adadelta_E_g2=0, scalar_adadelta_E_delta_x2=0, delta_x=0, delta_x_2=0;
		adadelta_E_g2.resize(n_size, 0.0), adadelta_E_delta_x2.resize(n_size, 0.0);
		adadelta_delta_x.resize(n_size, 0.0);
		delta_x = 0;
		delta_x_2 = 0;
            // for adagrad
		grad_squared.resize(n_size, 0.0);
		adagrad_eps = my_Adagrad_eps;
		adagrad_gamma = my_adagrad_gamma;
		m.resize(n_size, 0.0);
		v.resize(n_size, 0.0);
		alpha = my_adam_alpha;
		beta_1 = my_adam_beta_1;
		beta_2 = my_adam_beta_2;
		m_hat.resize(n_size,0.0);
		v_hat.resize(n_size,0.0);
		time = 0;
		adam_eps = my_adam_eps;

		rms_E_g2.resize(n_size,0.0);

	}


	void step_gradient_exact_adadelta(vector<double> & x,vector<double> & grad){
		double g_2 = 0;
		vector<double> adadelta_delta_x(n_size);
		for(int k=0; k<grad.size(); k++){
			g_2 += grad[k]*grad[k];
		}    
		scalar_adadelta_E_g2 = adadelta_rho*scalar_adadelta_E_g2+(1-adadelta_rho)*g_2;
		for (int i=0; i<n_size; i++){
			adadelta_delta_x[i] = -(std::sqrt(scalar_adadelta_E_delta_x2+adadelta_eps)/std::sqrt(scalar_adadelta_E_g2+adadelta_eps))*grad[i];
		}    
		double delta_x_2 = 0;
		for (int i=0; i<n_size; i++){
			delta_x_2 += adadelta_delta_x[i]*adadelta_delta_x[i];
		}
		scalar_adadelta_E_delta_x2 = adadelta_rho*scalar_adadelta_E_delta_x2+(1-adadelta_rho)*delta_x_2;  
		for (int i=0; i<n_size; i++){
			x[i] = x[i]+adadelta_delta_x[i];
		}
	}


	void step_gradient_element_wise_adadelta(vector<double> & x,vector<double> & grad){

		for (int i=0; i<x.size(); i++){
			adadelta_E_g2[i] = adadelta_rho*adadelta_E_g2[i]+(1-adadelta_rho)*grad[i]*grad[i];
		}
		for (int i=0; i<x.size(); i++){
			adadelta_delta_x[i] = -(std::sqrt((adadelta_E_delta_x2[i]*adadelta_E_delta_x2[i]+adadelta_eps)/(adadelta_E_g2[i]*adadelta_E_g2[i]+adadelta_eps)))*grad[i];
		}
		for (int i=0; i<x.size(); i++){
			adadelta_E_delta_x2[i] = adadelta_rho*adadelta_E_delta_x2[i]+(1-adadelta_rho)*adadelta_delta_x[i]*adadelta_delta_x[i];
		}
		for (int i=0; i<x.size(); i++){
			x[i] = x[i]+adadelta_delta_x[i];
		}
	}


	void step_gradient_adagrad(vector<double> & x, vector<double> & grad){
			//update-a grad_squared
		for (int i=0; i<n_size; i++){
			grad_squared[i] += grad[i]*grad[i];
		}

		for (int i=0; i<n_size; i++){
			x[i] = x[i]-(adagrad_gamma/(std::sqrt(grad_squared[i])+adagrad_eps))*grad[i];
		}

	}


	void step_gradient_Adam(vector<double> & x, vector<double> & grad){
		time += 1;
		for (int i=0; i<n_size; i++){
			grad_squared[i] = grad[i]*grad[i];
		}
		for (int i=0; i<n_size; i++){
			m[i] = beta_1*m[i]+(1-beta_1)*grad[i];
		}
		for (int i=0; i<n_size; i++){
			v[i] = beta_2*v[i]+(1-beta_2)*grad_squared[i];
		}
		for (int i=0; i<n_size; i++){
			m_hat[i] = m[i]/(1-pow(beta_1, time));
		}
		for (int i=0; i<n_size; i++){
			v_hat[i] = v[i]/(1-pow(beta_2, time));
		}
		for (int i=0; i<n_size; i++){
			x[i] = x[i]-alpha*(m_hat[i]/(std::sqrt(v_hat[i])+adam_eps));
		}
	}



	void step_gradient_Adamax(vector<double> & x, vector<double> & grad){
		time += 1;

		vector<double> pos_grad = grad;
		for(int i = 0; i < grad.size(); i++)
		{
			pos_grad[i] = abs(grad[i]);
		} 
		infinity_norm = *std::max_element(std::begin(pos_grad), std::end(pos_grad));

		for (int i=0; i<n_size; i++){
			m[i] = beta_1*m[i]+(1-beta_1)*grad[i];
		}
		for (int i=0; i<n_size; i++){
			v[i] = std::max(beta_2*v[i], infinity_norm);
		}
		for (int i=0; i<n_size; i++){
			x[i] = x[i]-(alpha/(1-pow(beta_1, time)))*(m[i]/v[i]);
		}
	}


	void step_gradient_RMS_prop(vector<double> & x, vector<double> & grad){
		for (int i=0; i<n_size; i++){
			rms_E_g2[i] = 0.9*rms_E_g2[i]+0.1*grad[i]*grad[i];
		}


		for (int i=0; i<n_size; i++){
			x[i] = x[i] - (0.001/(std::sqrt(rms_E_g2[i]+1.0e-9)))*grad[i];
		}

	}



};











