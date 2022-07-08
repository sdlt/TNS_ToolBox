void init_prediction(char pkfile[]);
void free_prediction(void);
int get_prediction(double in_s_min, double in_s_max, int in_n_rsd, double *out, 
		   double in_f, double in_b1, double in_b2, double in_bg, double in_bt, double in_sigv, double in_alpha_par, double in_alpha_per);
int get_prediction_LL(double in_s_min, double in_s_max, int in_n_rsd, double *out, 
		      double in_f, double in_b1, double in_b2, double in_sigv, double in_alpha_par, double in_alpha_per);
