#ifndef flopoco_random_utils_anneal_asa_hpp
#define flopoco_random_utils_anneal_asa_hpp
	
#define LONG_INT long int
#define ALLOC_INT int
	
  typedef struct {
    LONG_INT Limit_Acceptances;
    LONG_INT Limit_Generated;
    int Limit_Invalid_Generated_States;
    double Accepted_To_Generated_Ratio;

    double Cost_Precision;
    int Maximum_Cost_Repeat;
    int Number_Cost_Samples;
    double Temperature_Ratio_Scale;
    double Cost_Parameter_Scale_Ratio;
    double Temperature_Anneal_Scale;
/*#if USER_INITIAL_COST_TEMP
    double *User_Cost_Temperature;
#endif*/

    int Include_Integer_Parameters;
    int User_Initial_Parameters;
    ALLOC_INT Sequential_Parameters;
    double Initial_Parameter_Temperature;
/*#if RATIO_TEMPERATURE_SCALES
    double *User_Temperature_Ratio;
#endif
#if USER_INITIAL_PARAMETERS_TEMPS
    double *User_Parameter_Temperature;
#endif*/

    int Acceptance_Frequency_Modulus;
    int Generated_Frequency_Modulus;
    int Reanneal_Cost;
    int Reanneal_Parameters;

    double Delta_X;
/*#if DELTA_PARAMETERS
    double *User_Delta_Parameter;
#endif*/
    int User_Tangents;
    int Curvature_0;

/*#if QUENCH_PARAMETERS
    double *User_Quench_Param_Scale;
#endif
#if QUENCH_COST
    double *User_Quench_Cost_Scale;
#endif*/

    LONG_INT N_Accepted;
    LONG_INT N_Generated;
    int Locate_Cost;
    int Immediate_Exit;

    double *Best_Cost;
    double *Best_Parameters;
    double *Last_Cost;
    double *Last_Parameters;

/*
#if OPTIONAL_DATA_DBL
    ALLOC_INT Asa_Data_Dim_Dbl;
    double *Asa_Data_Dbl;
#endif
#if OPTIONAL_DATA_INT
    ALLOC_INT Asa_Data_Dim_Int;
    LONG_INT *Asa_Data_Int;
#endif*/
//#if OPTIONAL_DATA_PTR
    size_t Asa_Data_Dim_Ptr;
    void *Asa_Data_Ptr;
//#endif
/*
#if USER_ASA_OUT
    char *Asa_Out_File;
#endif
*/
    /* Keep OPTIONS_TMP in parameter lists in asa_usr.[ch] as they are
     * needed if using recursively, e.g., with SELF_OPTIMIZE=TRUE.
     * Make (USER_DEFINE *) casts explicit within functions. */
	 /*
#if USER_COST_SCHEDULE
#if HAVE_ANSI
    double (*Cost_Schedule) (double current_cost_temperature,
                             const void *OPTIONS_TMP);
#else                           // HAVE_ANSI 
    double (*Cost_Schedule) ();
#endif                          // HAVE_ANSI
#endif
#if USER_ACCEPT_ASYMP_EXP
    double Asymp_Exp_Param;
#endif
#if USER_ACCEPTANCE_TEST
#if HAVE_ANSI
    void (*Acceptance_Test) (double cost,
                             double *parameter_minimum,
                             double *parameter_maximum,
                             ALLOC_INT * number_parameters,
                             const void *OPTIONS_TMP);
#else                           // HAVE_ANSI 
    void (*Acceptance_Test) ();
#endif                          // HAVE_ANSI
    int User_Acceptance_Flag;
    int Cost_Acceptance_Flag;
    double Cost_Temp_Curr;
    double Cost_Temp_Init;
    double Cost_Temp_Scale;
    double Prob_Bias;
    LONG_INT *Random_Seed;
#endif
#if USER_GENERATING_FUNCTION
#if HAVE_ANSI
    double (*Generating_Distrib) (LONG_INT * seed,
                                  ALLOC_INT * parameter_dimension,
                                  ALLOC_INT index_v,
                                  double temperature_v,
                                  double init_param_temp_v,
                                  double temp_scale_params_v,
                                  double parameter_v,
                                  double parameter_range_v,
                                  double *last_saved_parameter,
                                  const void *OPTIONS_TMP);
#else                           // HAVE_ANSI
    double (*Generating_Distrib) ();
#endif                          // HAVE_ANSI
#endif
#if USER_REANNEAL_COST
#if HAVE_ANSI
    int (*Reanneal_Cost_Function) (double *cost_best,
                                   double *cost_last,
                                   double *initial_cost_temperature,
                                   double *current_cost_temperature,
                                   const void *OPTIONS_TMP);
#else                           // HAVE_ANSI
    int (*Reanneal_Cost_Function) ();
#endif                          // HAVE_ANSI 
#endif
#if USER_REANNEAL_PARAMETERS
#if HAVE_ANSI
    double (*Reanneal_Params_Function) (double current_temp,
                                        double tangent,
                                        double max_tangent,
                                        const void *OPTIONS_TMP);
#else                           // HAVE_ANSI
    double (*Reanneal_Params_Function) ();
#endif                          // HAVE_ANSI
#endif
#if ASA_SAMPLE
    double Bias_Acceptance;
    double *Bias_Generated;
    double Average_Weights;
    double Limit_Weights;
#endif
#if ASA_QUEUE
    ALLOC_INT Queue_Size;
    double *Queue_Resolution;
#endif
#if ASA_RESOLUTION
    double *Coarse_Resolution;
#endif
#if FITLOC
    int Fit_Local;
    int Iter_Max;
    double Penalty;
#endif
#if MULTI_MIN
    int Multi_Number;
    double *Multi_Cost;
    double **Multi_Params;
    double *Multi_Grid;
    int Multi_Specify;
#endif
#if ASA_PARALLEL
    int Gener_Mov_Avr;
    LONG_INT Gener_Block;
    LONG_INT Gener_Block_Max;
#endif
#if ASA_SAVE
    ALLOC_INT Random_Array_Dim;
    double *Random_Array;
#endif
*/
    int Asa_Recursive_Level;
/*
#if ASA_FUZZY
    int NoOfSamples;
    double ThresholdDeviation;
    double Threshold1;
    double Performance_Target;
    double Factor_a;
#endif
*/
  } asa_USER_DEFINES;	
	
#undef LONG_INT
#undef ALLOC_INT
  
typedef double (*asa_user_cost_function_t)(
  double *, double *, double *, double *, double *, int *, int *, int *, int *, asa_USER_DEFINES *)
);
typedef double (*user_random_generator_t)(long int *);
  
double asa (
	asa_user_cost_function_t cost_function,
	asa_user_random_generator_t user_random_generator,
	long int * rand_seed,
	 double *parameter_initial_final, double *parameter_minimum,
	 double *parameter_maximum, double *tangents, double *curvature,
	 int * number_parameters, int *parameter_type,
	 int *valid_state_generated_flag, int *exit_status,
	 asa_USER_DEFINES * OPTIONS
);
  
namespace flopoco
{
namespace random
{
namespace detail
{

/*
	struct TF{
		unsigned Arity() const;
		// Return (is_integer,should_anneal) for the parameter
		std::pair<bool,bool> ParameterProperties(unsigned param) const;
		// Return lower and upper bounds for the parameter
		std::pair<double,double> ParameterBounds(unsigned param) const;
		double operator()(const std::vector<double> &v) const;
	};
*/
template<class TF>
class AnnealASAImpl
{
	TF m_f;
	unsigned m_arity;
	
	static cost_function (
		double *x,
		   double *parameter_lower_bound,
		   double *parameter_upper_bound,
		   double *cost_tangents,
		   double *cost_curvature,
		    int* parameter_dimension,
		   int *parameter_int_real,
		   int *cost_flag, int *exit_code,
			asa_USER_DEFINES * USER_OPTIONS
	){
		AnnealASAImpl *me=(AnnealASAImpl *)((asa_USER_DEFINES->Asa_Data_Ptr);
		std::vector<double> point(x, x+m_arity);
		*cost_flag = true;
		return me->m_f(point);
	}
	
	double CostFunction(const std::vector<double> &x) const
	{ m_f(x); }
	
	std::vector<double> Execute(
		const std::vector<double> &init,
		const std::vector<std::pair<double,double> > &bounds
	){
		m_arity=m_f.Arity();
		std::vector<double> lower_bounds(m_arity), upper_bounds(m_arity);
		std::vector<int> parameter_props(m_arity);
		
		for(unsigned i=0;i<m_arity;i++){
			lower_bounds[i]=m_f.ParameterBounds(i).first;
			upper_bounds[i]=m_f.ParameterBounds(i).second;
			parameter_props[i]=m_f.ParameterProperties(i).first ? 1 : -1;
			parameter_props[i] *= m_f.ParameterProperties(i).second ? 1 : 2;
		}
		
		asa_USER_DEFINES defs;
		
		defs.Limit_Acceptances = 1000;
		defs.Limit_Generated = 99999;
		defs.Limit_Invalid_Generated_States = 1000;
		/* defs.Accepted_To_Generated_Ratio = 1.0E-6; */
		defs.Accepted_To_Generated_Ratio = 1.0E-4;

		defs.Cost_Precision = 1.0E-18;
		defs.Maximum_Cost_Repeat = 5;
		defs.Number_Cost_Samples = 5;
		defs.Temperature_Ratio_Scale = 1.0E-5;
		defs.Cost_Parameter_Scale_Ratio = 1.0;
		defs.Temperature_Anneal_Scale = 100.0;

		defs.Include_Integer_Parameters = FALSE;
		defs.User_Initial_Parameters = FALSE;
		defs.Sequential_Parameters = -1;
		defs.Initial_Parameter_Temperature = 1.0;

		defs.Acceptance_Frequency_Modulus = 100;
		defs.Generated_Frequency_Modulus = 10000;
		defs.Reanneal_Cost = 1;
		defs.Reanneal_Parameters = TRUE;

		defs.Delta_X = 0.001;
		defs.User_Tangents = FALSE;
		defs.Curvature_0 = FALSE;

		int exit_state, valid_state_generated_flag;
		asa (
			asa_user_cost_function_t cost_function,
			asa_user_random_generator_t user_random_generator,
			long int * rand_seed,
			 double *parameter_initial_final, double *parameter_minimum,
			 double *parameter_maximum,
			 double *tangents,
			 double *curvature,
			 &m_arity,
			 &parameter_props[0],
			 &valid_state_generated_flag,
			 &exit_status,
			 &defs /*asa_USER_DEFINES * OPTIONS*/
		);
	}
};	
	
};
}; // random
}; // flopoco

#endif
