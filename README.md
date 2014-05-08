nonlocal_qif
============

Neural field models are commonly used to investigate spatio-temporal dynamics in the cortex. The derivation of such models from microscopic networks of spiking neurons is not clear and typically involves a number of approximations which may not be fully correct. Here we present an exact derivation of a neural field model from a network of quadratic integrate-and-fire neurons, the canonical form of a Type 1 neuron. We analize the model by means of a linear stability analysis and  demonstrate the existence of a number pattern forming instabilities.


About Structures:
¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯

	Data structure contains two structures, one related to he theoretical eqs, 
	t_FR, and another related to the spìking neuron simulation, T_qif.

	Once a "data" structure is initialized with 
	     t_data *d;
	you can assign the correspondent structure to its elements, e.g.:
	     T_FR *FR;
	     FR = (T_FR*) calloc(N,sizeof(T_FR));
	     
	     /* Now we assign this structure to the respective element in "data" */
	     d->FR = &FR;

	If we want to access to this information we must use the next sintax:
	     (*(d->FR) + 1)->r = 0.3;
	or
	     (*(d->FR))[1].r = 0.3;
	which will assign the value of 0.3 to the element r of the structure FR,
	using the path to its memory from the "data" structure.