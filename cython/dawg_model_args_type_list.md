# List of dawg's model arguments, order of types matters when it gets # sent to CPP. Should match parameters in 'addModelArgument' in
# PyDawg.pyx
# The hard part is some of the strings are really lists,
# and some are just strings.
# Example? Tree is always a string,
# but 'substitution_params' can be a string or list
# This means we need to always do a find for a ','
# Even more complicated is the list could be of type double or # string (string of strings)

string name,
string inherits_from,
string substitution_model,
string substitution_params,
string substitution_freqs,
string substitution_rate_model,
string substitution_rate_params,

string indel_model_insertion,
string indel_params_insertion,
string indel_rate_insertion,
unsigned int indel_max_insertion,
string indel_model_deletion,
string indel_params_deletion,
string indel_rate_deletion,
unsigned int indel_max_deletion,

string tree,
string tree_model,
string tree_params,
double tree_scale,

unsigned int root_length,
string root_sequence,
string root_rates,
unsigned int root_code,
unsigned int root_segment,
bint root_gapoverlap,

bint output_rna,
bint output_lowercase,
bint output_keepempty,
bint output_markins
