/*  Dawg - DNA Assembly with Gaps - Simulating Sequence Evolution
    Copyright (c) 2018 Juan J. Garcia Mesa
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#pragma once
#ifndef DAWG_ERROR_H
#define DAWG_ERROR_H

#include <iostream>
#include <system_error>

// custom error code enum
enum class dawg_error {
    // i/o					//object needed to display error message
    open_input_file_fail = 10,			//const char*
    open_output_file_fail = 11,
    unknown_output_format = 12,			//std::string
    parse_input_file_fail = 13,			//const char*
    parsing_failed = 14,

    // indel model
    indel_maximum_size_out_range = 20,
    indel_model_no_type = 21,

    // rate model
    rate_no_model = 30,				//std::string&
    rate_gamma_empty_invariant = 31,
    rate_heterogeneous_creation_fail = 32,

    // subst model
    subst_model_gtr_negative_param = 40,	//int, It1
    subst_model_gtr_insufficient_params = 41,
    subst_model_aagtr_negative_param = 42,
    subst_model_aagtr_number_requirement = 43,
    subst_model_codgtr_negative_param = 44,	//int, It1
    subst_model_negative_frequency = 45,	//const char*, It1
    subst_model_no_name = 46,			//const char*
    subst_model_creation_fail = 47,

    // configuration
    bad_configuration = 50,
    configuration_section_fail = 51,		//first->name (?)

    // section
    section_not_found = 60,			//std::string, std::string
    section_already_specified = 61,		//std::string

    // miscellaneous
    invalid_genetic_code = 71,
    invalid_sequence_type = 72,
    tree_node_already_used = 73,		//it->label (?)
};

namespace std {
    template <>
	struct is_error_code_enum<dawg_error> : public true_type {};
}// end namespace std

namespace dawg {

    class dawg_error_category : public std::error_category {
    public:
	// return short descriptive name for the category
	virtual const char *name() const noexcept override final {return "Error custom category"; }
	// return enum description
	// TODO: add parameter(s) when needed for error message
	virtual std::string message(int ev) const override final {
	    switch (static_cast<dawg_error>(ev)) {
		case dawg_error::open_input_file_fail:
		    return "Unable to open input file __.";
		case dawg_error::open_output_file_fail:
		    return "Unable to open output file __.";
		case dawg_error::unknown_output_format:
		    return "Unknown output format __.";
		case dawg_error::parse_input_file_fail:
		    return "Unable to parse input __.";
		case dawg_error::parsing_failed:
		    return "Parsing failed.";
		case dawg_error::indel_maximum_size_out_range:
		    return "Maximum indel size is out of range.";
		case dawg_error::indel_model_no_type:
		    return "Invalid indel model; no model type specified.";
		case dawg_error::rate_no_model:
		    return "Invalid rate model; no model named __.";
		case dawg_error::rate_gamma_empty_invariant:
		    return "Invalid rate model; gamma-invariant requires at least 1 parameter.";
		case dawg_error::rate_heterogeneous_creation_fail:
		    return "Heterogeneous rate model could not be created.";
		case dawg_error::subst_model_gtr_negative_param:
		    return "Invalid subst model; gtr parameter __ is not >= 0.";
		case dawg_error::subst_model_gtr_insufficient_params:
		    return "Invalid subst model; gtr requires six parameters.";
		case dawg_error::subst_model_aagtr_negative_param:
		    return "Invalid subst model; aagtr parameter __ is not >= 0.";
		case dawg_error::subst_model_aagtr_number_requirement:
		    return "Invalid subst model; aagtr requires 190 parameters.";
		case dawg_error::subst_model_codgtr_negative_param:
		    return "Invalid subst model; codgtr parameter __ is not >= 0.";
		case dawg_error::subst_model_negative_frequency:
		    return "Invalid subst model; __ frequency __ is not >= 0.";
		case dawg_error::subst_model_no_name:
		    return "Invalid subst model; no model named __.";
		case dawg_error::subst_model_creation_fail:
		    return "Substitution model could not be created.";
		case dawg_error::bad_configuration:
		    return "Bad configuration.";
		case dawg_error::configuration_section_fail:
		    return "Configuration section __ failed to process.";
		case dawg_error::section_not_found:
		    return "Section __ not found (inherited by __.)";
		case dawg_error::section_already_specified:
		    return "Section __ specified more than once.";
		case dawg_error::invalid_genetic_code:
		    return "Invalid genetic code.";
		case dawg_error::invalid_sequence_type:
		    return "Invalid sequence type.";
		case dawg_error::tree_node_already_used:
		    return "Invalid tree; node label __ used more than once by Tree.Tree.";
		default:
		    return "Unrecognized error.";
	    }
	}
    }; // end of dawg_error_category class
} // end of namespace dawg

    #define THIS_MODULE_API_DECL extern inline
    THIS_MODULE_API_DECL const dawg::dawg_error_category &dawg_error_category() {
	static dawg::dawg_error_category c;
	return c;
    }

    inline std::error_code make_error_code(dawg_error e) {
	return {static_cast<int>(e), dawg_error_category()};
    }

#endif //DAWG_ERROR_H

