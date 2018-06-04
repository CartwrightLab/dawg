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

#include <string>
#include <iostream>
#include <system_error>

// custom error code enum
enum class dawg_error {
    // i/o
    open_input_file_fail = 10,
    open_output_file_fail = 11,
    unknown_output_format = 12,
    parse_input_file_fail = 13,
    parsing_failed = 14,

    // indel model
    indel_model_no_type = 20,

    // configuration
    bad_configuration = 30,
    configuration_section_fail = 31,

    // section
    section_not_found = 40,
    section_already_specified = 41,
    section_sequence_type_invalid = 42,

    // miscellaneous
    invalid_genetic_code = 51,
    invalid_sequence_type = 52,
    creation_fail = 53,
    invalid_value = 54,
    model_no_name = 55,
    segment_extra_root = 56,
    no_sequence = 57,
    invalid_tree = 58,
    unexpected = 59,
    param_missing = 60,
    invalid_user_sequence = 61,
};

namespace std {
    template <>
	struct is_error_code_enum<dawg_error> : public true_type {};
}// end namespace std

namespace dawg {

    // declare dawg_error_t struct
    struct dawg_error_t : std::exception{
	std::error_code ecode;
	std::string information;

	dawg_error_t(dawg_error e) {
	    ecode = e;
	    information = "";
	}
	dawg_error_t(dawg_error e, std::string str) {
	    ecode = e;
	    information = str;
	}
	std::string getInformation() {
	    return information;
	}

    };

    class dawg_error_category : public std::error_category {
    public:

	// return short descriptive name for the category
	virtual const char *name() const noexcept override final {return "Dawg error category"; }

	// return enum description
	virtual std::string message(int ev) const override final {
	    switch (static_cast<dawg_error>(ev)) {
		case dawg_error::open_input_file_fail:	// i/o
		    return "Unable to open input file ";
		case dawg_error::open_output_file_fail:
		    return "Unable to open output file ";
		case dawg_error::unknown_output_format:
		    return "Unknown output format ";
		case dawg_error::parse_input_file_fail:
		    return "Unable to parse input ";
		case dawg_error::parsing_failed:
		    return "Parsing failed.";
		case dawg_error::indel_model_no_type:
		    return "Invalid indel model; no model type specified.";
		case dawg_error::param_missing:
		    return "Param is missing; ";
		case dawg_error::bad_configuration: // configuration
		    return "Bad configuration.";
		case dawg_error::configuration_section_fail:
		    return "Configuration section failed to process ";
		case dawg_error::section_not_found: // section
		    return "Section not found ";
		case dawg_error::section_already_specified:
		    return "Section specified more than once ";
		case dawg_error::section_sequence_type_invalid:
		    return "The sequence type or format options of a section is\
			different that its segment.";
		case dawg_error::invalid_genetic_code:	// miscellaneous
		    return "Invalid genetic code.";
		case dawg_error::invalid_sequence_type:
		    return "Invalid sequence type.";
		case dawg_error::creation_fail:
		    return "Could not be created: ";
		case dawg_error::invalid_value:
		    return "Invalid parameter value; ";
		case dawg_error::model_no_name:
		    return "No model named ";
		case dawg_error::segment_extra_root:
		    return "Extra root; ";
		case dawg_error::no_sequence:
		    return "No sequences to align.";
		case dawg_error::invalid_tree:
		    return "Invalid tree; ";
		case dawg_error::unexpected:
		    return "Unexpected error.";
		case dawg_error::invalid_user_sequence:
		    return "Invalid user sequence; ";
		default:
		    return "Unrecognized error.";
	    }
	}
    }; // end of dawg_error_category class
} // end of namespace dawg

    inline const dawg::dawg_error_category &dawg_error_category() {
	static dawg::dawg_error_category c;
	return c;
    }

    //inline std::error_code make_error_code(dawg_error e) {
    inline std::error_code make_error_code(dawg_error e) {
	return {static_cast<int>(e), dawg_error_category()};
    }

#endif //DAWG_ERROR_H

