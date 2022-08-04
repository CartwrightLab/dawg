#pragma once
#ifndef DAWG_TRICK_PARSE_H
#define DAWG_TRICK_PARSE_H
/****************************************************************************
 *  Copyright (C) 2009-2010 Reed A. Cartwright, PhD <reed@scit.us>          *
 *  Copyright (C) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>   *
 ****************************************************************************/
#define RYML_SINGLE_HDR_DEFINE_NOW
#include <iterator>
#include <rapidyaml.hpp>

#include "dawg/trick.h"

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4127)
#endif

#ifdef _MSC_VER
#pragma warning(pop)
#endif
namespace dawg {

template <typename Char, typename Traits>
inline std::string slurp(std::basic_istream<Char, Traits> &is) {
    if(!is) {
        return {};
    }
    std::string content;
    is.seekg(0, std::ios::end);
    content.resize(is.tellg());
    is.seekg(0, std::ios::beg);
    is.read(&content[0], content.size());
    return content;
}

inline std::string to_string(ryml::csubstr substr) {
    return {substr.str, substr.len};
}

inline std::pair<std::string, trick::section::value_type> get_values(
    const ryml::NodeRef &node, std::string &label) {
    trick::section::value_type values;
    // create key
    std::string key = label + "." + dawg::to_string(node.key());
    if(node.num_children() == 0) {  // single value
        values = {dawg::to_string(node.val())};
    } else {  // array of values
        for(auto child : node.children()) {
            values.emplace_back(dawg::to_string(child.val()));
        }
    }

    return std::make_pair(key, values);
}

inline bool add_param(trick::section *sec, const ryml::NodeRef &node,
                      ryml::csubstr name = "") {
    ryml::NodeRef label;
    if(name == "") {
        label = node;
    } else {
        if(!node.has_child(name)) {
            return DAWG_ERROR(dawg::to_string(node.key())
                              << " '" << dawg::to_string(name)
                              << "' not found.");
        }
        label = node.find_child(name);
    }

    for(ryml::NodeRef const &subnode : label.children()) {
        if(subnode.is_map()) {
            for(ryml::NodeRef const &subsubnode : subnode.children()) {
                if(subsubnode.is_map()) {
                    for(ryml::NodeRef const &item : subsubnode.children()) {
                        // create key
                        std::string key{dawg::to_string(subnode.key()) + "." +
                                        dawg::to_string(subsubnode.key())};
                        sec->db.insert(get_values(item, key));
                    }
                } else {
                    // create key
                    std::string key;
                    if(dawg::to_string(node.key()) == "parts") {
                        key = "root";
                    } else {
                        key = dawg::to_string(subnode.key());
                    }
                    sec->db.insert(get_values(subsubnode, key));
                }
            }
        } else {
            // create key
            std::string key;
            if(dawg::to_string(node.key()) == "parts") {
                key = "root";
            } else {
                key = dawg::to_string(node.key());
            }
            sec->db.insert(get_values(subnode, key));
        }
    }
    return true;
}

inline void add_region(trick::section *sec, ryml::Tree &tree,
                       const ryml::NodeRef &region) {
    sec->name.assign(dawg::to_string(region.key()));

    // segment position (sequential)
    if(region.has_child("segment")) {
        trick::section::value_type seg{
            dawg::to_string(region["segment"].val())};
        sec->db.insert_or_assign("root.segment", seg);
    }

    // inherits
    if(region.has_child("inherits")) {
        sec->inherits = dawg::to_string(region["inherits"].val());
    }

    // trees
    if(region.has_child("tree")) {
        add_param(sec, tree["tree"], region["tree"].val());
    }

    // rules
    if(region.has_child("rule")) {
        add_param(sec, tree["rules"], region["rule"].val());
    }

    // parts
    if(region.has_child("part")) {
        add_param(sec, tree["parts"], region["part"].val());
    }
}

template <typename Char, typename Traits>
bool trick::parse_yaml(std::basic_istream<Char, Traits> &is) {
    // read and parse yaml config
    std::string content{slurp(is)};
    ryml::Tree tree = ryml::parse_in_place(ryml::to_substr(content));
    ryml::NodeRef root = tree.rootref();
    section *psec = &data.front();

    // first add simulation and output parameters
    if(root.has_child("sim")) {
        // add_global_param(psec, root["sim"]);
        add_param(psec, root["sim"]);
    }

    if(root.has_child("output")) {
        // add_global_param(psec, root["output"]);
        add_param(psec, root["output"]);
    }

    // iterate through all regions
    ryml::NodeRef regions = tree["regions"];
    for(ryml::NodeRef const &region : regions.children()) {
        if(region != regions.first_child()) {
            // add new empty section
            std::string previous{data.back().name};
            data.push_back(section());
            data.back().inherits = previous;
            psec = &data.back();
        }

        add_region(psec, tree, region);
    }
    return true;
}

}  // namespace dawg

#endif  // DAWG_TRICK_PARSE_H
