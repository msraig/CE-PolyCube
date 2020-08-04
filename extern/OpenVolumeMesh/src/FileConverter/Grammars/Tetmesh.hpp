/*
 * Tetmesh.hpp
 *
 *  Created on: Mar 20, 2012
 *      Author: kremer
 */

#include <boost/spirit/include/qi.hpp>

#include <boost/bind.hpp>

#include "../MeshGenerator.hpp"

#include <iostream>

#ifndef TETMESH_HPP_
#define TETMESH_HPP_

namespace qi = boost::spirit::qi;
namespace spirit = boost::spirit;

template <typename Iterator>
class tetmesh_grammar : public qi::grammar<Iterator, qi::space_type> {
public:
    tetmesh_grammar(MeshGenerator& _generator) :
        tetmesh_grammar::base_type(content),
        generator_(_generator) {

        content = node_section_header >> *node >> element_section_header >> *element;

        node_section_header = spirit::lit("Vertices") >> qi::int_;

        node = qi::double_[boost::bind(&MeshGenerator::add_vertex_component, &generator_, ::_1)] >>
               qi::double_[boost::bind(&MeshGenerator::add_vertex_component, &generator_, ::_1)] >>
               qi::double_[boost::bind(&MeshGenerator::add_vertex_component, &generator_, ::_1)] >>
               qi::double_;

        element_section_header = spirit::lit("Tetrahedra") >>
                qi::int_[boost::bind(&MeshGenerator::set_num_cells, &generator_, ::_1)];

        element = qi::int_[boost::bind(&MeshGenerator::add_cell_vertex, &generator_, ::_1)] >>
                  qi::int_[boost::bind(&MeshGenerator::add_cell_vertex, &generator_, ::_1)] >>
                  qi::int_[boost::bind(&MeshGenerator::add_cell_vertex, &generator_, ::_1)] >>
                  qi::int_[boost::bind(&MeshGenerator::add_cell_vertex, &generator_, ::_1)];
    }

private:

    qi::rule<Iterator, qi::space_type> node_section_header, element_section_header, node, element, content;

    MeshGenerator& generator_;
};

#endif /* TETMESH_HPP_ */
