/*
 * Netgen.hpp
 *
 *  Created on: Apr 16, 2012
 *      Author: kremer
 */

#ifndef NETGEN_HPP_
#define NETGEN_HPP_

#include <boost/spirit/include/qi.hpp>

#include <boost/bind.hpp>

#include "../MeshGenerator.hpp"

#include <iostream>

namespace qi = boost::spirit::qi;
namespace spirit = boost::spirit;

void print() {
    static int c = 0;
    std::cerr << "Hey there!" << c++ << std::endl;
}

template <typename Iterator>
class netgen_grammar : public qi::grammar<Iterator/*, qi::space_type*/> {
public:
    netgen_grammar(MeshGenerator& _generator) :
        netgen_grammar::base_type(content),
        generator_(_generator) {

        content = node_section_header >> *node >>
                  element_section_header >> *element >>
                  face_section_header >> *face;

        node_section_header = qi::int_ /* Number of vertices */ >> spirit::eol;

        node = *space >> qi::double_[boost::bind(&MeshGenerator::add_vertex_component, &generator_, ::_1)] >>
               *space >> qi::double_[boost::bind(&MeshGenerator::add_vertex_component, &generator_, ::_1)] >>
               *space >> qi::double_[boost::bind(&MeshGenerator::add_vertex_component, &generator_, ::_1)] >> spirit::eol;

        element_section_header = qi::int_[boost::bind(&MeshGenerator::set_num_cells, &generator_, ::_1)] /* Number of tetrahedra */ >>
                spirit::eol;;

        element = *space >> qi::int_ >>
                  *space >> qi::int_[boost::bind(&MeshGenerator::add_cell_vertex, &generator_, ::_1)] >>
                  *space >> qi::int_[boost::bind(&MeshGenerator::add_cell_vertex, &generator_, ::_1)] >>
                  *space >> qi::int_[boost::bind(&MeshGenerator::add_cell_vertex, &generator_, ::_1)] >>
                  *space >> qi::int_[boost::bind(&MeshGenerator::add_cell_vertex, &generator_, ::_1)] >>
                  spirit::eol;;

        face_section_header = qi::int_ /* Number of faces */ >> spirit::eol;

        face = *space >> qi::int_ >>
               *space >> qi::int_ >>
               *space >> qi::int_ >>
               *space >> qi::int_ >>
               spirit::eol;

        space = spirit::ascii::space - spirit::eol;
    }

private:

    qi::rule<Iterator/*, qi::space_type*/> node_section_header, element_section_header, face_section_header,
                                       node, element, face, space, content;

    MeshGenerator& generator_;
};

#endif /* NETGEN_HPP_ */
