//
// Created by tyler on 4/7/17.
//
#include "mesh/Mesh.hpp"
#include "utils/Range.hpp"
#include <fstream>
#include <string>
//#include <sstream>
#include <string_view>
//#include <charconv> //for std::from_chars, which is not yet in gcc as of v7.1

YAFEL_NAMESPACE_OPEN

static CellType gmsh_to_CellType(int);

/*
static void string_split(const std::string &s, char delim,
                         std::vector<std::string> &elems)
{
    std::stringstream ss(s);
    std::string word;
    while (std::getline(ss, word, delim)) {
        elems.push_back(word);
    }

}*/

namespace {

void string_split(std::string_view s, char delim,
                  std::vector<std::string_view> &elems) {

    auto pos = s.find_first_of(delim);
    std::string::size_type start = 0;
    while(pos != std::string::npos) {
        if(pos != start) {
            elems.push_back(s.substr(start, pos-start));
        }
        start = pos+1;

        pos = s.find_first_of(delim, start);
    }
    if(start < s.length()) {
        elems.push_back(s.substr(start, s.length()-start));
    }

}

auto svtoi(std::string_view s) {
    int val{0};
    //This is the right way to do it in C++17, but gcc hasn't implemented it
    //std::from_chars(s.data(), std::next(s.data(), s.size()), val);

    //Super-fragile integer parsing, cannot handle badly-formatted input
    //Evaluate using Horner's method for polynomial evaluation.
    auto pos = s.data();
    if(!pos) {
        return val;
    }
    while(pos < s.end()) {
        val = 10*val + (*pos - '0');
        ++pos;
    }
    return val;
}

auto svtod(std::string_view s) {
    double val{0};

    //This is the right way to do it in C++17, but gcc hasn't implemented it
    //std::from_chars(s.data(), std::next(s.data(), s.size()), val);

    //Shitty fallback. Praying to the gods of small-string optimization here...
    std::string S(s);
    val = std::stod(S);

    return val;
}
}


void Mesh::parse_gmsh(const std::string &fname)
{

    std::ifstream in(fname);

    if (!in.good()) {
        throw (std::runtime_error("Mesh::parse_gmsh: bad std::istream"));
    }

    enum class ParsingAction
    {
        ACTION_UNSET,
        PARSE_NODES,
        PARSE_ELEMENTS
    };

    ParsingAction currentAction = ParsingAction::ACTION_UNSET;

    int offset{0};
    std::string line;
    std::vector<std::string_view> words;
    while (!in.eof()) {

        std::getline(in, line);
        if (in.eof()) {
            break;
        }

        words.clear();
        string_split(line, ' ', words);

        if (words[0]=="$Nodes") {
            currentAction = ParsingAction::PARSE_NODES;
            std::getline(in, line);
            int nNodes = std::stoi(line);
            geometryNodes_.resize(nNodes);
            continue;
        } else if (words[0] == "$EndNodes") {
            currentAction = ParsingAction::ACTION_UNSET;
            continue;
        } else if (words[0] == "$Elements") {
            currentAction = ParsingAction::PARSE_ELEMENTS;
            std::getline(in, line);
            int nElems = std::stoi(line);
            cellNodes_.reserve(4 * nElems); //decent guess?
            cellOffsets_.resize(nElems + 1);
            cellTags_.resize(nElems);
            cellTypes_.resize(nElems);
            continue;
        } else if (words[0] == "$EndElements") {
            currentAction = ParsingAction::ACTION_UNSET;
            continue;
        }

        int id = -1;
        switch (currentAction) {
            case ParsingAction::ACTION_UNSET:
                break;

            case ParsingAction::PARSE_NODES: {
                id = svtoi(words[0]);
                coordinate<> node;
                for (auto i : IRange(0,3)) {
                    node(i) = svtod(words[i + 1]);
                }
                geometryNodes_[id - 1] = node;
                break;
            }
            case ParsingAction::PARSE_ELEMENTS: {
                id = svtoi(words[0]);
                cellTypes_[id - 1] = gmsh_to_CellType(svtoi(words[1]));
                int ntags = svtoi(words[2]);
                int nodes_in_el = words.size() - ntags - 3;
                std::vector<size_type> tag(ntags, 0);
                for (int i = 3; i < 3 + ntags; ++i) {
                    tag[i - 3] = svtoi(words[i]);
                }
                cellTags_[id - 1] = tag;
                for (size_type i = 3 + ntags; i < static_cast<int>(words.size()); ++i) {
                    cellNodes_.push_back(svtoi(words[i]) - 1); //using 0-based node numbering from gmsh 1-based
                }
                cellOffsets_[id - 1] = offset;
                offset += nodes_in_el;
                break;
            }
        }

    } //end while(!in.eof())
    cellOffsets_.back() = offset;


}


static CellType gmsh_to_CellType(int id)
{
    switch (id) {
        case 1:
            return CellType::Line2;
        case 2:
            return CellType::Tri3;
        case 3:
            return CellType::Quad4;
        case 4:
            return CellType::Tet4;
        case 5:
            return CellType::Hex8;
        case 15:
            return CellType::Point1;
        default:
            std::cerr << "Unknown gmsh type. (" << id << ")" << std::endl;
            return CellType::None;
    }
}


YAFEL_NAMESPACE_CLOSE