//
// Created by tyler on 4/7/17.
//
#include "mesh/Mesh.hpp"
#include "utils/Range.hpp"
#include <fstream>
#include <string>
#include <sstream>

YAFEL_NAMESPACE_OPEN

static CellType gmsh_to_CellType(int);

static void string_split(const std::string &s, char delim,
                         std::vector<std::string> &elems)
{
    std::stringstream ss(s);
    std::string word;
    while (std::getline(ss, word, delim)) {
        elems.push_back(word);
    }

    return;
}


void Mesh::parse_gmsh(const std::string &fname)
{

    std::ifstream in(fname);

    if (!in.good()) {
        throw (std::runtime_error("Mesh::parse_gmsh: bad std::istream"));
    }

    typedef enum action
    {
        ACTION_UNSET,
        PARSE_NODES,
        PARSE_ELEMENTS
    } action_t;

    action_t currentAction = ACTION_UNSET;

    int offset{0};
    std::string line;
    std::vector<std::string> words;
    while (!in.eof()) {

        std::getline(in, line);
        if (in.eof()) {
            break;
        }

        words.clear();
        string_split(line, ' ', words);

        if (words[0]=="$Nodes") {
            currentAction = PARSE_NODES;
            std::getline(in, line);
            int nNodes = std::stoi(line);
            geometryNodes_.resize(nNodes);
            continue;
        } else if (words[0] == "$EndNodes") {
            currentAction = ACTION_UNSET;
            continue;
        } else if (words[0] == "$Elements") {
            currentAction = PARSE_ELEMENTS;
            std::getline(in, line);
            int nElems = std::stoi(line);
            cellNodes_.reserve(4 * nElems); //decent guess?
            cellOffsets_.resize(nElems + 1);
            cellTags_.resize(nElems);
            cellTypes_.resize(nElems);
            continue;
        } else if (words[0] == "$EndElements") {
            currentAction = ACTION_UNSET;
            continue;
        }

        int id = -1;
        switch (currentAction) {
            case ACTION_UNSET:
                break;

            case PARSE_NODES: {
                id = std::stoi(words[0]);
                coordinate<> node;
                for (auto i : IRange(0,3)) {
                    node(i) = std::stof(words[i + 1]);
                }
                geometryNodes_[id - 1] = node;
                break;
            }
            case PARSE_ELEMENTS: {
                id = std::stoi(words[0]);
                cellTypes_[id - 1] = gmsh_to_CellType(std::stoi(words[1]));
                int ntags = std::stoi(words[2]);
                int nodes_in_el = words.size() - ntags - 3;
                std::vector<size_type> tag(ntags, 0);
                for (int i = 3; i < 3 + ntags; ++i) {
                    tag[i - 3] = std::stoi(words[i]);
                }
                cellTags_[id - 1] = tag;
                for (size_type i = 3 + ntags; i < static_cast<int>(words.size()); ++i) {
                    cellNodes_.push_back(std::stoi(words[i]) - 1); //using 0-based node numbering from gmsh 1-based
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