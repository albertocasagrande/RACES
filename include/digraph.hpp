/**
 * @file digraph.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Represent directed graphs
 * @version 0.1
 * @date 2023-06-03
 * 
 * @copyright Copyright (c) 2023
 * 
 * MIT License
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef __RACES_DIGRAPH__
#define __RACES_DIGRAPH__

#include <map>

namespace Races {

typedef unsigned int VertexId;

/**
 * @brief Labelled directed graphs
 * 
 * This class represents labelled directed graphs
 */
template<typename LABEL>
struct DiGraph
{
    typedef std::map<VertexId, LABEL> DstEdgeMap;

protected:
    std::map<VertexId, DstEdgeMap> edge_map; //!< Adjacency map

public:
    /**
     * @brief The empty directed graph constructor
     */
    DiGraph();

    /**
     * @brief The copy constructor
     * 
     * @param orig is the template object
     */
    DiGraph(const DiGraph<LABEL>& orig);

    /**
     * @brief Add a vertex to the directed graph
     * 
     * No changes occur when the edge is not part of the graph.
     * 
     * @param vertex is a vertex
     * @return the updated directed graph
     */
    DiGraph<LABEL>& add_vertex(const VertexId& vertex);

    /**
     * @brief Add an edge to the directed graph
     * 
     * This method adds an edge to the directed graph.
     * It raises an exception whenever the edge is already present.
     * 
     * @param src is the edge source
     * @param dst is the edge destination
     * @param label is the edge label
     * @return the updated directed graph
     */
    DiGraph<LABEL>& add_edge(const VertexId& src, const VertexId& dst, const LABEL& label);

    /**
     * @brief Remove an edge from the directed graph
     * 
     * No changes occur when the edge is not part of the graph.
     * 
     * @param src is the edge source
     * @param dst is the edge destination
     * @return the updated directed graph
     */
    DiGraph<LABEL>& remove_edge(const VertexId& src, const VertexId& dst);

    /**
     * @brief Test whether a vertex is in the graph
     * 
     * @param vertex is the vertex to be checked
     * @return true if and only if the vertex is in the graph
     */
    bool has_vertex(const VertexId& vertex) const;

    /**
     * @brief Test whether an edge is part of the graph 
     * 
     * @param src is the source of the edge to be checked
     * @param dst is the destination of the edge to be checked
     * @return true if and only if the edge is in the graph
     */
    bool has_edge(const VertexId& src, const VertexId& dst) const;

    /**
     * @brief Get the label of the edge
     * 
     * @param src is the edge source
     * @param dst is the edge destination
     * @return a constant reference to the edge label
     */
    const LABEL& get_edge_label(const VertexId& src, const VertexId& dst) const;

    /**
     * @brief Get the destination-label map of the edges leaving a vertex
     * 
     * @param src is the source of the edges 
     * @return the destination-label map of the edges leaving `src`
     */
    const DstEdgeMap& get_outgoing_edge_map(const VertexId& src) const;
};


/* Implementation */

template<typename LABEL>
DiGraph<LABEL>::DiGraph()
{}

template<typename LABEL>
DiGraph<LABEL>::DiGraph(const DiGraph<LABEL>& orig):
    edge_map(orig.edge_map)
{}

template<typename LABEL>
DiGraph<LABEL>& DiGraph<LABEL>::add_vertex(const VertexId& vertex)
{
    if (!has_vertex(vertex)) {

        // if the vertex is not in the graph, add it to the 
        // vertex map
        edge_map[vertex] = DiGraph<LABEL>::DstEdgeMap();
    }

    return *this;
}

template<typename LABEL>
DiGraph<LABEL>& DiGraph<LABEL>::add_edge(const VertexId& src, const VertexId& dst, const LABEL& label)
{
    auto& outgoing_edges = edge_map[src];

    if (outgoing_edges.count(dst)>0) {
        throw std::domain_error("The edge is already in the graph");
    }

    // add the labelled edge
    outgoing_edges[dst] = label;

    return *this;
}

template<typename LABEL>
DiGraph<LABEL>& DiGraph<LABEL>::remove_edge(const VertexId& src, const VertexId& dst)
{
    auto vertex_it = edge_map.find(src);

    // if the vertex is not in the graph
    if (vertex_it == std::end(edge_map)) {
        return *this;
    }

    auto edge_it = vertex_it.second.find(dst);

    // if the vertex does not contain the edge
    if (edge_it == std::end(vertex_it.second)) {
        return *this;
    }

    // remove the edge
    vertex_it.second.erase(edge_it);

    return *this;
}

template<typename LABEL>
inline bool DiGraph<LABEL>::has_vertex(const VertexId& vertex) const
{
    return edge_map.count(vertex)>0;
}

template<typename LABEL>
bool DiGraph<LABEL>::has_edge(const VertexId& src, const VertexId& dst) const
{
    const auto vertex_it = edge_map.find(src);

    if (vertex_it == std::end(edge_map)) {
        return false;
    }

    return vertex_it.second.count(dst)>0;
};

template<typename LABEL>
inline const LABEL& DiGraph<LABEL>::get_edge_label(const VertexId& src, const VertexId& dst) const
{
    return edge_map.at(src).at(dst);
}

template<typename LABEL>
inline const typename DiGraph<LABEL>::DstEdgeMap& DiGraph<LABEL>::get_outgoing_edge_map(const VertexId& src) const
{
    return edge_map.at(src);
}

}; // namespace Races end

#endif // __RACES_DIGRAPH__