/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | www.cfmesh.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    labelListList2D is a List of independently-growable rows
    (List<DynamicList<label>>), used as a direct OpenFOAM-native
    replacement for VRWGraph at call sites that build their graph
    incrementally (row count known upfront, per-row contents discovered
    and appended one element at a time). The free functions below cover
    the handful of VRWGraph member functions that have no direct
    List/DynamicList equivalent.

SourceFiles
    None (header-only)

\*---------------------------------------------------------------------------*/

#ifndef labelListListOps_H
#define labelListListOps_H

#include "DynamicList.H"
#include "List.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{

typedef List<DynamicList<label>> labelListList2D;

#ifndef forAllRow
//- Works for VRWGraph (via its graphRow/graphConstRow view, whose size()
//- delegates to sizeOfRow()) and for List<DynamicList<label>> rows alike,
//- so callers do not need to change as VRWGraph instances are migrated.
//- Defined here (not just in VRWGraph.H) so files that only need
//- labelListList2D infrastructure, with no other VRWGraph usage, don't
//- have to separately include VRWGraph.H just for this macro.
#define forAllRow(graph, rowI, index)                                         \
    for (Foam::label index = 0; index < (graph)[rowI].size(); ++index)
#endif

//- Resize the outer list and set each row's size (mirrors
//- VRWGraph::setSizeAndRowSize / VRWGraphSMPModifier::setSizeAndRowSize)
template<class ListType>
inline void setSizeAndRowSize(labelListList2D& g, const ListType& rowSizes)
{
    const label nRows = rowSizes.size();
    g.setSize(nRows);

    for (label i = 0; i < nRows; ++i)
    {
        g[i].setSize(rowSizes[i]);
    }
}

//- Check if val is present in row rowI (mirrors VRWGraph::found/find)
inline bool found(const labelListList2D& g, const label rowI, const label val)
{
    const DynamicList<label>& row = g[rowI];
    forAll(row, i)
    {
        if (row[i] == val)
        {
            return true;
        }
    }
    return false;
}

//- Return the index within row rowI where val is found, or -1 if not
//- present (mirrors VRWGraph::find)
inline label find(const labelListList2D& g, const label rowI, const label val)
{
    const DynamicList<label>& row = g[rowI];
    forAll(row, i)
    {
        if (row[i] == val)
        {
            return i;
        }
    }
    return -1;
}

//- Set row rowI's contents to l, resizing the row first (mirrors
//- VRWGraph::setRow)
template<class ListType>
inline void setRow(labelListList2D& g, const label rowI, const ListType& l)
{
    g[rowI].setSize(l.size());
    forAll(l, i)
    {
        g[rowI][i] = l[i];
    }
}

//- Append a new row to the graph, containing a copy of l (mirrors
//- VRWGraph::appendList). WARNING: labelListList2D is a plain List, not
//- a DynamicList - g.setSize(g.size()+1) reallocates and copies the
//- *entire* outer list on every call, unlike VRWGraph's original
//- LongList-backed storage which only ever pages (O(1) amortized).
//- Calling this in a loop that grows the graph to N rows is O(N^2) -
//- only safe when the total number of appends stays small/bounded
//- (e.g. partTetMesh's Gauss-Seidel colour-ordering, a handful of
//- iterations). For incremental row-by-row construction up to
//- mesh-element scale (e.g. rebuilding a per-face graph over the whole
//- mesh), use the DynamicList<DynamicList<label>> overload below as the
//- local builder instead, and assign/transfer into the labelListList2D
//- member only once, at the end.
template<class ListType>
inline void appendList(labelListList2D& g, const ListType& l)
{
    g.setSize(g.size() + 1);
    setRow(g, g.size() - 1, l);
}

//- As above, but for a DynamicList<DynamicList<label>> outer container -
//- O(1) amortized per append, since DynamicList (unlike List) grows its
//- capacity geometrically instead of reallocating to the exact size
//- every time. Use this as the local incremental-build type when the
//- final row count scales with mesh size; the result can still be
//- indexed/assigned into a labelListList2D the same way since both are
//- lists of DynamicList<label> rows.
template<class ListType>
inline void appendList(DynamicList<DynamicList<label>>& g, const ListType& l)
{
    // Note: DynamicList(const label) allocates that much *capacity* but
    // sets the logical size back to 0 - must use setSize() (not the
    // constructor) to actually grow the row to l.size() elements before
    // writing into it.
    DynamicList<label> row;
    row.setSize(l.size());
    forAll(l, i)
    {
        row[i] = l[i];
    }
    g.append(row);
}

//- As setRow above, but for a DynamicList<DynamicList<label>> outer
//- container - needed for members like refineBoundaryLayers::newFaces_
//- that use the DynamicList<DynamicList<label>> overload of appendList
//- above (mesh-scale incremental row-count growth) and also need to
//- rewrite an EXISTING row's contents in place (e.g. reversing a face's
//- point order) rather than only ever appending new rows.
template<class ListType>
inline void setRow
(
    DynamicList<DynamicList<label>>& g,
    const label rowI,
    const ListType& l
)
{
    g[rowI].setSize(l.size());
    forAll(l, i)
    {
        g[rowI][i] = l[i];
    }
}

//- As found above, but for a DynamicList<DynamicList<label>> outer
//- container - see setRow's doc comment above for why this overload
//- exists alongside the labelListList2D one.
inline bool found
(
    const DynamicList<DynamicList<label>>& g,
    const label rowI,
    const label val
)
{
    const DynamicList<label>& row = g[rowI];
    forAll(row, i)
    {
        if (row[i] == val)
        {
            return true;
        }
    }
    return false;
}

//- Append val to row if not already present (mirrors VRWGraph::appendUniq)
inline void appendUniq(DynamicList<label>& row, const label val)
{
    forAll(row, i)
    {
        if (row[i] == val)
        {
            return;
        }
    }

    row.append(val);
}

//- Set result to the reverse addressing of origGraph, i.e. result[el]
//- lists all rowI for which el appears in origGraph[rowI]. nRows must be
//- at least (max entry in origGraph) + 1. Serial equivalent of
//- VRWGraph::reverseAddressing(nRows, origGraph) /
//- VRWGraphSMPModifier::reverseAddressing(origGraph) - same two-pass
//- count-then-fill algorithm, so insertion order (increasing rowI within
//- each result row) is identical to the original.
//- GraphType only needs size() and an operator[] returning something
//- with size()/operator[] per row - satisfied by both VRWGraph (via its
//- graphRow/graphConstRow proxy) and labelListList2D itself.
template<class GraphType>
inline void reverseAddressing
(
    labelListList2D& result,
    const label nRows,
    const GraphType& origGraph
)
{
    const label origSize = origGraph.size();
    labelList nElmtsInRow(nRows, 0);

    for (label rowI = 0; rowI < origSize; ++rowI)
    {
        const label rowSize = origGraph[rowI].size();
        for (label i = 0; i < rowSize; ++i)
        {
            ++nElmtsInRow[origGraph[rowI][i]];
        }
    }

    setSizeAndRowSize(result, nElmtsInRow);
    nElmtsInRow = 0;

    for (label rowI = 0; rowI < origSize; ++rowI)
    {
        const label rowSize = origGraph[rowI].size();
        for (label i = 0; i < rowSize; ++i)
        {
            const label el = origGraph[rowI][i];
            result[el][nElmtsInRow[el]++] = rowI;
        }
    }
}

//- As above, but nRows is derived from the data itself (one more than the
//- largest entry in origGraph) rather than supplied by the caller -
//- mirrors VRWGraph::reverseAddressing(origGraph) /
//- VRWGraphSMPModifier::reverseAddressing(origGraph)'s no-nRows overload.
template<class GraphType>
inline void reverseAddressing
(
    labelListList2D& result,
    const GraphType& origGraph
)
{
    const label origSize = origGraph.size();
    label maxValue = -1;

    for (label rowI = 0; rowI < origSize; ++rowI)
    {
        const label rowSize = origGraph[rowI].size();
        for (label i = 0; i < rowSize; ++i)
        {
            maxValue = Foam::max(maxValue, origGraph[rowI][i]);
        }
    }

    reverseAddressing(result, maxValue + 1, origGraph);
}

//- As above, but each entry read from origGraph is first passed through
//- mapper before being used as the target row index - mirrors
//- VRWGraphSMPModifier::reverseAddressing(mapper, origGraph). Used when
//- origGraph's entries are indices into a different numbering space than
//- the desired result rows (e.g. global point indices that need mapping
//- to boundary-point indices).
template<class ListType, class GraphType>
inline void reverseAddressing
(
    labelListList2D& result,
    const label nRows,
    const ListType& mapper,
    const GraphType& origGraph
)
{
    const label origSize = origGraph.size();
    labelList nElmtsInRow(nRows, 0);

    for (label rowI = 0; rowI < origSize; ++rowI)
    {
        const label rowSize = origGraph[rowI].size();
        for (label i = 0; i < rowSize; ++i)
        {
            ++nElmtsInRow[mapper[origGraph[rowI][i]]];
        }
    }

    setSizeAndRowSize(result, nElmtsInRow);
    nElmtsInRow = 0;

    for (label rowI = 0; rowI < origSize; ++rowI)
    {
        const label rowSize = origGraph[rowI].size();
        for (label i = 0; i < rowSize; ++i)
        {
            const label el = mapper[origGraph[rowI][i]];
            result[el][nElmtsInRow[el]++] = rowI;
        }
    }
}

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
