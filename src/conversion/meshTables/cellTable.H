/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::cellTable

Description
    The cellTable persistent data saved as a Map<dictionary>.

    The meshReader supports cellTable information.

    The <tt>constant/cellTable</tt> file is an \c IOMap<dictionary> that is
    used to save the information persistently. It contains the cellTable
    information of the following form:

    \verbatim
        (
            ID
            {
                Label           WORD;
                MaterialType    WORD;
                MaterialId      INT;
                PorosityId      INT;
                ColorIdx        INT;
                ...
            }
        ...
        )
    \endverbatim

    If the \a Label is missing, a value <tt>cellTable_{ID}</tt> will be
    inferred. If the \a MaterialType is missing, the value @a fluid will
    be inferred.

SourceFiles
    cellTable.C

\*---------------------------------------------------------------------------*/

#ifndef cellTable_H
#define cellTable_H

#include "polyMesh.H"
#include "Map.H"
#include "dictionary.H"
#include "labelList.H"
#include "wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class cellTable Declaration
\*---------------------------------------------------------------------------*/

class cellTable
:
    public Map<dictionary>
{
    // Private data

       static const char* const defaultMaterial_;


    // Private Member Functions

        //- Map from cellTable ID => zone number
        Map<label> zoneMap() const;

        //- A contiguous list of cellTable names
        List<word> namesList() const;

        //- Add required entries - MaterialType
        void addDefaults();

        void setEntry(const label id, const word& keyWord, const word& value);

        //- Disallow default bitwise copy construct
        cellTable(const cellTable&);


public:

    // Constructors

        //- Construct null
        cellTable();

        //- Construct read from registry, name. instance
        cellTable
        (
            const objectRegistry&,
            const word& name = "cellTable",
            const fileName& instance = "constant"
        );


    //- Destructor
    ~cellTable();


    // Member Functions

        //- Append to the end, return index
        label append(const dictionary&);

        //- Return index corresponding to name
        //  returns -1 if not found
        label findIndex(const word& name) const;

        //- Return the name corresponding to id
        //  returns cellTable_ID if not otherwise defined
        word name(const label id) const;

        //- Return a Map of (id => name)
        Map<word> names() const;

        //- Return a Map of (id => names) selected by patterns
        Map<word> names(const UList<wordRe>& patterns) const;

        //- Return a Map of (id => name) for materialType
        //  (fluid | solid | shell)
        Map<word> selectType(const word& materialType) const;

        //- Return a Map of (id => name) for fluids
        Map<word> fluids() const;

        //- Return a Map of (id => name) for shells
        Map<word> shells() const;

        //- Return a Map of (id => name) for solids
        Map<word> solids() const;

        //- Return a Map of (id => fluid|solid|shell)
        Map<word> materialTypes() const;

        //- Assign material Type
        void setMaterial(const label, const word&);

        //- Assign name
        void setName(const label, const word&);

        //- Assign default name if not already set
        void setName(const label);

        //- Read constant/cellTable
        void readDict
        (
            const objectRegistry&,
            const word& name = "cellTable",
            const fileName& instance = "constant"
        );

        //- Write constant/cellTable for later reuse
        void writeDict
        (
            const objectRegistry&,
            const word& name = "cellTable",
            const fileName& instance = "constant"
        ) const;


    // Member Operators

        //- Assignment
        void operator=(const cellTable&);

        //- Assign from Map<dictionary>
        void operator=(const Map<dictionary>&);

        //- Assign from cellZones
        void operator=(const polyMesh&);


    // Friend Functions

        //- Classify tableIds into cellZones according to the cellTable
        void addCellZones(polyMesh&, const labelList& tableIds) const;

        //- Combine tableIds together
        //  each dictionary entry is a wordList
        void combine(const dictionary&, labelList& tableIds);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
