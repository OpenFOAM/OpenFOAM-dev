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
    Foam::ensightPart

Description
    Base class for ensightPartCells and ensightPartFaces

SourceFiles
    ensightPart.C
    ensightPartIO.C
    ensightPartTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef ensightPart_H
#define ensightPart_H

#include "ensightFile.H"
#include "ensightGeoFile.H"
#include "typeInfo.H"
#include "labelList.H"
#include "polyMesh.H"
#include "Field.H"
#include "IOPtrList.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class ensightPart;

Ostream& operator<<(Ostream&, const ensightPart&);
ensightGeoFile& operator<<(ensightGeoFile&, const ensightPart&);


/*---------------------------------------------------------------------------*\
                         Class ensightPart Declaration
\*---------------------------------------------------------------------------*/

class ensightPart
{
    // Private data

        // Static data members
        static const List<word> elemTypes_;


protected:

    // Protected data

        //- Part number
        label number_;

        //- Part name (or description)
        string name_;

        //- Simple labelList with a name
        labelListList elemLists_;

        //- Start offset for elemLists_
        label offset_;

        //- Number of elements in this part
        label size_;

        //- Cell or face data
        bool isCellData_;

        //- Material id (numeric)
        label matId_;

        //- pointField referenced
        const pointField& points_;


    // Protected Classes

        //- Track the points used by the part and map global to local indices
        class localPoints
        {
        public:
            //- Number of points used
            label nPoints;

            //- Map global to local indices
            labelList list;

            //- Null constructor
            localPoints()
            :
                nPoints(0),
                list(0)
            {}

            //- Construct for mesh points
            localPoints(const pointField& pts)
            :
                nPoints(0),
                list(pts.size(), -1)
            {}
        };


    // Protected Member Functions

        //- Reconstruct part characteristics (eg, element types) from Istream
        //  A part reconstructed in this manner can be used when writing fields,
        //  but cannot be used to write a new geometry
        void reconstruct(Istream&);

        //- Check for fully defined fields
        bool isFieldDefined(const List<scalar>&) const;

        //- Write the part header
        void writeHeader(ensightFile&, bool withDescription=false) const;

        //- Write a scalar field for idList
        //  A null reference for idList writes the perNode values
        void writeFieldList
        (
            ensightFile& os,
            const List<scalar>& field,
            const labelUList& idList
        ) const;

        //- Track points used
        virtual localPoints calcLocalPoints() const
        {
            return localPoints();
        }

        //- Write connectivities
        virtual void writeConnectivity
        (
            ensightGeoFile&,
            const word& key,
            const labelUList& idList,
            const labelUList& pointMap
        ) const
        {}


public:

    //- Runtime type information
    TypeName("ensightPart");


    // Constructors

        //- Construct null
        ensightPart();

        //- Construct empty part with number and description
        ensightPart(label partNumber, const string& partDescription);

        //- Construct part with number, description and points reference
        ensightPart
        (
            label partNumber,
            const string& partDescription,
            const pointField& points
        );

        //- Construct as copy
        ensightPart(const ensightPart&);


    // Selectors

        // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            ensightPart,
            istream,
            (
                Istream& is
            ),
            (is)
        );

        //- Construct and return clone
        autoPtr<ensightPart> clone() const
        {
            return autoPtr<ensightPart>(new ensightPart(*this));
        };

        //- Reconstruct part characteristics on freestore from Istream
        //  \sa reconstruct
        static autoPtr<ensightPart> New(Istream&);


    //- Destructor
    virtual ~ensightPart();


    // Static members

        virtual const List<word>& elementTypes() const
        {
            return elemTypes_;
        }


    // Access

        //- Number of elements in this part
        label size() const
        {
            return size_;
        }

        //- Represents cell data
        bool isCellData() const
        {
            return isCellData_;
        }

        //- Represents face data
        bool isFaceData() const
        {
            return !isCellData_;
        }

        //- Part number
        label number() const
        {
            return number_;
        }

        //- Part name or description
        const string& name() const
        {
            return name_;
        }

        //- Material id
        label materialId() const
        {
            return matId_;
        }

        //- non-const access to part name or description
        void name(const string& value)
        {
            name_ = value;
        }

        //- non-const access to material id
        void materialId(const label value)
        {
            matId_ = value;
        }

        //- Simple labelList with a name
        const labelListList& elemLists() const
        {
            return elemLists_;
        }

        //- Offset for element ids
        label offset() const
        {
            return offset_;
        }


    // Edit

        //- Renumber elements
        void renumber(const labelUList&);

        //- Write summary information about the object
        bool writeSummary(Ostream&) const;

        //- Write reconstruction information for the object
        bool writeData(Ostream&) const;

        //- Write geometry
        virtual void writeGeometry(ensightGeoFile&) const
        {}

        //- Helper: write geometry given the pointField
        void writeGeometry(ensightGeoFile&, const pointField&) const;

        //- Write scalar field
        //  optionally write data per node
        void writeScalarField
        (
            ensightFile&,
            const List<scalar>& field,
            const bool perNode = false
        ) const;

        //- Write vector field components
        //  optionally write data per node
        void writeVectorField
        (
            ensightFile&,
            const List<scalar>& field0,
            const List<scalar>& field1,
            const List<scalar>& field2,
            const bool perNode = false
        ) const;


        //- Write generalized field components
        //  optionally write data per node
        template<class Type>
        void writeField
        (
            ensightFile&,
            const Field<Type>&,
            const bool perNode = false
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ensightPart&)
        {
            NotImplemented;
        }


    // IOstream Operators

        //- Write data (reconstruction information)
        friend Ostream& operator<<(Ostream&, const ensightPart&);

        //- Write geometry
        friend ensightGeoFile& operator<<(ensightGeoFile&, const ensightPart&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensightPartTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
