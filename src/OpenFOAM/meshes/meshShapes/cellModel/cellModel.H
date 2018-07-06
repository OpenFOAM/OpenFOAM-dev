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
    Foam::cellModel

Description
    Maps a geometry to a set of cell primitives, which enables
    geometric cell data to be calculated without access to the primitive
    geometric level.  This means mapping a 3D geometry to a set of
    pyramids which are each described by a cell face and the cell centre
    point.

SourceFiles
    cellModelI.H

\*---------------------------------------------------------------------------*/

#ifndef cellModel_H
#define cellModel_H

#include "pointField.H"
#include "edgeList.H"
#include "faceList.H"
#include "InfoProxy.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class cellModel;
inline bool operator==(const cellModel&, const cellModel&);
inline bool operator!=(const cellModel&, const cellModel&);
Ostream& operator<<(Ostream&, const cellModel&);


/*---------------------------------------------------------------------------*\
                           Class cellModel Declaration
\*---------------------------------------------------------------------------*/

class cellModel
{
    // Private data

        //- Name
        word name_;

        //- Label in the model list
        label index_;

        //- Number of points in the model which determines the geometry
        label nPoints_;

        //- Faces of the model
        faceList faces_;

        //- Edges of the model
        edgeList edges_;


public:

    // Constructors

        //- Construct from Istream
        cellModel(Istream&);

        //- Return a new cellModel on free-store created from Istream
        static autoPtr<cellModel> New(Istream& is)
        {
            return autoPtr<cellModel>(new cellModel(is));
        }

        //- Return clone
        autoPtr<cellModel> clone() const
        {
            return autoPtr<cellModel>(new cellModel(*this));
        }


    // Member functions

        // Access

            //- Return model name
            inline const word& name() const;

            //- Return index of model in the model list
            inline label index() const;

            //- Return number of points
            inline label nPoints() const;

            //- Return number of edges
            inline label nEdges() const;

            //- Return number of faces
            inline label nFaces() const;

            //- Return list of edges
            inline edgeList edges(const labelList& pointLabels) const;

            //- Return a raw list of model faces
            inline const faceList& modelFaces() const;

            //- Return list of faces
            inline faceList faces(const labelList& pointLabels) const;


        //- Vector centroid
        vector centre
        (
            const labelList& pointLabels,
            const pointField& points
        ) const;

        //- Cell volume
        scalar mag
        (
            const labelList& pointLabels,
            const pointField& points
        ) const;

        //- Return info proxy.
        //  Used to print token information to a stream
        InfoProxy<cellModel> info() const
        {
            return *this;
        }

        //- WriteData member function required by regIOobject
        bool writeData(Ostream& os) const
        {
            os << *this;
            return os.good();
        }


    // Friend operators

        //- Equality operator: true => ptr to models are equal !
        friend bool operator==(const cellModel&, const cellModel&);

        //- Inequality operator: true => ptr to models are not equal !
        friend bool operator!=(const cellModel&, const cellModel&);


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const cellModel&);
};


template<>
Ostream& operator<<(Ostream& os, const InfoProxy<cellModel>& ip);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
