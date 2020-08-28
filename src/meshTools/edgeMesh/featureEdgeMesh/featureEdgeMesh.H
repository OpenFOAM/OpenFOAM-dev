/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::featureEdgeMesh

Description
    edgeMesh + IO.

See also
    Foam::extendedFeatureEdgeMesh

SourceFiles
    featureEdgeMesh.C

\*---------------------------------------------------------------------------*/

#ifndef featureEdgeMesh_H
#define featureEdgeMesh_H

#include "edgeMesh.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class featureEdgeMesh Declaration
\*---------------------------------------------------------------------------*/

class featureEdgeMesh
:
    public regIOobject,
    public edgeMesh
{

public:

    TypeName("featureEdgeMesh");


    // Constructors

        //- Construct (read) given an IOobject
        featureEdgeMesh(const IOobject&);

        //- Construct from featureEdgeMesh data
        featureEdgeMesh
        (
            const IOobject&,
            const pointField&,
            const edgeList&
        );

        //- Copy constructor
        featureEdgeMesh(const IOobject&, const featureEdgeMesh&);


    // IO

        //- ReadData function required for regIOobject read operation
        virtual bool readData(Istream&);

        //- WriteData function required for regIOobject write operation
        virtual bool writeData(Ostream&) const;

        //- Is object global
        virtual bool global() const
        {
            return true;
        }

        //- Return complete path + object name if the file exists
        //  either in the case/processor or case otherwise null
        virtual fileName filePath() const
        {
            return globalFilePath(type());
        }
};


//- Template function for obtaining global status
template<>
inline bool typeGlobal<featureEdgeMesh>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
