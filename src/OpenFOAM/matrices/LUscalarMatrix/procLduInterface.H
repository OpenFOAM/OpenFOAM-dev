/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::procLduInterface

Description
    IO interface for processorLduInterface

SourceFiles
    procLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef procLduInterface_H
#define procLduInterface_H

#include "labelList.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class lduInterfaceField;


// Forward declaration of friend functions and operators

class procLduInterface;

Ostream& operator<<(Ostream&, const procLduInterface&);


/*---------------------------------------------------------------------------*\
                           Class procLduInterface Declaration
\*---------------------------------------------------------------------------*/

class procLduInterface
{
    // Private Data

        labelList faceCells_;
        scalarField coeffs_;
        label myProcNo_;
        label neighbProcNo_;
        label tag_;
        label comm_;


    // Private Member Functions

        //- Dissallow construction as copy
        procLduInterface(const procLduInterface&);


public:

    friend class LUscalarMatrix;


    // Constructors

        procLduInterface
        (
            const lduInterfaceField& interface,
            const scalarField& coeffs
        );

        procLduInterface(Istream& is);

        autoPtr<procLduInterface> clone()
        {
            NotImplemented;
            return autoPtr<procLduInterface>(nullptr);
        }

        static autoPtr<procLduInterface> New(Istream& is)
        {
            return autoPtr<procLduInterface>(new procLduInterface(is));
        }


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const procLduInterface&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
