/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::lduMesh

Description
    Abstract base class for meshes which provide LDU addressing for the
    construction of lduMatrix and LDU-solvers.

\*---------------------------------------------------------------------------*/

#ifndef lduMesh_H
#define lduMesh_H

#include "lduAddressing.H"
#include "lduInterfacePtrsList.H"
#include "typeInfo.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class objectRegistry;

/*---------------------------------------------------------------------------*\
                           Class lduMesh Declaration
\*---------------------------------------------------------------------------*/

class lduMesh
{

public:

    //- Runtime type information
    TypeName("lduMesh");


    // Constructors

    //- Destructor
    virtual ~lduMesh()
    {}


    // Member Functions

        // Access

            //- Return the object registry
            virtual const objectRegistry& thisDb() const;

            //- Return ldu addressing
            virtual const lduAddressing& lduAddr() const = 0;

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            virtual lduInterfacePtrsList interfaces() const = 0;

            //- Return communicator used for parallel communication
            virtual int comm() const = 0;

            //- Helper: reduce with current communicator
            template<class T, class BinaryOp>
            void reduce
            (
                T& Value,
                const BinaryOp& bop
            ) const;

        // Info

            //- Return info proxy.
            //  Used to print mesh information to a stream
            InfoProxy<lduMesh> info() const
            {
                return *this;
            }

    // Ostream operator

        friend Ostream& operator<<(Ostream&, const InfoProxy<lduMesh>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "lduMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
