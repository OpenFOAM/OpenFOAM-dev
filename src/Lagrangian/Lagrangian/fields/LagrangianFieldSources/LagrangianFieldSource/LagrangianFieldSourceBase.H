/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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
    Foam::LagrangianFieldSourceBase

Description
    Base class for Lagrangian source conditions

SourceFiles
    LagrangianFieldSourceBase.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianFieldSourceBase_H
#define LagrangianFieldSourceBase_H

#include "LagrangianFieldsFwd.H"
#include "LagrangianSubFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class regIOobject;
class dimensionSet;
class objectRegistry;
class dictionary;
class LagrangianSubMesh;
class LagrangianModel;
class LagrangianSource;
class LagrangianInjection;

/*---------------------------------------------------------------------------*\
                    Class LagrangianFieldSourceBase Declaration
\*---------------------------------------------------------------------------*/

class LagrangianFieldSourceBase
{
protected:

    // Protected Data

        //- Reference to the internal IO object
        const regIOobject& internalIo_;


public:

    //- Runtime type information
    TypeName("LagrangianFieldSourceBase");


    // Constructors

        //- Construct from internal field
        LagrangianFieldSourceBase(const regIOobject&);


    //- Destructor
    virtual ~LagrangianFieldSourceBase();


    // Member Functions

        // Access

            //- Return local objectRegistry
            const objectRegistry& db() const;

            //- Return internal field name
            const word& internalName() const;

            //- Return internal field group
            word internalGroup() const;

            //- Return internal dimensions reference
            virtual const dimensionSet& internalDimensions() const = 0;


        // Evaluation

            //- Cast a to the specified type. Handle errors.
            template<class OtherFieldSourceType>
            const OtherFieldSourceType& fieldSourceCast
            (
                const LagrangianModel&
            ) const;

            //- Lookup and return another field source
            template<class OtherType>
            const LagrangianFieldSource<OtherType>& fieldSource
            (
                const word& name,
                const LagrangianModel&
            ) const;

            //- Lookup and return another field source and then cast to a
            //  specified type. Handle errors.
            template<class OtherType, class OtherFieldSourceType>
            const OtherFieldSourceType& fieldSourceCast
            (
                const word& name,
                const LagrangianModel&
            ) const;

            //- Cast a model to the specified type. Handle errors.
            template<class OtherModelType>
            const OtherModelType& modelCast
            (
                const LagrangianModel&
            ) const;

            //- Lookup and return the value of another field source
            template<class OtherType>
            tmp<LagrangianSubField<OtherType>> value
            (
                const word& name,
                const LagrangianSource&,
                const LagrangianSubMesh&
            ) const;

            //- Lookup and return the value of another field source
            template<class OtherType>
            tmp<LagrangianSubField<OtherType>> value
            (
                const word& name,
                const LagrangianInjection&,
                const LagrangianSubMesh&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LagrangianFieldSourceBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
