/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::moleFractions

Description
    This function object calculates mole-fraction fields from the mass-fraction
    fields of the psi/rhoReactionThermo and caches them for output and further
    post-processing.

    The names of the mole-fraction fields are obtained from the corresponding
    mass-fraction fields prepended by "X_"

    Example of function object specification:
    \verbatim
    moleFractions
    {
        type psiReactionThermoMoleFractions;
    }
    \endverbatim
    or
    \verbatim
    moleFractions
    {
        type rhoReactionThermoMoleFractions;
    }
    \endverbatim
    depending on the thermodynamics package used in the solver.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    moleFractions.C

\*---------------------------------------------------------------------------*/

#ifndef moleFractions_H
#define moleFractions_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class moleFractions Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class moleFractions
:
    public functionObjects::fvMeshFunctionObject
{
    // Private data

        //- Species mole fractions
        PtrList<volScalarField> X_;


    // Private Member Functions

        //- Calculate the mole fraction fields
        virtual void calculateMoleFractions();

        //- Disallow default bitwise copy construct
        moleFractions(const moleFractions&);

        //- Disallow default bitwise assignment
        void operator=(const moleFractions&);


public:

    //- Runtime type information
    TypeName("moleFractions");


    // Constructors

        //- Construct from Time and dictionary
        moleFractions
        (
            const word& name,
            const Time& t,
            const dictionary& dict
        );


    //- Destructor
    virtual ~moleFractions();


    // Member Functions

        //- Read the moleFractions data
        virtual bool read(const dictionary&);

        //- Calculate the mole-fraction fields
        virtual bool execute();

        //- The mole-fraction fields auto-write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "moleFractions.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
