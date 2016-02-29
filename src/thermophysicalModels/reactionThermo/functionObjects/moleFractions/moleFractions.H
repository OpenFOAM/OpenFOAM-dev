/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

SeeAlso
    Foam::functionObject
    Foam::OutputFilterFunctionObject

SourceFiles
    moleFractions.C

\*---------------------------------------------------------------------------*/

#ifndef moleFractions_H
#define moleFractions_H

#include "functionObjectFile.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class moleFractions Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class moleFractions
:
    public functionObjectFile
{
    // Private data

        //- Name of moleFractions functionObject
        word name_;

        //- Reference to the mesh
        const fvMesh& mesh_;

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

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        moleFractions
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~moleFractions();


    // Member Functions

        //- Return name of the moleFractions functionObject
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the moleFractions data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the moleFractions and write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
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
