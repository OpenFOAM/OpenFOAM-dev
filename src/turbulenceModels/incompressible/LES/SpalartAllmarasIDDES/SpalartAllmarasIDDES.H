/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::incompressible::LESModels::SpalartAllmarasIDDES

Group
    grpIcoDESTurbulence

Description
    SpalartAllmarasIDDES LES turbulence model for incompressible flows

SourceFiles
    SpalartAllmarasIDDES.C

\*---------------------------------------------------------------------------*/

#ifndef SpalartAllmarasIDDES_H
#define SpalartAllmarasIDDES_H

#include "SpalartAllmaras.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                        Class SpalartAllmarasIDDES Declaration
\*---------------------------------------------------------------------------*/

class SpalartAllmarasIDDES
:
    public SpalartAllmaras
{
    // Private data

        // Model constants

            autoPtr<Foam::LESdelta> hmax_;
            autoPtr<Foam::LESdelta> IDDESDelta_;
            dimensionedScalar fwStar_;
            dimensionedScalar cl_;
            dimensionedScalar ct_;


    // Private Member Functions

        tmp<volScalarField> alpha() const;
        tmp<volScalarField> ft(const volScalarField& S) const;
        tmp<volScalarField> fl(const volScalarField& S) const;

        tmp<volScalarField> rd
        (
            const volScalarField& visc,
            const volScalarField& S
        ) const;

        //- Delay function
        tmp<volScalarField> fd(const volScalarField& S) const;

        // Disallow default bitwise copy construct and assignment
        SpalartAllmarasIDDES(const SpalartAllmarasIDDES&);
        SpalartAllmarasIDDES& operator=(const SpalartAllmarasIDDES&);


protected:

    // Protected Member Functions

        //- Length scale
        virtual tmp<volScalarField> dTilda(const volScalarField& S) const;


public:

    //- Runtime type information
    TypeName("SpalartAllmarasIDDES");


    // Constructors

        //- Construct from components
        SpalartAllmarasIDDES
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~SpalartAllmarasIDDES()
    {}


    // Member Functions

        //- Access function to filter width
        virtual const volScalarField& delta() const
        {
            return IDDESDelta_();
        }

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
