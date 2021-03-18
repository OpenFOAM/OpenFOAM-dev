/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::fv::heatTransfer

Description
    Model for heat exchange. Requires specification of an ambient temperature
    with which to exchange heat, and a model for the heat transfer coefficient
    (htc) and the area per unit volume (AoV). These are then used to apply the
    following source to the energy equation:

    \f[
        -htc*AoV*(T_a - T)
    \f]

    If the semiImplicit option is set, then this becomes:

    \f[
        -htc*AoV*(T_a - T) + htc*AoV/Cp*h - Sp(htc*AoV/Cp, h);
    \f]

Usage
    Example usage:
    \verbatim
    heatTransfer
    {
        type            heatTransfer;

        heatTransferCoeffs
        {
            selectionMode   cellSet;
            cellSet         c0;

            semiImplicit    no;

            Ta              300;

            type            constant;

            AoV             200;
            htc             10;
        }
    }
    \endverbatim

See also
    fv::heatTransferModel

SourceFiles
    heatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransfer_H
#define heatTransfer_H

#include "fvModel.H"
#include "fvCellSet.H"
#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                Class heatTransfer Declaration
\*---------------------------------------------------------------------------*/

class heatTransfer
:
    public fvModel
{
    // Private data

        //- The set of cells the model applies to
        fvCellSet set_;

        //- Flag to activate semi-implicit coupling
        bool semiImplicit_;

        //- Name of temperature field; default = "T"
        word TName_;

        //- Ambient temperature
        dimensionedScalar Ta_;

        //- The heat transfer model
        autoPtr<heatTransferModel> heatTransferModel_;


    // Private member functions

        //- Non-virtual read
        void readCoeffs();


public:

    //- Runtime type information
    TypeName("heatTransfer");


    // Constructors

        //- Construct from dictionary
        heatTransfer
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~heatTransfer();


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Sources

            //- Source term to energy equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;

            //- Source term to compressible energy equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;


        // Correction

            //- Correct the model
            virtual void correct();


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
