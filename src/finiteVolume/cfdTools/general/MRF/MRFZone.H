/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::MRFZone

Description
    MRF zone definition based on cell zone and parameters
    obtained from a control dictionary constructed from the given stream.

    The rotation of the MRF region is defined by an origin and axis of
    rotation and an angular speed.

SourceFiles
    MRFZone.C

\*---------------------------------------------------------------------------*/

#ifndef MRFZone_H
#define MRFZone_H

#include "fvCellSet.H"
#include "volFieldsFwd.H"
#include "surfaceFields.H"
#include "omega.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class MRFZone Declaration
\*---------------------------------------------------------------------------*/

class MRFZone
{
    // Private Data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Name of the MRF region
        const word name_;

        //- Coefficients dictionary
        dictionary coeffs_;

        //- MRF zone cell set
        fvCellSet cellSet_;

        //- Internal faces that are in the MRF region
        labelList internalFaces_;

        //- Patch faces that are in the MRF region
        labelListList patchFaces_;

        //- Origin of the axis
        const vector origin_;

        //- Axis vector
        vector axis_;

        //- Angular velocity (rad/sec or rpm)
        Function1s::omega omega_;


    // Private Member Functions

        //- Divide faces in frame according to patch
        void setMRFFaces();

        //- Check that the case has been updated with correct MRF BCs
        void checkMRFBCs(const volVectorField& U) const;

        //- Make the given absolute mass/vol flux relative within the MRF region
        template<class RhoFieldType>
        void makeRelativeRhoFlux
        (
            const RhoFieldType& rho,
            surfaceScalarField& phi
        ) const;

        //- Make the given absolute mass/vol flux relative within the MRF region
        template<class RhoFieldType>
        void makeRelativeRhoFlux
        (
            const RhoFieldType& rho,
            FieldField<fvsPatchField, scalar>& phi
        ) const;

        //- Make the given absolute mass/vol flux relative within the MRF region
        template<class RhoFieldType>
        void makeRelativeRhoFlux
        (
            const RhoFieldType& rho,
            Field<scalar>& phi,
            const label patchi
        ) const;

        //- Make the given relative mass/vol flux absolute within the MRF region
        template<class RhoFieldType>
        void makeAbsoluteRhoFlux
        (
            const RhoFieldType& rho,
            surfaceScalarField& phi
        ) const;


public:

    // Declare name of the class and its debug switch
    ClassName("MRFZone");


    // Constructors

        //- Construct from fvMesh
        MRFZone
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        MRFZone(const MRFZone&) = delete;

        //- Return clone
        autoPtr<MRFZone> clone() const
        {
            NotImplemented;
            return autoPtr<MRFZone>(nullptr);
        }


    // Member Functions

        //- Return const access to the MRF region name
        inline const word& name() const;

        //- Return the current Omega vector
        vector Omega() const;

        //- Add the Coriolis force contribution to the acceleration field
        void addCoriolis
        (
            const volVectorField& U,
            volVectorField& ddtU
        ) const;

        //- Add the centrifugal acceleration
        void addCentrifugalAcceleration
        (
            volVectorField& centrifugalAcceleration
        ) const;

        //- Make the given absolute velocity relative within the MRF region
        void makeRelative(volVectorField& U) const;

        //- Make the given absolute flux relative within the MRF region
        void makeRelative(surfaceScalarField& phi) const;

        //- Make the given absolute boundary flux relative
        //  within the MRF region
        void makeRelative(FieldField<fvsPatchField, scalar>& phi) const;

        //- Make the given absolute patch flux relative
        //  within the MRF region
        void makeRelative(Field<scalar>& phi, const label patchi) const;

        //- Make the given absolute mass-flux relative within the MRF region
        void makeRelative
        (
            const surfaceScalarField& rho,
            surfaceScalarField& phi
        ) const;

        //- Make the given relative patch velocity relative
        //  within the MRF region
        void makeRelative(Field<vector>& Up, const label patchi) const;

        //- Make the given relative velocity absolute within the MRF region
        void makeAbsolute(volVectorField& U) const;

        //- Make the given relative flux absolute within the MRF region
        void makeAbsolute(surfaceScalarField& phi) const;

        //- Make the given relative mass-flux absolute within the MRF region
        void makeAbsolute
        (
            const surfaceScalarField& rho,
            surfaceScalarField& phi
        ) const;

        //- Make the given relative patch velocity absolute
        //  within the MRF region
        void makeAbsolute(Field<vector>& Up, const label patchi) const;

        //- Zero the MRF region of the given field
        template<class Type>
        void zero(SurfaceField<Type>& phi) const;

        //- Update MRFZone faces if the mesh topology changes
        void update();


    // I-O

        //- Read MRF dictionary
        bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const MRFZone&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MRFZoneTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MRFZoneI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
