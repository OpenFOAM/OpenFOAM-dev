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

Application
    surfaceInertia

Description
    Calculates the inertia tensor, principal axes and moments of a
    command line specified triSurface. Inertia can either be of the
    solid body or of a thin shell.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "ListOps.H"
#include "triSurface.H"
#include "OFstream.H"
#include "meshTools.H"
#include "Random.H"
#include "transform.H"
#include "IOmanip.H"
#include "Pair.H"
#include "momentOfInertia.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    argList::addNote
    (
        "Calculates the inertia tensor and principal axes and moments "
        "of the specified surface.\n"
        "Inertia can either be of the solid body or of a thin shell."
    );

    argList::validArgs.append("surface file");
    argList::addBoolOption
    (
        "shellProperties",
        "inertia of a thin shell"
    );

    argList::addOption
    (
        "density",
        "scalar",
        "Specify density, "
        "kg/m3 for solid properties, kg/m2 for shell properties"
    );

    argList::addOption
    (
        "referencePoint",
        "vector",
        "Inertia relative to this point, not the centre of mass"
    );

    argList args(argc, argv);

    const fileName surfFileName = args[1];
    const scalar density = args.optionLookupOrDefault("density", 1.0);

    vector refPt = Zero;
    bool calcAroundRefPt = args.optionReadIfPresent("referencePoint", refPt);

    triSurface surf(surfFileName);

    scalar m = 0.0;
    vector cM = Zero;
    tensor J = Zero;

    if (args.optionFound("shellProperties"))
    {
        momentOfInertia::massPropertiesShell(surf, density, m, cM, J);
    }
    else
    {
        momentOfInertia::massPropertiesSolid(surf, density, m, cM, J);
    }

    if (m < 0)
    {
        WarningInFunction
            << "Negative mass detected, the surface may be inside-out." << endl;
    }

    vector eVal = eigenValues(J);

    tensor eVec = eigenVectors(J);

    label pertI = 0;

    Random rand(57373);

    while ((magSqr(eVal) < vSmall) && pertI < 10)
    {
        WarningInFunction
            << "No eigenValues found, shape may have symmetry, "
            << "perturbing inertia tensor diagonal" << endl;

        J.xx() *= 1.0 + small*rand.scalar01();
        J.yy() *= 1.0 + small*rand.scalar01();
        J.zz() *= 1.0 + small*rand.scalar01();

        eVal = eigenValues(J);

        eVec = eigenVectors(J);

        pertI++;
    }

    bool showTransform = true;

    if
    (
        (mag(eVec.x() ^ eVec.y()) > (1.0 - small))
     && (mag(eVec.y() ^ eVec.z()) > (1.0 - small))
     && (mag(eVec.z() ^ eVec.x()) > (1.0 - small))
    )
    {
        // Make the eigenvectors a right handed orthogonal triplet
        eVec = tensor
        (
            eVec.x(),
            eVec.y(),
            eVec.z() * sign((eVec.x() ^ eVec.y()) & eVec.z())
        );

        // Finding the most natural transformation.  Using Lists
        // rather than tensors to allow indexed permutation.

        // Cartesian basis vectors - right handed orthogonal triplet
        List<vector> cartesian(3);

        cartesian[0] = vector(1, 0, 0);
        cartesian[1] = vector(0, 1, 0);
        cartesian[2] = vector(0, 0, 1);

        // Principal axis basis vectors - right handed orthogonal
        // triplet
        List<vector> principal(3);

        principal[0] = eVec.x();
        principal[1] = eVec.y();
        principal[2] = eVec.z();

        scalar maxMagDotProduct = -great;

        // Matching axis indices, first: cartesian, second:principal

        Pair<label> match(-1, -1);

        forAll(cartesian, cI)
        {
            forAll(principal, pI)
            {
                scalar magDotProduct = mag(cartesian[cI] & principal[pI]);

                if (magDotProduct > maxMagDotProduct)
                {
                    maxMagDotProduct = magDotProduct;

                    match.first() = cI;

                    match.second() = pI;
                }
            }
        }

        scalar sense = sign
        (
            cartesian[match.first()] & principal[match.second()]
        );

        if (sense < 0)
        {
            // Invert the best match direction and swap the order of
            // the other two vectors

            List<vector> tPrincipal = principal;

            tPrincipal[match.second()] *= -1;

            tPrincipal[(match.second() + 1) % 3] =
                principal[(match.second() + 2) % 3];

            tPrincipal[(match.second() + 2) % 3] =
                principal[(match.second() + 1) % 3];

            principal = tPrincipal;

            vector tEVal = eVal;

            tEVal[(match.second() + 1) % 3] = eVal[(match.second() + 2) % 3];

            tEVal[(match.second() + 2) % 3] = eVal[(match.second() + 1) % 3];

            eVal = tEVal;
        }

        label permutationDelta = match.second() - match.first();

        if (permutationDelta != 0)
        {
            // Add 3 to the permutationDelta to avoid negative indices

            permutationDelta += 3;

            List<vector> tPrincipal = principal;

            vector tEVal = eVal;

            for (label i = 0; i < 3; i++)
            {
                tPrincipal[i] = principal[(i + permutationDelta) % 3];

                tEVal[i] = eVal[(i + permutationDelta) % 3];
            }

            principal = tPrincipal;

            eVal = tEVal;
        }

        label matchedAlready = match.first();

        match =Pair<label>(-1, -1);

        maxMagDotProduct = -great;

        forAll(cartesian, cI)
        {
            if (cI == matchedAlready)
            {
                continue;
            }

            forAll(principal, pI)
            {
                if (pI == matchedAlready)
                {
                    continue;
                }

                scalar magDotProduct = mag(cartesian[cI] & principal[pI]);

                if (magDotProduct > maxMagDotProduct)
                {
                    maxMagDotProduct = magDotProduct;

                    match.first() = cI;

                    match.second() = pI;
                }
            }
        }

        sense = sign
        (
            cartesian[match.first()] & principal[match.second()]
        );

        if (sense < 0 || (match.second() - match.first()) != 0)
        {
            principal[match.second()] *= -1;

            List<vector> tPrincipal = principal;

            tPrincipal[(matchedAlready + 1) % 3] =
                principal[(matchedAlready + 2) % 3]*-sense;

            tPrincipal[(matchedAlready + 2) % 3] =
                principal[(matchedAlready + 1) % 3]*-sense;

            principal = tPrincipal;

            vector tEVal = eVal;

            tEVal[(matchedAlready + 1) % 3] = eVal[(matchedAlready + 2) % 3];

            tEVal[(matchedAlready + 2) % 3] = eVal[(matchedAlready + 1) % 3];

            eVal = tEVal;
        }

        eVec = tensor(principal[0], principal[1], principal[2]);

        // {
        //     tensor R = rotationTensor(vector(1, 0, 0), eVec.x());

        //     R = rotationTensor(R & vector(0, 1, 0), eVec.y()) & R;

        //     Info<< "R = " << nl << R << endl;

        //     Info<< "R - eVec.T() " << R - eVec.T() << endl;
        // }
    }
    else
    {
        WarningInFunction
            << "Non-unique eigenvectors, cannot compute transformation "
            << "from Cartesian axes" << endl;

        showTransform = false;
    }

    // calculate the total surface area

    scalar surfaceArea = 0;

    forAll(surf, facei)
    {
        const labelledTri& f = surf[facei];

        if (f[0] == f[1] || f[0] == f[2] || f[1] == f[2])
        {
            WarningInFunction
               << "Illegal triangle " << facei << " vertices " << f
               << " coords " << f.points(surf.points()) << endl;
        }
        else
        {
            surfaceArea += triPointRef
            (
                surf.points()[f[0]],
                surf.points()[f[1]],
                surf.points()[f[2]]
            ).mag();
        }
    }

    Info<< nl << setprecision(12)
        << "Density: " << density << nl
        << "Mass: " << m << nl
        << "Centre of mass: " << cM << nl
        << "Surface area: " << surfaceArea << nl
        << "Inertia tensor around centre of mass: " << nl << J << nl
        << "eigenValues (principal moments): " << eVal << nl
        << "eigenVectors (principal axes): " << nl
        << eVec.x() << nl << eVec.y() << nl << eVec.z() << endl;

    if (showTransform)
    {
        Info<< "Transform tensor from reference state (orientation):" << nl
            << eVec.T() << nl
            << "Rotation tensor required to transform "
               "from the body reference frame to the global "
               "reference frame, i.e.:" << nl
            << "globalVector = orientation & bodyLocalVector"
            << endl;

        Info<< nl
            << "Entries for sixDoFRigidBodyDisplacement boundary condition:"
            << nl
            << "        mass            " << m << token::END_STATEMENT << nl
            << "        centreOfMass    " << cM << token::END_STATEMENT << nl
            << "        momentOfInertia " << eVal << token::END_STATEMENT << nl
            << "        orientation     " << eVec.T() << token::END_STATEMENT
            << endl;
    }

    if (calcAroundRefPt)
    {
        Info<< nl << "Inertia tensor relative to " << refPt << ": " << nl
            << momentOfInertia::applyParallelAxisTheorem(m, cM, J, refPt)
            << endl;
    }

    OFstream str("axes.obj");

    Info<< nl << "Writing scaled principal axes at centre of mass of "
        << surfFileName <<  " to " << str.name() << endl;

    scalar scale = mag(cM - surf.points()[0])/eVal.component(findMin(eVal));

    meshTools::writeOBJ(str, cM);
    meshTools::writeOBJ(str, cM + scale*eVal.x()*eVec.x());
    meshTools::writeOBJ(str, cM + scale*eVal.y()*eVec.y());
    meshTools::writeOBJ(str, cM + scale*eVal.z()*eVec.z());

    for (label i = 1; i < 4; i++)
    {
        str << "l " << 1 << ' ' << i + 1 << endl;
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
