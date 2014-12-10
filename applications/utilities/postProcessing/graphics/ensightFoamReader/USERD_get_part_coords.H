// Note: coord_array is 1-based.

int USERD_get_part_coords
(
    int part_number,
    float **coord_array
)
{

#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_part_coords" << endl <<
        "part_number = " << part_number << endl << flush;
#endif

    if (part_number == 1)
    {

        //#       include "checkForNewMesh.H"

        const vectorField& points = meshPtr->points();
        label nPoints = points.size();

        for (label indx=0; indx<nPoints; indx++)
        {
            coord_array[0][indx+1] = (float)points[indx].x();
            coord_array[1][indx+1] = (float)points[indx].y();
            coord_array[2][indx+1] = (float)points[indx].z();
        }
    }
    else if (part_number < nPatches+2)
    {

        //#       include "checkForNewMesh.H"

        label patchi = part_number-2;
        const polyBoundaryMesh& bMesh = meshPtr->boundaryMesh();
        const vectorField& points = bMesh[patchi].points();
        label nPoints = points.size();

        for (label indx=0; indx<nPoints; indx++)
        {
            coord_array[0][indx+1] = (float)points[indx].x();
            coord_array[1][indx+1] = (float)points[indx].y();
            coord_array[2][indx+1] = (float)points[indx].z();
        }

    }
    else if (part_number == nPatches+2)
    {
        label indx = 1;

        forAllConstIter(Cloud<passiveParticle>, *sprayPtr, iter)
        {
            coord_array[0][indx] = (float)iter().position().x();
            coord_array[1][indx] = (float)iter().position().y();
            coord_array[2][indx] = (float)iter().position().z();
            indx++;
        }

    }
    else
    {
        return Z_ERR;
    }

#ifdef ENSIGHTDEBUG
    Info<< "Leaving: USERD_get_part_coords" << endl << flush;
#endif

    return Z_OK;
}
