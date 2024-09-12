word meshPath;
if (args.optionFound("mesh"))
{
    meshPath = "meshes"/args.optionRead<word>("mesh");
}
