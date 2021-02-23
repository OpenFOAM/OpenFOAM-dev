fv::options& fvOptions(fv::options::New(mesh));

if (!fvOptions.optionList::size())
{
    Info << "No finite volume options present" << endl;
}
