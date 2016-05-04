BEGIN {
    FS = "|";
    data = "";
    flag = 0;
    firstRow = 0;
}


/\\table/ {
    flag = 1;
    firstRow = 1;
    next;
}


/\\endtable/ {
    if (data != "")
    {
        printf "<table class=\"OFTable\">\n";
        printf data;
        printf "</table>\n";
    }

    data = "";
    flag = 0;
    next;
}


/\\vartable/ {
    flag = 2;
    firstRow = 1;
    next;
}


/\\endvartable/ {
    if (data != "")
    {
        printf "<table border="0">\n";
        printf data;
        printf "</table>\n";
    }

    data = "";
    flag = 0;
    next;
}


/\\plaintable/ {
    flag = 3;
    firstRow = 1;
    next;
}


/\\endplaintable/ {
    if (data != "")
    {
        printf "<table border="0">\n";
        printf data;
        printf "</table>\n";
    }

    data = "";
    flag = 0;
    next;
}

{
    if (flag > 0)
    {
        data = (data "<tr>\n");
        if (flag == 1)
        {
            for (i = 0; i <= NF; i++)
            {
                if ((i != 0) && (firstRow == 1))
                {
                    data = (data "    <th align=\"center\"><b>"$i"</b></th>\n");
                }
                else
                {
                    if (i == 1)
                    {
                        data = (data "    <td>\\c "$i"</td>\n");
                    }
                    else if (i > 1)
                    {
                        data = (data "    <td>"$i"</td>\n");
                    }
                }
            }
        }
        else if (flag == 2)
        {
            for (i = 0; i <= NF; i++)
            {
                if (i == 1)
                {
                    data = (data "    <td style=\"padding-left: 10px\">\\f$"$i"\\f$</td>\n");
                    data = (data "    <td style=\"padding-left: 10px; padding-right: 10px;\">=</td>\n");
                }
                else if (i > 1)
                {
                    data = (data "    <td>"$i"</td>\n");
                }
            }
        }
        else if (flag == 3)
        {
            for (i = 0; i <= NF; i++)
            {
                if (i == 1)
                {
                    data = (data "    <td style=\"padding-left: 10px\">"$i"</td>\n");
                    data = (data "    <td style=\"padding-left: 10px; padding-right: 10px;\">:</td>\n");
                }
                else if (i > 1)
                {
                    data = (data "    <td>"$i"</td>\n");
                }
            }
        }
        data = (data "</tr>\n");
        firstRow = 0;
    }
    else
    {
        print $0
    }
}
