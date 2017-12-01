#!/usr/bin/python


from osgeo import gdal
import os
import sys
from collections import namedtuple
import json
from subprocess import call



AUTO_CONVERT_COG = False
AUTO_BUILD_OVERVIEWS = False
CATALOG_FILE = "/opt/eo_webtools/sentinel_webview/www/catalog.json"
S2VRT_PATH = "/usr/bin/s2vrt"
VRT_PATH = "/opt/VRT"


def process_file(filename):
    if not (filename.endswith(".zip") or filename.endswith(".ZIP")):
        print("WARNING: ignoring '" + filename + "'")
    print("WARNING: Sentinel 2 zip archives are currently not supported, ignoring '" + filename + "'")
    return None
    # if old format ...
    # elif new compact format ...




def process_dir(dirname):
    tmpname = dirname
    if tmpname.endswith(".SAFE"):
        tmpname = tmpname[:-5]
    tmpname = os.path.basename(tmpname);
    parts = tmpname.split("_")

    if not (parts[0] == "S2A" or parts[0] == "S2B"):

        return None


    if len(parts[1]) == 4:
        # old format
        print("WARNING: Sentinel 2 old format is currently not supported, ignoring '" + tmpname + "'")
        return None
    elif len(parts[1]) == 6:
        if parts[1] == "MSIL1C":
            # level 1c
            id = tmpname
            date = parts[2]
            path = os.path.join(dirname,"MTD_MSIL1C.xml")
            footprint = ""
            return {"id" : id, "path" : path, "date" : date, "footprint" : ""}

        elif parts[1] == "MSIL2A":
            # level 2a
            print("WARNING: invalid scene identifier, ignoring '" + tmpname + "'")
            return None
        else:
            return None
        # new format
    else:
        print("WARNING: invalid scene identifier, ignoring '" + tmpname + "'")
        return None

    return None







gdal.AllRegister()
print ("Using GDAL " + gdal.VersionInfo("RELEASE_NAME"))


assert(len(sys.argv) == 2)
root_dir = sys.argv[1]
assert(os.path.exists(root_dir) and os.path.isdir(root_dir))


paths = os.listdir(root_dir)


if not os.path.exists(CATALOG_FILE):
    f = open(CATALOG_FILE, "w")
    f.write("[]")
    f.close()


catalog = json.load(open(CATALOG_FILE))
print(catalog)

# remove missing scenes
missing_idx = []
for i in reversed(range(0, len(catalog))):
    if not os.path.exists(catalog[i]["path"]):
        missing_idx.append(i)

for i in range(0, len(missing_idx)):
    del catalog[missing_idx[i]]




for i in range(0, len(paths)):
    s = None
    if os.path.isfile(os.path.join(root_dir,paths[i])):
        process_file(os.path.join(root_dir,paths[i]))
    else:
       s  = process_dir(os.path.join(root_dir,paths[i]))

    if s is None:
        continue

    # check whether scene with id is already in catalog
    in_catalog = False
    for j in range(0, len(catalog)):
        if catalog[j]["id"] == s["id"]:
            in_catalog = True
            break


    if not in_catalog:
        catalog.append(s)


# iterate over all datasets and find non VRT scenes
for i in range(0, len(catalog)):
    if catalog[i]["path"].endswith(".VRT") or catalog[i]["path"].endswith(".vrt"):
        continue
    s2vrtcall = [S2VRT_PATH]
    if AUTO_CONVERT_COG:
        s2vrtcall.append("--cog")
    if AUTO_BUILD_OVERVIEWS:
        s2vrtcall.append("-o")
        s2vrtcall.append("2 4 8 16 32")
    s2vrtcall.append(catalog[i]["path"])
    outpath = os.path.join(VRT_PATH, catalog[i]["id"] + ".vrt")
    s2vrtcall.append(outpath)

    print("Running s2vrt for scene '" + catalog[i]["id"] + "' ...")

    #print (" ".join(s2vrtcall))
    if (call(s2vrtcall) == 0):
        catalog[i]["path"] = outpath
    open(CATALOG_FILE, "w").write(json.dumps(catalog, indent=4)) # TODO: append only current item?



open(CATALOG_FILE,"w").write(json.dumps(catalog, indent=4))






