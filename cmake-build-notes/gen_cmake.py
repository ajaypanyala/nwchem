#!/usr/bin/python

import os
import glob

nwc_top = os.getcwd()
nwc_src = nwc_top+"/src"

#print "TOP=" + nwc_src

all_inc = "set(NWC_LIBRARIES ${NWCHEM_INSTALL_DIR}/external/ga/lib/libga.a)\n"
all_inc += "set(NWC_INCLUDES\n${NWC_MPI_INCLUDE}\n${GA_ROOT_DIR}/include\n"

for subdirs,_,_ in os.walk(nwc_src):
    if os.path.basename(subdirs) != 'tools':
        all_inc += subdirs + "\n"
all_inc += ")"
print all_inc

all_sd = []
globaldirs = [x for x in os.listdir(nwc_src) if os.path.isdir(os.path.join(nwc_src, x))]

# develop nwxc scfaux rimp2_grad python argos diana uccsdt rism geninterfac\
# e smd nbo leps lucia
exclude_mods = ['java','config','include','nbo','tools','nwxc','rism','python','smd','uccsdt','lucia','argos','rimp2_grad']

for subdirs in globaldirs:
    if os.getcwd() == nwc_src: continue
    if subdirs in exclude_mods: continue
 #   print subdirs

    os.chdir(nwc_src)
    os.chdir(subdirs)
    if not os.path.exists("GNUmakefile"): continue

    all_sd.append(subdirs)
    d = os.path.basename(subdirs)

    all_src = []
    #all_src = glob.glob("*.F")
    #all_src.extend(glob.glob("*.f"))
    for sdd,_,sdfiles in os.walk("."):
        if 'test' in sdd or 'detci_sym' in sdd: continue
        for sdf in sdfiles:
            _,ext = os.path.splitext(sdf)
            if ext==".F" or ext==".f" or ext==".c": all_src.append(sdd+"/"+sdf)

    cmakestr = ""
    if all_src:
        cmakestr = "\nset(" + d.upper() + "_SOURCE_FILES\n"
        for ff in all_src:
            cmakestr += ff +"\n"
        cmakestr += ")\n\n"
        cmakestr += "add_library(" + d + " STATIC ${" + d.upper() + "_SOURCE_FILES} ${NWC_LIBRARIES})\n"
        cmakestr += "target_include_directories(" + d + " PUBLIC ${NWC_INCLUDES})\n\n"


    with open("CMakeLists.txt", 'w') as cmf:
        cmf.write(cmakestr)
        #os.system("git add CMakeLists.txt")

os.chdir(nwc_src)
with open("CMakeLists.txt", 'w') as cmf:
    for modules in all_sd:
        cmf.write("add_subdirectory(" + modules + ")\n")

# nwc_src = nwc_top+"/tce"
# exclude_mods = ['java','config','include','nbo','tools','nwxc','leps','rism','python','smd','uccsdt','lucia','argos','rimp2_grad']
# for subdirs,_,_ in os.walk(nwc_src):
#     #print subdirs
#     os.chdir(nwc_src)
#     os.chdir(subdirs)
#     if not os.path.exists("GNUmakefile"): continue
#
#     d = os.path.basename(subdirs)
#     if d in exclude_mods: continue
#     #print d
#
#
#     libname = subdirs.replace("/","_")[1:]
#     print libname
#     cmakestr = ""
#
#     all_src = glob.glob("*.F")
#     all_src.extend(glob.glob("*.f"))
#     all_src.extend(glob.glob("*.c"))
#     sdd = [x for x in os.listdir(subdirs) if os.path.isdir(os.path.join(subdirs, x)) and x not in exclude_mods]
#
#     for x in sdd:
#         if os.path.exists(x+"/GNUmakefile"):
#             cmakestr += "add_subdirectory(" + x + ")\n"
#
#     if all_src:
#         cmakestr += "\nset(" + d.upper() + "_SOURCE_FILES\n"
#         for ff in all_src:
#             cmakestr += ff +"\n"
#         cmakestr += ")\n\n"
#         cmakestr += "add_library(" + libname + " STATIC ${" + d.upper() + "_SOURCE_FILES} ${NWC_LIBRARIES})\n"
#         cmakestr += "target_include_directories(" + libname + " PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})\n\n"
#
#
#     with open("CMakeLists.txt", 'w') as cmf:
#         cmf.write(cmakestr)
