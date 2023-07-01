import os
HOME=os.environ["HOME"]
config = {
    "all":{
        "OBABEL_EXE": f"{HOME}/mambaforge/envs/basic/bin/obabel",
        "PLOP_EXE": f"{HOME}/opt/plop/25.1/plop32-static",
        "plop_data": f"{HOME}/opt/plop/25.1/data",
        "DOCKBASE": f"{HOME}/opt/ucsfdock/3.7/bela30c",
        "SCHUTILS": f"{HOME}/soft/schrodinger2021-2/utilities"
    }
}