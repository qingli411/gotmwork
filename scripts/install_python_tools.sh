#!/bin/bash
#
# This script installs python tools for GOTM
#
# Qing Li, 20171024

function if_python_module() {
    # Return 1 if the module (name as argument) is installed in the current
    # Python environment, 0 otherwise
    python -c "import pkgutil; print(1 if pkgutil.find_loader(\"$1\") else 0)"
}

# root directory
gotm_root="${HOME}/models/gotm"
mkdir -p ${gotm_root}
cd ${gotm_root}

# github repositories
gitrep_editscenario="https://github.com/BoldingBruggeman/editscenario.git"
gitrep_xmlstore="https://github.com/BoldingBruggeman/xmlstore.git"
gitrep_xmlplot="https://github.com/BoldingBruggeman/xmlplot.git"
gitrep_gotmgui="https://github.com/BoldingBruggeman/gotmgui.git"

# clone tools
tools_root="${gotm_root}/tools"
mkdir -p ${tools_root}
cd ${tools_root}
for f in editscenario xmlstore xmlplot gotmgui
do
    # check if installed
    if [[ $(if_python_module ${f}) == 1 ]]; then
        echo "${f} installed. Skip."
    else
        # clone source code from github
        mkdir -p ${f}
        gvar="gitrep_${f}"
        gitrep=$(eval echo \$${gvar})
        echo "Cloning source code from ${gitrep} ..."
        git clone ${gitrep} ${f}
        # install tool
        cd ${f}
        pip wheel .
        pip install ${f}-*.whl
        cd ${tools_root}
    fi
done
