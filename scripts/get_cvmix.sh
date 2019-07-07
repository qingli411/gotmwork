#!/bin/bash
# get CVMix from Github
#
# Qing Li, 20190706

# CVMix tag for build
# v0.94b-beta or later is required for Langmuir options.
repo_tag=v0.95-beta

# address for acquiring CVMix source code
repo_git_ssh_address=git@github.com:CVMix/CVMix-src.git
repo_git_http_address=https://github.com/CVMix/CVMix-src.git

# gotmwork environment file
gotmwork_env_file="${HOME}/.gotmwork_env.sh"
if [[ -f ${gotmwork_env_file} ]]; then
    source ${gotmwork_env_file}
else
    echo "** GOTMWORK environment not set. Use set_gotmwork_env.sh to set it up."
    exit 1
fi

# root directory
repo_root=${CVMIX_ROOT}

# check if git exist
git_exe=`which git`
if [[ "${git_exe}" = "" ]]; then
    echo -e "** Git not found. Stop.\n"
    exit 1
fi

# repo exists, check to see if it is the correct version.
# otherwise, clean up the directory to ensure it's updated.
if [[ -d ${repo_root} ]]; then
	if [[ -d ${repo_root}/.git ]]; then
		cd ${repo_root}
		CURR_TAG=$(git describe --tags)
		cd - &> /dev/null
		if [[ "${CURR_TAG}" == "${repo_tag}" ]]; then
			echo -e "** CVMix is updated. Skip update.\n"
		else
			rm -rf ${repo_root}
		fi
	else
		rm -rf ${repo_root}
	fi
fi

# repo does not exist, need to acquire source code
if [[ ! -d ${repo_root} ]]; then
    echo -e "** Acquiring CVMix source code...\n"
    git clone ${repo_git_ssh_address} ${repo_root}
    if [[ -d ${repo_root} ]]; then
        cd ${repo_root}
        git checkout ${repo_tag}
        cd - &> /dev/null
    else
        git clone ${repo_git_http_address} ${repo_root}
        if [[ -d ${repo_root} ]]; then
            cd ${repo_root}
            git checkout ${repo_tag}
            cd - &> /dev/null
        fi
    fi
fi
