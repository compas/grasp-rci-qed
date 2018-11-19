#!/bin/bash
if ! [ "${TRAVIS}" == "true" ]; then
	>&2 echo "ERROR: Not running on Travis (TRAVIS='${TRAVIS}')."
	exit 1
fi

# cd into to root directory of the repository
cd "${TRAVIS_BUILD_DIR}" || exit

# Run Doxygen, which should create a doc/html directory
doxygen || {
	>&2 echo "ERROR: Doxygen failed with $?"
	exit 1
}
if ! [ -d "doc/html" ]; then
	>&2 echo "ERROR: doc/html missing."
	ls -Alh doc/ # for debugging
	exit 1
fi

# If this build is a non-PR build of the master branch, try to deploy to gh-pages
if [ "${TRAVIS_BRANCH}" == "master" ] && [ "${TRAVIS_PULL_REQUEST}" == "false" ]; then
	# Write the SSH deploy key to a file.
	if [ -z ${GITHUB_DEPLOY+x} ]; then >&2 echo "ERROR: \$GITHUB_DEPLOY variable is unset."; exit 1; fi
	base64 -d > .ci/github_deploy <<-EOF
		${GITHUB_DEPLOY}
	EOF
	chmod og-rwx .ci/github_deploy # the SSH key needs conservative permissions
	ssh-add .ci/github_deploy # make the key available to ssh-agent

	# Pushing the documentation to the gh-pages branch. First, we'll create a
	# temporary directory and clone the repository again, and cd there:
	gh_pages_dir=`mktemp -d`
	echo "INFO: Cloning gh-pages into ${gh_pages_dir}"
	git clone -b gh-pages git@github.com:mortenpi/grasp-rci-qed.git "${gh_pages_dir}" || exit
	cd "${gh_pages_dir}" || exit
	# Remove all the old files (but not stuff under .git)
	find . -type f -not -regex "\./\.git/.*" -delete
	# Copy over the new docs:
	rsync -a ${TRAVIS_BUILD_DIR}/doc/html/ . || exit
	# Git add, commit and push the updated documentation:
	git add -A || exit
	git commit -m "Deploy documentation" || {
		echo "git commit returned a non-zero exit code ($?). This is expected if there are no changes to commit."
		exit 0
	}
	git push --set-upstream origin gh-pages || exit

	echo "Documentation deployment complete."
else
	echo "Skipping deployment."
	echo "  TRAVIS_BRANCH=${TRAVIS_BRANCH}"
	echo "  TRAVIS_PULL_REQUEST=${TRAVIS_PULL_REQUEST}"
fi
