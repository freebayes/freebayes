#!/bin/sh

BASE_DIR=$1

RELEASED_VERSION_FILE=${BASE_DIR}/src/version_release.txt
if [ -d "${BASE_DIR}/.git" ] && which git > /dev/null
then
	DETECTED_VERSION=$(git describe --always --tags --dirty)
else
	DETECTED_VERSION=$(grep -v "^#" $RELEASED_VERSION_FILE)
fi

VERSION_FILE=${BASE_DIR}/src/version_git.h
CURRENT_VERSION=""
[ -e "$VERSION_FILE" ] && CURRENT_VERSION=$(grep "define VERSION_GIT " "$VERSION_FILE" | cut -f3 -d" " | sed 's/"//g')
echo "DETECTED_VERSION = $DETECTED_VERSION"
echo "CURRENT_VERSION  = $CURRENT_VERSION"
if [ "${DETECTED_VERSION}" != "${CURRENT_VERSION}" ]
then
	echo "Updating version file."
	echo "#ifndef VERSION_GIT_H" > $VERSION_FILE
	echo "#define VERSION_GIT_H" >> $VERSION_FILE
	echo "#define VERSION_GIT \"${DETECTED_VERSION}\"" >> $VERSION_FILE
	echo "#endif /* VERSION_GIT_H */" >> $VERSION_FILE
fi
