#!/usr/bin/env bash

# make sure this is being run from the course_files directory
if [ ! -f "download_course_data.sh" ]; then
    echo "Please run this script from the course_files directory"
    exit 1
fi

# Download the data for the course
wget -O data.zip "https://www.dropbox.com/scl/fo/9x1rg6qxqw5crq2vtb1ho/AMawguf1kqRYQQs-qZPhFZA?rlkey=6y4w1skyjzpq36zfyocis24t6&st=uwfhd31u&dl=1"
unzip -o data.zip
rm data.zip