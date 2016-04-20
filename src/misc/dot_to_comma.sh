#!/bin/bash

### script for converting commas (german decimal) to dots

new_dir="./converted"

mkdir "$new_dir"

for FILE in "$@"
do
    echo "processing $FILE"
    cat $FILE | sed 's/\./\,/g' > $new_dir/$FILE
done
