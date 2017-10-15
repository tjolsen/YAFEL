#!/bin/bash
gitstr=$(echo $(git log | head -n1 | cut -d' ' -f2))
echo \"$gitstr\"
