#!/bin/bash
echo $(git log | head -n1 | cut -d' ' -f2)
