#!/bin/bash

sort -k1,1 -k3,3n -k4,4n -k8,8n -k9,9n "$1" > "$2"
